//! Feeds an emitted `ArithmeticCircuit` through ecdsafail-v2's own code paths and
//! prints one JSON receipt. The circuit is read from a path, or stdin for `-`.

use std::io::{BufReader, Read};

use num_bigint::BigUint;
use serde::Deserialize;
use serde_json::json;
use v2_contracts::AffinePoint;
use v2_ir::{ArithmeticCircuit, CircuitValidator, IngestionLimits, RegisterRole};
use v2_qualification::{controlled_exact_profile, EngineContext, QualificationBundle, QualificationGraph, SubmittedMethod};
use v2_runtime::{BitSlicedBlock, BitSlicedExecutor};

#[derive(Deserialize)]
struct Lane {
    c: u8,
    x: String,
    y: String,
    ex: String,
    ey: String,
}

#[derive(Deserialize)]
struct Cases {
    public_key: [String; 2],
    lanes: Vec<Lane>,
}

fn hex(s: &str) -> BigUint {
    BigUint::parse_bytes(s.trim_start_matches("0x").as_bytes(), 16).expect("hex")
}

fn arg(name: &str) -> Option<String> {
    let args: Vec<String> = std::env::args().collect();
    args.iter().position(|a| a == name).map(|i| args[i + 1].clone())
}

fn role_wires(c: &ArithmeticCircuit, role: RegisterRole) -> Vec<u32> {
    c.registers.iter().filter(|r| r.role == role).flat_map(|r| r.wires.clone()).collect()
}

/// Runs every lane on the repository's bit-sliced executor, 64 lanes per block.
fn execute(c: &ArithmeticCircuit, cases: &Cases) -> serde_json::Value {
    let ctrl = role_wires(c, RegisterRole::Control);
    let acc = role_wires(c, RegisterRole::Accumulator);
    let n = acc.len() / 2;
    let live: std::collections::HashSet<u32> = ctrl.iter().chain(acc.iter()).copied().collect();
    let (mut failures, mut phase_lanes, mut dirty_lanes) = (Vec::new(), 0u32, 0u32);
    for (b, chunk) in cases.lanes.chunks(64).enumerate() {
        let mask = if chunk.len() == 64 { u64::MAX } else { (1u64 << chunk.len()) - 1 };
        let mut blk = BitSlicedBlock::new(c.num_wires as usize, mask);
        for (j, lane) in chunk.iter().enumerate() {
            let bit = 1u64 << j;
            if lane.c == 1 {
                blk.wires[ctrl[0] as usize] |= bit;
            }
            let (x, y) = (hex(&lane.x), hex(&lane.y));
            for i in 0..n {
                if x.bit(i as u64) {
                    blk.wires[acc[i] as usize] |= bit;
                }
                if y.bit(i as u64) {
                    blk.wires[acc[n + i] as usize] |= bit;
                }
            }
        }
        BitSlicedExecutor::execute_block(c, &mut blk);
        let mut dirty = 0u64;
        for w in 0..c.num_wires {
            if !live.contains(&w) {
                dirty |= blk.wires[w as usize];
            }
        }
        dirty_lanes += (dirty & mask).count_ones();
        phase_lanes += (blk.phase_signs & mask).count_ones();
        for (j, lane) in chunk.iter().enumerate() {
            let bit = 1u64 << j;
            let (ex, ey) = (hex(&lane.ex), hex(&lane.ey));
            let ok_ctrl = ((blk.wires[ctrl[0] as usize] & bit) != 0) == (lane.c == 1);
            let ok_acc = (0..n).all(|i| {
                ((blk.wires[acc[i] as usize] & bit) != 0) == ex.bit(i as u64)
                    && ((blk.wires[acc[n + i] as usize] & bit) != 0) == ey.bit(i as u64)
            });
            if !(ok_ctrl && ok_acc) || dirty & bit != 0 || blk.phase_signs & bit != 0 {
                failures.push(b * 64 + j);
            }
        }
    }
    json!({
        "executor": "v2_runtime::BitSlicedExecutor::execute_block",
        "lanes": cases.lanes.len(),
        "failures": failures.len(),
        "first_failures": &failures[..failures.len().min(5)],
        "lanes_with_dirty_scratch": dirty_lanes,
        "lanes_with_phase": phase_lanes,
    })
}

fn main() {
    let path = arg("--circuit").expect("--circuit <path|->");
    let mut text = Vec::new();
    if path == "-" {
        std::io::stdin().lock().read_to_end(&mut text).unwrap();
    } else {
        BufReader::new(std::fs::File::open(&path).unwrap()).read_to_end(&mut text).unwrap();
    }
    let circuit_bytes = text.len();
    let circuit: ArithmeticCircuit = serde_json::from_slice(&text).expect("ArithmeticCircuit JSON");
    drop(text);

    let service = IngestionLimits {
        max_wires: v2_capacity::FrozenLimits::official_secp256k1().max_qubits,
        max_ops: v2_capacity::MAX_DECLARED_RECORDS,
        ..IngestionLimits::default()
    };
    let verdict = |l: &IngestionLimits| match CircuitValidator::validate_arithmetic(&circuit, l, false) {
        Ok(()) => "PASS".to_string(),
        Err(e) => format!("REJECT: {}", e),
    };
    let ccx = circuit.operations.iter().filter(|o| matches!(o, v2_ir::ArithmeticOp::Ccx { .. })).count();
    let mut receipt = json!({
        "circuit_json_bytes": circuit_bytes,
        "num_wires": circuit.num_wires,
        "operations": circuit.operations.len(),
        "ccx": ccx,
        "validator_service_caps": verdict(&service),
        "validator_default_ingestion_limits": verdict(&IngestionLimits::default()),
    });

    if let Some(cases_path) = arg("--cases") {
        let cases: Cases = serde_json::from_str(&std::fs::read_to_string(cases_path).unwrap()).unwrap();
        receipt["execution"] = execute(&circuit, &cases);
        let digest = v2_composition::compute_operations_hash(&circuit.operations);
        let method = SubmittedMethod {
            method_id: "qec-exact-dialog-controlled-add".into(),
            public_key: AffinePoint { x: cases.public_key[0].clone(), y: cases.public_key[1].clone() },
            emitted_circuit_digest: digest,
            circuit: Some(circuit.clone()),
            universal_certificate: None,
            wrapper_certificate: None,
            schedule: None,
            observations: None,
            retry_model: None,
            statistical_binding: None,
            claims_family: false,
            claims_secp256k1: true,
        };
        let bytes = serde_json::to_vec(&method).unwrap();
        let reparsed = serde_json::from_slice::<SubmittedMethod>(&bytes).map(|_| "PARSED").unwrap_or("REJECTED");
        let limit = v2_api_contracts::DEFAULT_MAX_UPLOAD_BYTES;
        receipt["submission"] = json!({
            "bytes": bytes.len(),
            "upload_limit_bytes": limit,
            "within_upload_limit": bytes.len() <= limit,
            "service_parse": reparsed,
            "emitted_circuit_digest": format!("{:?}", digest),
        });
        if let Some(out) = arg("--write-submission") {
            std::fs::write(out, &bytes).unwrap();
        }
        let decision = QualificationGraph::evaluate(
            &QualificationBundle::assemble(method, EngineContext::unrevoked(controlled_exact_profile())),
            None,
        );
        receipt["qualification"] = json!({
            "is_qualified": decision.is_qualified,
            "board_status": decision.board_status,
            "failure_reasons": decision.failure_reasons.iter().map(|e| e.to_string()).collect::<Vec<_>>(),
        });
    }

    if std::env::args().any(|a| a == "--lower") {
        let scratch = role_wires(&circuit, RegisterRole::Scratch);
        match v2_lowering::CircuitDecomposer::lower_circuit(circuit.num_wires, &circuit.operations, &[], &scratch) {
            Ok(l) => {
                receipt["lowering"] = json!({
                    "t": l.t_count, "tdg": l.tdg_count, "clifford": l.clifford_count,
                    "measurements": l.measurement_count, "t_plus_tdg_per_ccx": (l.t_count + l.tdg_count) as f64 / ccx.max(1) as f64,
                })
            }
            Err(e) => receipt["lowering"] = json!({ "error": e.to_string() }),
        }
    }
    println!("{}", serde_json::to_string_pretty(&receipt).unwrap());
}
