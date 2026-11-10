"""Puts a circuit file into the board's upload format.

    python3 wrap.py my-circuit.json mine.method.json
"""
import json
import sys

G = {"x": "0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798",
     "y": "0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8"}

circuit_path, method_path = sys.argv[1:3]
method = {"method_id": "method_shor_secp256k1_canonical", "public_key": G, "emitted_circuit_digest": "0" * 64,
          "circuit": json.load(open(circuit_path)), "universal_certificate": None, "wrapper_certificate": None, "schedule": None,
          "observations": None, "retry_model": None, "statistical_binding": None, "claims_family": False, "claims_secp256k1": True}
json.dump(method, open(method_path, "w"), separators=(",", ":"))
print("wrote", method_path)
