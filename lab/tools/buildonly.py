"""Build only (no emission): wires and the flat CCX count from the hierarchical build, under flags."""
import ast, importlib, sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
for kv in sys.argv[2:]:
    k, v = kv.split("=", 1)
    if k == "two_bit.compare_blocks.max_blocks":
        import qec.two_bit as tb; tb.compare_blocks.__defaults__ = (int(v),); continue
    mod, flag = k.rsplit(".", 1)
    m = importlib.import_module("qec." + mod); assert hasattr(m, flag), k; setattr(m, flag, ast.literal_eval(v))
from qec.curve import SECP256K1, SECP256K1_G
from qec.point_add import ControlledPointAdd
from qec.builder import CCX
gen = ControlledPointAdd(SECP256K1, SECP256K1_G)
C = gen.build(True)
memo = {}
def ccx(block):
    if id(block) not in memo:
        memo[id(block)] = block.ops[0::4].count(CCX) + sum(ccx(c.block) for c in block.calls)
    return memo[id(block)]
top = C.prims[0::4].count(CCX) + sum(ccx(c.block) for c in C.calls)
print(sys.argv[2:], "wires", C.num_wires, "ccx", top, "T", 7 * top)
