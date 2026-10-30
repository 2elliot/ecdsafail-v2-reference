"""Reference short-Weierstrass arithmetic, independent of the circuit code."""


class Curve:
    def __init__(self, q, a, b, name=None):
        self.q, self.a, self.b, self.name = q, a, b, name or f"y^2=x^3+{a}x+{b} mod {q}"

    def on_curve(self, P):
        if P is None:
            return True
        x, y = P
        return 0 <= x < self.q and 0 <= y < self.q and (y * y - x * x * x - self.a * x - self.b) % self.q == 0

    def neg(self, P):
        return None if P is None else (P[0], -P[1] % self.q)

    def add(self, P, Q):
        q = self.q
        if P is None:
            return Q
        if Q is None:
            return P
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if (y1 + y2) % q == 0:
                return None
            lam = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, q) % q
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, q) % q
        x3 = (lam * lam - x1 - x2) % q
        return x3, (lam * (x1 - x3) - y1) % q

    def mul(self, k, P):
        R = None
        while k:
            if k & 1:
                R = self.add(R, P)
            P = self.add(P, P)
            k >>= 1
        return R

    def points(self):
        q = self.q
        roots = {}
        for y in range(q):
            roots.setdefault(y * y % q, []).append(y)
        pts = [None]
        for x in range(q):
            for y in roots.get((x ** 3 + self.a * x + self.b) % q, []):
                pts.append((x, y))
        return pts


SECP256K1 = Curve(2**256 - 2**32 - 977, 0, 7, "secp256k1")
SECP256K1_G = (
    0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
    0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8,
)
SECP256K1_R = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141


def is_prime(m):
    if m < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if m % p == 0:
            return m == p
    d, s = m - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, m)
        if x in (1, m - 1):
            continue
        for _ in range(s - 1):
            x = x * x % m
            if x == m - 1:
                break
        else:
            return False
    return True


def toy_secp_curves(bits):
    """Primes q of the given bit size where y^2 = x^3 + 7 has prime order."""
    out = []
    for q in range(1 << (bits - 1) | 1, 1 << bits, 2):
        if is_prime(q):
            c = Curve(q, 0, 7)
            if is_prime(len(c.points())):
                out.append(c)
    return out
