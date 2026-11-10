"""Board tick model on a builder op array (X=0, CX=1, CCX=2, SWAP=3), all wires free at 0."""
X, CX, CCX, SWAP = 0, 1, 2, 3


def ticks(ops):
    r = {}
    ccx = 0
    for i in range(0, len(ops), 4):
        k, a, b, c = ops[i:i + 4]
        if k == X:
            r[a] = r.get(a, 0) + 1
        elif k == CX:
            t = max(r.get(a, 0), r.get(b, 0)) + 1
            r[a] = r[b] = t
        elif k == SWAP:
            t = max(r.get(a, 0), r.get(b, 0)) + 3
            r[a] = r[b] = t
        else:
            ccx += 1
            ra, rb, rt = r.get(a, 0), r.get(b, 0), r.get(c, 0)
            rt += 1
            x = max(rb, rt) + 1; rb = rt = x
            rt += 1
            x = max(ra, rt) + 1; ra = rt = x
            rt += 1
            x = max(rb, rt) + 1; rb = rt = x
            rt += 1
            x = max(ra, rt) + 1; ra = rt = x
            rb += 1; rt += 1
            rt += 1
            x = max(ra, rb) + 1; ra = rb = x
            ra += 1; rb += 1
            x = max(ra, rb) + 1; ra = rb = x
            r[a], r[b], r[c] = ra, rb, rt
    return ccx, max(r.values(), default=0)


