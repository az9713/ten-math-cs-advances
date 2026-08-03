"""Brute-force edge-expansion data for the ch03 explainer (Fig: expander
vs bottleneck).  Two 16-vertex 3-regular multigraphs:

  * dumbbell: two 8-vertex 'blobs' (cycle + chords) joined by ONE edge
    (the two endpoints drop one internal edge to stay 3-regular);
  * circulant C16(1, 7): vertices 0..15, i ~ i+1 and i ~ i+7 (mod 16)
    -- i+7 pairing i<->i+7... 7*2=14 != 0 mod 16 so (1,7) gives degree 4;
    use C16(2, 3)?  Simpler: Moebius-Kantor graph, the generalized
    Petersen graph GP(8,3), 16 vertices, 3-regular, a good expander.

Prints h(G) = min_{0<|U|<=n/2} |boundary(U)|/|U| for both, plus the
worst set of the dumbbell.  Run: python figgen/ch03_expander_data.py
"""
from itertools import combinations


def h_const(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    best = None
    best_set = None
    for k in range(1, n//2 + 1):
        for U in combinations(range(n), k):
            Uset = set(U)
            bnd = sum(1 for a, b in edges
                      if (a in Uset) != (b in Uset))
            r = bnd/len(U)
            if best is None or r < best:
                best, best_set = r, U
    return best, best_set


# Moebius-Kantor = generalized Petersen GP(8,3):
# outer 8-cycle 0..7, inner vertices 8..15, spokes i - (i+8),
# inner edges (8+i) - (8+(i+3) mod 8)
mk = [(i, (i+1) % 8) for i in range(8)]
mk += [(i, i+8) for i in range(8)]
mk += [(8+i, 8+((i+3) % 8)) for i in range(8)]

# dumbbell: two blobs {0..7}, {8..15}; each blob = 8-cycle + 3 chords
# to make everyone degree 3 except the two bridge endpoints (0 and 8),
# which have degree 2 inside and get the bridge.
blob = [(i, (i+1) % 8) for i in range(8)] + [(1, 5), (2, 6), (3, 7)]
db = list(blob)
db += [(8+a, 8+b) for a, b in blob]
db += [(0, 8)]
deg = [0]*16
for a, b in db:
    deg[a] += 1
    deg[b] += 1
print('dumbbell degrees:', sorted(set(deg)))

h_mk, _ = h_const(16, mk)
h_db, worst = h_const(16, db)
print('Moebius-Kantor h =', h_mk)
print('dumbbell h =', h_db, 'worst set =', worst)
