#!/usr/bin/env python3
"""
Test Corollary 2.1: Do eigenvector-threshold cuts satisfy the Hoffman condition?

For each graph with χ > ω:
1. Compute the λ_min-eigenvector
2. Sweep threshold τ over all distinct eigenvector values
3. For each cut (S_τ, T_τ): compute L_Hof for both parts
4. Report whether any τ satisfies Lemma 2

This is the critical test: if spectral cuts work empirically,
it motivates the proof of Lemma 3.
"""

import numpy as np
import networkx as nx
from itertools import combinations

PHI = (1 + np.sqrt(5)) / 2
PHI_INV = PHI - 1
PHI_SQ = PHI + 1
SQRT5 = np.sqrt(5)


# ─── Graph Library (same as explore script) ──────────────────────────────────

def build_test_graphs():
    graphs = {}
    graphs["C5"] = nx.cycle_graph(5)
    graphs["C7"] = nx.cycle_graph(7)
    graphs["C9"] = nx.cycle_graph(9)
    graphs["Petersen"] = nx.petersen_graph()
    graphs["Grötzsch"] = nx.mycielski_graph(4)
    graphs["Mycielski5"] = nx.mycielski_graph(5)
    graphs["C̄7"] = nx.complement(nx.cycle_graph(7))
    graphs["Icosahedron"] = nx.icosahedral_graph()
    graphs["Dodecahedron"] = nx.dodecahedral_graph()
    return graphs


# ─── Known chromatic numbers (greedy can be wrong) ──────────────────────────

KNOWN_CHI = {
    "C5": 3, "C7": 3, "C9": 3,
    "Petersen": 3, "Grötzsch": 4, "Mycielski5": 5,
    "C̄7": 4, "Icosahedron": 4, "Dodecahedron": 4,  # dodec is 3 actually...
}


# ─── Spectral tools ─────────────────────────────────────────────────────────

def get_eigendata(G):
    """Return sorted eigenvalues and eigenvectors of adjacency matrix."""
    A = nx.adjacency_matrix(G).toarray().astype(float)
    eigenvalues, eigenvectors = np.linalg.eigh(A)
    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    return eigenvalues[idx], eigenvectors[:, idx]


def hoffman_bound(G):
    """Compute L_Hof(G) = 1 - λ_max/λ_min. Returns 0 if λ_min >= 0."""
    A = nx.adjacency_matrix(G).toarray().astype(float)
    eigs = np.linalg.eigvalsh(A)
    lam_max = eigs[-1]  # numpy sorts ascending by default for eigvalsh
    lam_min = eigs[0]

    # Actually eigvalsh returns ascending, so:
    lam_min_val = np.min(eigs)
    lam_max_val = np.max(eigs)

    if lam_min_val >= -1e-10:
        return 1.0  # graph has no negative eigenvalue → bipartite or edgeless
    return 1 - lam_max_val / lam_min_val


def clique_number(G):
    if G.number_of_nodes() == 0:
        return 0
    return max(len(c) for c in nx.find_cliques(G))


def chromatic_number_greedy(G):
    if G.number_of_nodes() == 0:
        return 0
    if G.number_of_edges() == 0:
        return 1
    coloring = nx.coloring.greedy_color(G, strategy="largest_first")
    return max(coloring.values()) + 1


# ─── Spectral cut testing ───────────────────────────────────────────────────

def test_spectral_cuts(name, G, chi_G):
    """
    Test all eigenvector-threshold cuts for the λ_min eigenvector.
    Also test other eigenvectors (golden eigenvalues).
    """
    n = G.number_of_nodes()
    nodes = list(G.nodes())
    omega_G = clique_number(G)

    if chi_G <= omega_G:
        return None

    eigenvalues, eigenvectors = get_eigendata(G)
    lam_min = eigenvalues[-1]

    # All valid (s,t) pairs
    st_pairs = []
    for s in range(2, chi_G):
        t = chi_G - s + 1
        if t >= 2:
            st_pairs.append((s, t))

    print(f"\n{'='*70}")
    print(f"  {name}: n={n}, χ={chi_G}, ω={omega_G}, λ_min={lam_min:.4f}")
    print(f"  Valid (s,t) pairs: {st_pairs}")
    print(f"{'='*70}")

    # Test cuts from multiple eigenvectors
    eigenvectors_to_test = []

    # Always test λ_min eigenvector
    eigenvectors_to_test.append(("λ_min", eigenvalues[-1], eigenvectors[:, -1]))

    # Also test second-most-negative
    if n >= 3:
        eigenvectors_to_test.append(("λ_{n-1}", eigenvalues[-2], eigenvectors[:, -2]))

    # Test eigenvectors near golden values
    for i, lam in enumerate(eigenvalues):
        for gname, gval in [("φ", PHI), ("-φ", -PHI), ("φ⁻¹", PHI_INV),
                            ("-φ⁻¹", -PHI_INV), ("√5", SQRT5), ("-√5", -SQRT5),
                            ("φ²", PHI_SQ), ("-φ²", -PHI_SQ)]:
            if abs(lam - gval) < 0.02 and i not in [idx for _, _, v in eigenvectors_to_test
                                                       for idx in [np.argmin(np.abs(eigenvalues - _)) 
                                                                   for _ in [] ]]:
                eigenvectors_to_test.append((f"λ≈{gname}", lam, eigenvectors[:, i]))
                break

    # Deduplicate (keep unique eigenvectors)
    seen = set()
    unique_evecs = []
    for label, lam, vec in eigenvectors_to_test:
        key = tuple(np.round(np.abs(vec), 4))
        if key not in seen:
            seen.add(key)
            unique_evecs.append((label, lam, vec))

    results = {}

    for evec_label, evec_lam, x in unique_evecs:
        print(f"\n  --- Eigenvector: {evec_label} (λ={evec_lam:.4f}) ---")
        print(f"  Components: [{', '.join(f'{v:.3f}' for v in x[:min(12,n)])}{'...' if n>12 else ''}]")

        # Get unique thresholds (midpoints between sorted unique values)
        sorted_vals = np.sort(np.unique(np.round(x, 10)))
        thresholds = []
        for i in range(len(sorted_vals) - 1):
            thresholds.append((sorted_vals[i] + sorted_vals[i+1]) / 2)

        # Also test cutting at each unique value (vertex on the boundary goes to S)
        for v in sorted_vals:
            thresholds.append(v - 1e-12)

        thresholds = sorted(set(np.round(thresholds, 12)))

        for s, t in st_pairs:
            best_tau = None
            best_margins = None

            for tau in thresholds:
                S_nodes = [nodes[i] for i in range(n) if x[i] >= tau]
                T_nodes = [nodes[i] for i in range(n) if x[i] < tau]

                if len(S_nodes) < 1 or len(T_nodes) < 1:
                    continue

                GS = G.subgraph(S_nodes)
                GT = G.subgraph(T_nodes)

                # Skip if either part has no edges (Hoffman undefined / trivial)
                if GS.number_of_edges() == 0 or GT.number_of_edges() == 0:
                    continue

                L_S = hoffman_bound(GS)
                L_T = hoffman_bound(GT)

                if L_S > s - 1 and L_T > t - 1:
                    margin = min(L_S - (s-1), L_T - (t-1))
                    if best_tau is None or margin > best_margins:
                        best_tau = tau
                        best_margins = margin
                        best_LS = L_S
                        best_LT = L_T
                        best_S = S_nodes
                        best_T = T_nodes

            if best_tau is not None:
                print(f"  (s={s},t={t}): ✅ HOFFMAN CUT FOUND at τ={best_tau:.4f}")
                print(f"    |S|={len(best_S)}, |T|={len(best_T)}")
                print(f"    L_Hof(G[S])={best_LS:.4f} > {s-1} ✓")
                print(f"    L_Hof(G[T])={best_LT:.4f} > {t-1} ✓")
                print(f"    margin={best_margins:.4f}")

                # Verify with actual chromatic number
                chi_S = chromatic_number_greedy(G.subgraph(best_S))
                chi_T = chromatic_number_greedy(G.subgraph(best_T))
                tihany_valid = chi_S >= s and chi_T >= t
                print(f"    χ(G[S])={chi_S}, χ(G[T])={chi_T} → Tihany: {'✅' if tihany_valid else '❌'}")

                results[(name, evec_label, s, t)] = {
                    "found": True, "tau": best_tau,
                    "L_S": best_LS, "L_T": best_LT,
                    "chi_S": chi_S, "chi_T": chi_T,
                    "tihany_valid": tihany_valid,
                }
            else:
                print(f"  (s={s},t={t}): ❌ No Hoffman cut found via {evec_label}")
                results[(name, evec_label, s, t)] = {"found": False}

    # Also do brute-force check: does ANY partition satisfy Hoffman?
    if n <= 15:
        print(f"\n  --- Brute-force Hoffman partition search (n={n}) ---")
        for s, t in st_pairs:
            found_bf = False
            for size_S in range(1, n):
                if found_bf:
                    break
                for S_idx in combinations(range(n), size_S):
                    S_nodes = [nodes[i] for i in S_idx]
                    T_nodes = [nodes[i] for i in range(n) if i not in S_idx]

                    GS = G.subgraph(S_nodes)
                    GT = G.subgraph(T_nodes)

                    if GS.number_of_edges() == 0 or GT.number_of_edges() == 0:
                        continue

                    L_S = hoffman_bound(GS)
                    L_T = hoffman_bound(GT)

                    if L_S > s - 1 and L_T > t - 1:
                        found_bf = True
                        chi_S = chromatic_number_greedy(GS)
                        chi_T = chromatic_number_greedy(GT)
                        print(f"  (s={s},t={t}): ✅ Brute-force Hoffman partition found")
                        print(f"    S={S_nodes}, T={T_nodes}")
                        print(f"    L_Hof(G[S])={L_S:.4f}, L_Hof(G[T])={L_T:.4f}")
                        print(f"    χ(G[S])={chi_S}, χ(G[T])={chi_T}")
                        break

            if not found_bf:
                print(f"  (s={s},t={t}): ❌ No Hoffman partition exists (brute-force)")
                print(f"    → Hoffman bound too weak for this graph/pair")

    return results


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  Lemma 2 / Corollary 2.1: Spectral Cut Testing                    ║")
    print("║  Testing: does λ_min-eigenvector threshold give Tihany partitions? ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")

    graphs = build_test_graphs()

    # Override dodecahedron χ (greedy might give 4, actual is 3)
    KNOWN_CHI["Dodecahedron"] = 3  # χ(dodecahedron) = 3 (not 4)

    all_results = {}

    for name, G in graphs.items():
        chi = KNOWN_CHI.get(name, chromatic_number_greedy(G))
        omega = clique_number(G)
        if chi > omega:
            r = test_spectral_cuts(name, G, chi)
            if r:
                all_results.update(r)

    # Summary
    print(f"\n{'='*70}")
    print(f"  SUMMARY: Spectral Cut Results")
    print(f"{'='*70}")

    spectral_found = sum(1 for v in all_results.values() if v.get("found"))
    spectral_total = len(all_results)
    print(f"  Spectral cuts found: {spectral_found}/{spectral_total}")

    print(f"\n  Detailed results:")
    for key, val in sorted(all_results.items()):
        name, evec, s, t = key
        status = "✅" if val.get("found") else "❌"
        print(f"    {name:15s} {evec:10s} (s={s},t={t}): {status}", end="")
        if val.get("found"):
            print(f"  τ={val['tau']:.4f}  L_S={val['L_S']:.3f}  L_T={val['L_T']:.3f}", end="")
            if val.get("tihany_valid"):
                print(f"  → Tihany ✅")
            else:
                print(f"  → Tihany ❌ (Hoffman too weak)")
        else:
            print()


if __name__ == "__main__":
    main()
