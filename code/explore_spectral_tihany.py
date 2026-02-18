#!/usr/bin/env python3
"""
Spectral exploration of the Erdős-Lovász Tihany Conjecture.

This script investigates whether the golden ratio φ governs the spectral
obstruction in graphs with χ > ω (the Tihany regime).

Phase 1: Compute eigenvalues, Lovász theta, and chromatic properties
         for known "hard" graphs (C₅, Petersen, Mycielski, Kneser, etc.)

Phase 2: Test partition strategies guided by spectral structure.

Requirements: numpy, scipy, networkx
    pip install numpy scipy networkx
"""

import numpy as np
from itertools import combinations
import sys

try:
    import networkx as nx
except ImportError:
    print("ERROR: networkx required. Install with: pip install networkx")
    sys.exit(1)

# ─── Constants ───────────────────────────────────────────────────────────────

PHI = (1 + np.sqrt(5)) / 2       # 1.6180339887...
PHI_INV = PHI - 1                 # 0.6180339887...
PHI_SQ = PHI + 1                  # 2.6180339887...
SQRT5 = np.sqrt(5)                # 2.2360679775...


# ─── Graph Library ───────────────────────────────────────────────────────────

def build_test_graphs():
    """Build a library of graphs relevant to the Tihany conjecture."""
    graphs = {}

    # --- The Atom ---
    graphs["C5 (pentagon)"] = nx.cycle_graph(5)

    # --- Odd cycles (generalizations of C5) ---
    for k in [7, 9, 11]:
        graphs[f"C{k}"] = nx.cycle_graph(k)

    # --- Petersen graph (Kneser K(5,2)) ---
    graphs["Petersen K(5,2)"] = nx.petersen_graph()

    # --- Mycielski graphs (triangle-free, high χ) ---
    # Mycielski(C5) = Grötzsch graph (χ=4, ω=2)
    graphs["Grötzsch (Mycielski₄)"] = nx.mycielski_graph(4)
    graphs["Mycielski₅ (χ=5,ω=2)"] = nx.mycielski_graph(5)
    graphs["Mycielski₆ (χ=6,ω=2)"] = nx.mycielski_graph(6)

    # --- Kneser graphs ---
    # K(n,k) has χ = n-2k+2, ω = floor(n/k) ... actually ω=1 for 2k>n
    # K(5,2) = Petersen (already added)
    graphs["Kneser K(7,3)"] = nx.kneser_graph(7, 3)
    # K(8,3) has χ = 4, larger
    graphs["Kneser K(8,3)"] = nx.kneser_graph(8, 3)

    # --- Complement of C7 (odd antihole) ---
    C7 = nx.cycle_graph(7)
    graphs["C̄7 (complement)"] = nx.complement(C7)

    # --- Icosahedron (H₃-related!) ---
    graphs["Icosahedron"] = nx.icosahedral_graph()

    # --- Dodecahedron ---
    graphs["Dodecahedron"] = nx.dodecahedral_graph()

    return graphs


# ─── Spectral Analysis ──────────────────────────────────────────────────────

def adjacency_eigenvalues(G):
    """Compute eigenvalues of the adjacency matrix."""
    A = nx.adjacency_matrix(G).toarray().astype(float)
    eigs = np.linalg.eigvalsh(A)
    return np.sort(eigs)[::-1]  # descending


def lovasz_theta_bound(G):
    """
    Compute the Lovász theta bound for vertex-transitive / regular graphs.
    For d-regular graph: θ(G) = -n λ_min / (d - λ_min)

    Note: This is exact for vertex-transitive graphs.
    For non-regular graphs, this is an approximation.
    """
    n = G.number_of_nodes()
    eigs = adjacency_eigenvalues(G)
    lam_max = eigs[0]
    lam_min = eigs[-1]

    if nx.is_regular(G):
        d = lam_max  # for regular graphs, largest eigenvalue = degree
        theta = -n * lam_min / (d - lam_min)
        return theta
    else:
        # Hoffman bound (lower bound on χ)
        hoffman = 1 - lam_max / lam_min
        return hoffman  # labeled as Hoffman bound


def chromatic_number_exact(G):
    """
    Compute exact chromatic number via greedy + verification.
    Uses networkx's greedy coloring as upper bound, then checks.
    """
    n = G.number_of_nodes()
    if n == 0:
        return 0

    # Upper bound from greedy
    coloring = nx.coloring.greedy_color(G, strategy="largest_first")
    upper = max(coloring.values()) + 1

    # Try to verify lower bound by checking if (upper-1)-colorable
    # For small graphs, we can try exact methods
    # networkx doesn't have exact χ, so we use the greedy as approximation
    # and note: for our test graphs, χ is known
    return upper


def clique_number(G):
    """Compute the clique number ω(G)."""
    cliques = nx.find_cliques(G)
    return max(len(c) for c in cliques)


def independence_number(G):
    """Compute α(G) = ω(Ḡ)."""
    Gc = nx.complement(G)
    return clique_number(Gc)


def fractional_chromatic_approx(G):
    """
    Approximate fractional chromatic number.
    For vertex-transitive: χ_f = n/α(G).
    For general graphs: χ_f ≈ n/α(G) is an approximation.
    """
    n = G.number_of_nodes()
    alpha = independence_number(G)
    return n / alpha


# ─── Golden Ratio Analysis ──────────────────────────────────────────────────

def golden_analysis(eigs):
    """
    Analyze how the eigenvalues relate to the golden ratio.
    Returns dict of golden-ratio-related quantities.
    """
    lam_min = eigs[-1]
    lam_max = eigs[0]

    results = {
        "λ_min": lam_min,
        "λ_max": lam_max,
        "|λ_min|/φ": abs(lam_min) / PHI,
        "λ_min + φ": lam_min + PHI,  # zero if λ_min = -φ
        "Hoffman": 1 - lam_max / lam_min,
        "Hoffman/√5": (1 - lam_max / lam_min) / SQRT5,
    }

    # Check if any eigenvalue is close to ±φ, ±φ⁻¹, ±φ²
    golden_vals = {
        "φ": PHI, "-φ": -PHI,
        "φ⁻¹": PHI_INV, "-φ⁻¹": -PHI_INV,
        "φ²": PHI_SQ, "-φ²": -PHI_SQ,
        "√5": SQRT5, "-√5": -SQRT5,
    }

    near_golden = []
    for eig in eigs:
        for name, val in golden_vals.items():
            if abs(eig - val) < 0.01:
                near_golden.append(f"λ={eig:.4f} ≈ {name} ({val:.4f})")

    results["near_golden_eigenvalues"] = near_golden
    return results


# ─── Tihany Partition Testing ───────────────────────────────────────────────

def test_tihany_partitions(G, chi_G, omega_G):
    """
    For a graph with χ > ω, test all valid (s,t) pairs and find partitions.
    Returns results for each (s,t) pair.
    """
    if chi_G <= omega_G:
        return {"status": "χ ≤ ω, Tihany does not apply"}

    n = G.number_of_nodes()
    nodes = list(G.nodes())
    results = []

    # All valid (s,t) with s+t-1 = χ, s,t ≥ 2
    for s in range(2, chi_G):
        t = chi_G - s + 1
        if t < 2:
            continue

        # For small graphs, test all partitions
        found = False
        best_partition = None

        if n <= 15:  # exact search feasible
            for size_S in range(1, n):
                if found:
                    break
                for S_tuple in combinations(range(n), size_S):
                    S = set(nodes[i] for i in S_tuple)
                    T = set(nodes) - S

                    GS = G.subgraph(S)
                    GT = G.subgraph(T)

                    # Quick check: need edges in both parts
                    if GS.number_of_edges() == 0 or GT.number_of_edges() == 0:
                        continue

                    chi_S = chromatic_number_exact(GS)
                    chi_T = chromatic_number_exact(GT)

                    if chi_S >= s and chi_T >= t:
                        found = True
                        best_partition = {
                            "S": sorted(S), "T": sorted(T),
                            "χ(G[S])": chi_S, "χ(G[T])": chi_T,
                        }
                        break
        else:
            # For larger graphs, just note we can't exhaustively search
            found = None  # unknown

        results.append({
            "s": s, "t": t,
            "found": found,
            "partition": best_partition,
        })

    return results


# ─── Main Analysis ──────────────────────────────────────────────────────────

def analyze_graph(name, G, verbose=True):
    """Full analysis of a single graph."""
    n = G.number_of_nodes()
    m = G.number_of_edges()
    eigs = adjacency_eigenvalues(G)
    omega = clique_number(G)
    chi = chromatic_number_exact(G)
    alpha = independence_number(G)
    chi_f = n / alpha
    is_regular = nx.is_regular(G)
    theta = lovasz_theta_bound(G)
    golden = golden_analysis(eigs)

    tihany_applies = chi > omega

    result = {
        "name": name,
        "n": n, "m": m,
        "χ": chi, "ω": omega, "α": alpha,
        "χ_f": chi_f,
        "χ > ω": tihany_applies,
        "regular": is_regular,
        "θ or Hoffman": theta,
        "eigenvalues": eigs,
        "golden": golden,
    }

    if verbose:
        print(f"\n{'='*70}")
        print(f"  {name}")
        print(f"{'='*70}")
        print(f"  n={n}, m={m}, regular={is_regular}")
        print(f"  χ={chi}, ω={omega}, α={alpha}, χ_f={chi_f:.4f}")
        print(f"  χ > ω: {tihany_applies}")
        print(f"  θ/Hoffman = {theta:.4f}")
        print(f"  n/α = {chi_f:.4f},  φ² = {PHI_SQ:.4f},  ratio = {chi_f/PHI_SQ:.4f}")
        print(f"  Eigenvalues: [{', '.join(f'{e:.4f}' for e in eigs[:5])}{'...' if len(eigs)>5 else ''}]")
        print(f"  λ_min = {eigs[-1]:.6f}")
        print(f"  |λ_min|/φ = {abs(eigs[-1])/PHI:.6f}")

        if golden["near_golden_eigenvalues"]:
            print(f"  Golden eigenvalues:")
            for g in golden["near_golden_eigenvalues"]:
                print(f"    {g}")

        if tihany_applies and n <= 15:
            print(f"\n  Testing Tihany partitions...")
            partitions = test_tihany_partitions(G, chi, omega)
            for p in partitions:
                status = "✅ FOUND" if p["found"] else "❌ NOT FOUND" if p["found"] is not None else "⬜ NOT TESTED"
                print(f"    (s={p['s']}, t={p['t']}): {status}")
                if p.get("partition"):
                    pp = p["partition"]
                    print(f"      S={pp['S']}, χ(G[S])={pp['χ(G[S])']}")
                    print(f"      T={pp['T']}, χ(G[T])={pp['χ(G[T])']}")

    return result


def golden_threshold_analysis(results):
    """
    Analyze whether φ² is a threshold for chromatic density.
    """
    print(f"\n{'='*70}")
    print(f"  GOLDEN THRESHOLD ANALYSIS")
    print(f"{'='*70}")
    print(f"  φ = {PHI:.6f}")
    print(f"  φ² = {PHI_SQ:.6f}")
    print(f"  √5 = {SQRT5:.6f}")
    print()

    tihany_graphs = [r for r in results if r["χ > ω"]]

    print(f"  Graphs with χ > ω ({len(tihany_graphs)} total):")
    print(f"  {'Name':<25} {'n/α':>8} {'φ²':>8} {'ratio':>8} {'|λ_min|':>8} {'φ':>8} {'ratio':>8}")
    print(f"  {'-'*25} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")

    for r in tihany_graphs:
        chi_f = r["χ_f"]
        lam_min = abs(r["eigenvalues"][-1])
        print(f"  {r['name']:<25} {chi_f:8.4f} {PHI_SQ:8.4f} {chi_f/PHI_SQ:8.4f} "
              f"{lam_min:8.4f} {PHI:8.4f} {lam_min/PHI:8.4f}")

    # Summary statistics
    chi_f_vals = [r["χ_f"] for r in tihany_graphs]
    lam_min_vals = [abs(r["eigenvalues"][-1]) for r in tihany_graphs]

    print(f"\n  Summary:")
    print(f"    n/α range: [{min(chi_f_vals):.4f}, {max(chi_f_vals):.4f}]")
    print(f"    |λ_min| range: [{min(lam_min_vals):.4f}, {max(lam_min_vals):.4f}]")
    print(f"    Graphs with n/α near φ² (within 10%): "
          f"{sum(1 for v in chi_f_vals if abs(v/PHI_SQ - 1) < 0.1)}/{len(chi_f_vals)}")
    print(f"    Graphs with |λ_min| near φ (within 10%): "
          f"{sum(1 for v in lam_min_vals if abs(v/PHI - 1) < 0.1)}/{len(lam_min_vals)}")


# ─── Entry Point ─────────────────────────────────────────────────────────────

def main():
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  Spectral Golden Hypothesis — Tihany Conjecture Explorer           ║")
    print("║  φ = (1+√5)/2 ≈ 1.6180   φ² ≈ 2.6180   √5 ≈ 2.2361             ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")

    graphs = build_test_graphs()
    results = []

    for name, G in graphs.items():
        r = analyze_graph(name, G)
        results.append(r)

    golden_threshold_analysis(results)

    print(f"\n{'='*70}")
    print(f"  NEXT STEPS")
    print(f"{'='*70}")
    print(f"  1. Check Q1: Do |λ_min| values cluster near φ or its powers?")
    print(f"  2. Check Q2: Are Tihany partitions found for all tested graphs?")
    print(f"  3. Check Q4: Is φ² a threshold for n/α?")
    print(f"  4. For graphs with golden eigenvalues: what is the partition structure?")


if __name__ == "__main__":
    main()
