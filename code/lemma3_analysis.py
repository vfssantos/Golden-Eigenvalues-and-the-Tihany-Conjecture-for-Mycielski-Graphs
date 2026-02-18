"""
Lemma 3 Analysis: φ-band cuts, Rayleigh quotients, and Kron reduction
Tests the three Golden Selection structures for Tihany:
  1. φ-ladder band cuts (Chapter 18)
  2. Kron reduction / Schur complement (Chapter 10) 
  3. Cheeger conductance (Chapter 12)
"""

import numpy as np
import networkx as nx
from itertools import combinations

PHI = (1 + np.sqrt(5)) / 2
GOLDEN_VALUES = [PHI, 1/PHI, PHI**2, 1/PHI**2, PHI**3, np.sqrt(5)]

def is_golden(val, tol=0.01):
    """Check if |val| is a golden value."""
    for g in GOLDEN_VALUES:
        if abs(abs(val) - g) < tol:
            return True
    return False

def hoffman_bound(A):
    """Hoffman chromatic bound: 1 - lmax/lmin."""
    eigs = np.linalg.eigvalsh(A)
    lmax, lmin = eigs[-1], eigs[0]
    if lmin >= -1e-10:
        return 1.0  # No negative eigenvalue
    return 1.0 - lmax / lmin

def build_mycielski_family():
    """Build Mycielski graphs M^k(C5) for k=0,1,2."""
    graphs = {}
    graphs['C5'] = (nx.cycle_graph(5), 3, 2)
    graphs['Grotzsch'] = (nx.mycielski_graph(4), 4, 2)
    graphs['Mycielski5'] = (nx.mycielski_graph(5), 5, 2)
    graphs['Icosahedron'] = (nx.icosahedral_graph(), 4, 3)
    return graphs

def golden_eigenvectors(A, tol=0.01):
    """Find all golden eigenvalues and their eigenvectors."""
    eigs, vecs = np.linalg.eigh(A)
    results = []
    for i, (e, v) in enumerate(zip(eigs, vecs.T)):
        if is_golden(e, tol):
            results.append((e, v))
    return results

def F_tau(A, x, tau, s, t):
    """Compute F(τ) = min{L_Hof(G[S]) - (s-1), L_Hof(G[T]) - (t-1)}."""
    n = len(x)
    S = [i for i in range(n) if x[i] >= tau]
    T = [i for i in range(n) if x[i] < tau]
    if len(S) < 2 or len(T) < 2:
        return -999, S, T
    
    A_S = A[np.ix_(S, S)]
    A_T = A[np.ix_(T, T)]
    
    hof_S = hoffman_bound(A_S)
    hof_T = hoffman_bound(A_T)
    
    return min(hof_S - (s-1), hof_T - (t-1)), S, T

def rayleigh_analysis(A, eigval, eigvec, S, T):
    """Analyze the Rayleigh quotient and correction term after cutting."""
    n = len(eigvec)
    x_S = eigvec[S]
    x_T = eigvec[T]
    
    A_SS = A[np.ix_(S, S)]
    A_ST = A[np.ix_(S, T)]
    A_TT = A[np.ix_(T, T)]
    
    # Rayleigh quotient on S
    if np.dot(x_S, x_S) > 1e-12:
        rq_S = np.dot(x_S, A_SS @ x_S) / np.dot(x_S, x_S)
        correction_S = np.dot(x_S, A_ST @ x_T) / np.dot(x_S, x_S)
    else:
        rq_S = 0
        correction_S = 0
    
    # Rayleigh quotient on T
    if np.dot(x_T, x_T) > 1e-12:
        rq_T = np.dot(x_T, A_TT @ x_T) / np.dot(x_T, x_T)
        correction_T = np.dot(x_T, A_ST.T @ x_S) / np.dot(x_T, x_T)
    else:
        rq_T = 0
        correction_T = 0
    
    return {
        'rayleigh_S': rq_S,
        'rayleigh_T': rq_T,
        'correction_S': correction_S,
        'correction_T': correction_T,
        'theoretical': eigval,  # Should be rq_S + correction_S
    }

def kron_reduction_spectrum(A, S, T):
    """Compute Kron reduction of Laplacian onto S, return spectrum."""
    n = A.shape[0]
    D = np.diag(A.sum(axis=1))
    L = D - A
    
    L_SS = L[np.ix_(S, S)]
    L_ST = L[np.ix_(S, T)]
    L_TT = L[np.ix_(T, T)]
    
    # Check invertibility of L_TT (add small regularization if needed)
    try:
        L_TT_inv = np.linalg.inv(L_TT)
        L_eff = L_SS - L_ST @ L_TT_inv @ L_ST.T
        kron_eigs = sorted(np.linalg.eigvalsh(L_eff))
    except np.linalg.LinAlgError:
        kron_eigs = None
    
    # Compare with induced subgraph Laplacian
    A_SS = A[np.ix_(S, S)]
    D_SS = np.diag(A_SS.sum(axis=1))
    L_induced = D_SS - A_SS
    induced_eigs = sorted(np.linalg.eigvalsh(L_induced))
    
    return kron_eigs, induced_eigs

def edge_cut_ratio(G, S, T):
    """Fraction of edges crossing the cut."""
    cut_edges = sum(1 for u, v in G.edges() if (u in set(S)) != (v in set(S)))
    return cut_edges / G.number_of_edges()

def conductance(G, S):
    """Conductance of the cut (S, V\S)."""
    S_set = set(S)
    cut = sum(1 for u, v in G.edges() if (u in S_set) != (v in S_set))
    vol_S = sum(dict(G.degree())[v] for v in S)
    vol_T = sum(dict(G.degree())[v] for v in G.nodes() if v not in S_set)
    return cut / min(vol_S, vol_T) if min(vol_S, vol_T) > 0 else float('inf')

# ═══════════════════════════════════════════════════════
# MAIN ANALYSIS
# ═══════════════════════════════════════════════════════

graphs = build_mycielski_family()

for name, (G, chi, omega) in graphs.items():
    n = G.number_of_nodes()
    A = nx.adjacency_matrix(G).toarray().astype(float)
    
    print(f"\n{'='*70}")
    print(f"  {name}: n={n}, χ={chi}, ω={omega}")
    print(f"{'='*70}")
    
    # All eigenvalues
    all_eigs = sorted(np.linalg.eigvalsh(A))
    golden_eigs = [(i, e) for i, e in enumerate(all_eigs) if is_golden(e)]
    print(f"Spectrum: {[f'{e:.3f}' for e in all_eigs]}")
    print(f"Golden eigenvalues: {[(f'{e:.4f}') for _, e in golden_eigs]}")
    
    # Target (s,t) pairs
    st_pairs = [(s, chi - s + 1) for s in range(2, chi) if chi - s + 1 >= 2]
    print(f"Tihany targets: {st_pairs}")
    
    # Get golden eigenvectors
    g_eigvecs = golden_eigenvectors(A)
    if not g_eigvecs:
        print("  No golden eigenvectors found.")
        continue
    
    for s, t in st_pairs:
        print(f"\n  --- (s,t) = ({s},{t}) ---")
        best_F = -999
        best_info = None
        
        for idx, (eigval, eigvec) in enumerate(g_eigvecs):
            # Sort eigenvector values to find all thresholds
            sorted_vals = sorted(set(np.round(eigvec, 10)))
            thresholds = []
            for i in range(len(sorted_vals) - 1):
                thresholds.append((sorted_vals[i] + sorted_vals[i+1]) / 2)
            
            for tau in thresholds:
                f_val, S, T = F_tau(A, eigvec, tau, s, t)
                if f_val > best_F:
                    best_F = f_val
                    best_info = {
                        'eigval': eigval,
                        'eigvec_idx': idx,
                        'tau': tau,
                        'S': S,
                        'T': T,
                        'f_val': f_val,
                        'S_size': len(S),
                        'T_size': len(T),
                    }
        
        if best_info is None:
            print(f"    No valid threshold found.")
            continue
            
        S = best_info['S']
        T = best_info['T']
        eigval = best_info['eigval']
        eigvec = g_eigvecs[best_info['eigvec_idx']][1]
        
        # Compute detailed analysis at optimal cut
        A_S = A[np.ix_(S, S)]
        A_T = A[np.ix_(T, T)]
        hof_S = hoffman_bound(A_S)
        hof_T = hoffman_bound(A_T)
        
        rayleigh = rayleigh_analysis(A, eigval, eigvec, S, T)
        kron_eigs, induced_eigs = kron_reduction_spectrum(A, S, T)
        phi_cond = conductance(G, S)
        cut_ratio = edge_cut_ratio(G, S, T)
        
        print(f"    Best F(τ) = {best_F:.4f} {'✅' if best_F >= 0 else '❌'}")
        print(f"    Eigenvector of λ = {eigval:.4f}, τ = {best_info['tau']:.4f}")
        print(f"    Split: |S|={len(S)}, |T|={len(T)}")
        print(f"    L_Hof(G[S]) = {hof_S:.4f} (need > {s-1})")
        print(f"    L_Hof(G[T]) = {hof_T:.4f} (need > {t-1})")
        print(f"    Rayleigh quotient on S: {rayleigh['rayleigh_S']:.4f} (eigenvalue: {eigval:.4f})")
        print(f"    Correction term S: {rayleigh['correction_S']:.4f}")
        print(f"    Rayleigh quotient on T: {rayleigh['rayleigh_T']:.4f}")
        print(f"    Correction term T: {rayleigh['correction_T']:.4f}")
        print(f"    Cut conductance: {phi_cond:.4f}")
        print(f"    Edge cut ratio: {cut_ratio:.4f}")
        
        # Golden eigenvalue preservation in induced subgraphs
        eigs_S = sorted(np.linalg.eigvalsh(A_S))
        eigs_T = sorted(np.linalg.eigvalsh(A_T))
        golden_S = [e for e in eigs_S if is_golden(e)]
        golden_T = [e for e in eigs_T if is_golden(e)]
        print(f"    Golden eigs in G[S]: {[f'{e:.4f}' for e in golden_S]}")
        print(f"    Golden eigs in G[T]: {[f'{e:.4f}' for e in golden_T]}")
        
        # Kron reduction comparison
        if kron_eigs is not None:
            golden_kron = [e for e in kron_eigs if is_golden(e, tol=0.05)]
            print(f"    Kron reduction golden eigs: {[f'{e:.4f}' for e in golden_kron]}")
            
            # Key comparison: induced vs Kron eigenvalues
            print(f"    Induced subgraph L eigs: {[f'{e:.3f}' for e in induced_eigs[:5]]}...")
            print(f"    Kron reduction L eigs:    {[f'{e:.3f}' for e in kron_eigs[:5]]}...")

