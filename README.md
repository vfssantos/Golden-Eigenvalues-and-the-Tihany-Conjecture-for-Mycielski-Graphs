# A Golden-Ratio Spectral Approach to the Erdős-Lovász Tihany Conjecture

## Status: 🟡 Active Research

**Conjecture (Erdős-Lovász, 1968)**: For every graph $G$ and integers $s, t \geq 2$ with $\chi(G) = s + t - 1 > \omega(G)$, there exists a vertex partition $(S, T)$ with $\chi(G[S]) \geq s$ and $\chi(G[T]) \geq t$.

**Our approach**: The golden ratio $\phi = (1+\sqrt{5})/2$ governs the spectral structure of graphs where $\chi > \omega$. We develop a spectral partition method where eigenvectors of golden eigenvalues provide valid Tihany partitions, certified by the Hoffman bound.

---

## Contents

| File | Description |
|---|---|
| `main.md` | Full working paper: context, lemmas, evidence, strategy |
| `code/explore_spectral_tihany.py` | Phase 1: spectral survey of test graphs |
| `code/test_spectral_cuts.py` | Phase 2: eigenvector-threshold partition testing |
| `code/requirements.txt` | Python dependencies |

## Quick Results

- **Lemma 1** (Golden Propagation): Mycielski functor propagates $\mathbb{Q}(\sqrt{5})$ eigenvalues. Empirically verified; algebraic proof in progress.
- **Lemma 2** (Cut-to-Chromatic Transfer): Hoffman/Lovász-theta bounds certify Tihany partitions. **Proven.**
- **Lemma 3** (Golden Band Cut): Golden eigenvector thresholds provide valid cuts. Computationally verified for C₅, Petersen, Mycielski₅, icosahedron, dodecahedron. Proof strategy via φ-ladder bands + Kron reduction.

## Key Empirical Finding

The $\lambda_{\min}$-eigenvector does **not** give valid Tihany cuts. But eigenvectors of **golden eigenvalues** ($\phi, \phi^{-1}, \sqrt{5}$) do — consistently across all tested graphs.

## Origin

This work grew from the [Golden Selection](../The_Golden_Selection/v2/) theory, which derives Standard Model physics from a quasicrystalline framework where $\phi$ governs spectral structure. The mathematical tools developed there (Kron reduction, φ-ladder bands, Cheeger certification) suggested the spectral approach to Tihany.
