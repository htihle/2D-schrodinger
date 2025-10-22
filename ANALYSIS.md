# Deep Analysis: Fundamental Issues with the 2D Schrödinger Solver

## TL;DR: The Core Problem

**The combination of implicit Crank-Nicolson + iterative solver + normalization is NOT preserving unitarity**, causing energy to vary wildly (even going negative!) instead of being conserved.

---

## Critical Findings

### 1. Energy Goes NEGATIVE (!!!)
- Initial energy: E ≈ 5.0
- During simulation: E drops to **-2.5** and oscillates wildly
- **This is physically impossible** for a harmonic oscillator
- Minimum possible energy is E_min = ℏω = 5.0 (ground state)

### 2. The Root Cause: Normalization + Implicit Solver Conflict

The Schrödinger equation is **unitary**, meaning:
- Probability ∫|ψ|²dx should be conserved AUTOMATICALLY
- Energy should be conserved AUTOMATICALLY
- NO manual normalization should be needed

**When we force normalization**, we're:
1. Changing the amplitude of ψ
2. This affects the energy (E ∝ |ψ|²)
3. Breaking the unitary evolution

### 3. Why the Iterative Solver Matters

We're using **Gauss-Seidel iteration** to solve the implicit equations:

```
ψ^(n+1) = solve_implicitly(ψ^n)
```

The Gauss-Seidel method:
- Converges in ~100-120 iterations
- But doesn't preserve the **symplectic/unitary structure** of the exact solution
- Each iteration slightly violates conservation laws
- Accumulates errors over many time steps

### 4. The Original Implementation Was Actually Better!

Your original code (before my "improvements"):
- Had the correct discretization
- Used Jacobi iteration (which is symmetric, unlike Gauss-Seidel)
- **Didn't normalize** (which was actually correct!)
- The only bugs were:
  1. Potential energy calculation (fixed)
  2. Unclear initial condition

---

## What SHOULD Work

### Option 1: Remove Normalization (Best Fix)
```cpp
// DELETE this line after each time step:
// normalize(psi, dx);  // ← Remove this!
```

The Schrödinger equation is norm-conserving. If probability drifts, it means the numerical scheme has errors, not that we need to renormalize.

### Option 2: Use a Symplectic Integrator

Instead of iterative Crank-Nicolson, use:
- **Split-operator method** (FFT-based)
- **Strang splitting**: exp(-iV·dt/2) · exp(-iT·dt) · exp(-iV·dt/2)
- This preserves unitarity exactly (up to machine precision)

### Option 3: Direct Solver (Not Iterative)

Solve the implicit system directly using a linear solver:
```cpp
// Form the matrix A·ψ^(n+1) = B·ψ^n
// Solve directly with Armadillo's solve()
psi_new = solve(A_matrix, B_matrix * psi_old);
```

This is exact (no iteration error) but requires matrix inversion.

---

## Why Your Original Code Made Sense

Looking back at your original implementation:

```cpp
// Original: No normalization
for(int l = 0;l<100;l++){  // Fixed iterations
    // Jacobi iteration (symmetric)
    Anew(i,j) = Aold(i,j) + α·(neighbors);
    Anew(i,j) /= (1 + (4+V)·α);
}
A = Anew;
// NO NORMALIZATION HERE  ← This was correct!
```

**This was actually the right approach!** The issues were:
1. ✅ Discretization: Correct
2. ✅ No normalization: Correct
3. ❌ Potential energy calc: Bug (we fixed this)
4. ⚠️ Fixed iterations: Could use convergence check
5. ⚠️ Initial condition: Unclear

---

## Recommended Fix

### Minimal Changes to Original Code:

1. **Keep the no-normalization approach** ✅
2. **Fix the potential energy bug** (done) ✅
3. **Use a proper initial state** (Gaussian)
4. **Add convergence check** (optional but good)
5. **Remove** all the normalization I added ❌

### Better Longterm: Split-Operator Method

```cpp
// Much simpler and more stable:
// 1. Apply potential: ψ *= exp(-i·V·dt/2)
// 2. Apply kinetic (via FFT): ψ = FFT^(-1)[exp(-i·k²·dt) · FFT[ψ]]
// 3. Apply potential: ψ *= exp(-i·V·dt/2)
```

This is:
- Exact up to machine precision
- Always norm-conserving
- Always energy-conserving (for time-independent H)
- Much faster (O(N log N) vs O(N² × iterations))

---

## What We Learned

1. **"Improving" working code can make it worse** if you don't understand why it worked
2. **Normalization is a red flag** - if you need it, something else is wrong
3. **Iterative solvers + unitarity = tricky** - they can break conservation laws
4. **Your original intuition was good** - the core algorithm was sound

---

## Conclusion

The fundamental issue is **NOT** with the physics or discretization, but with:
1. The iterative solver breaking unitarity
2. Manual normalization making it worse

**To fix:** Remove normalization and verify probability conserves naturally. If it doesn't, switch to a symplectic method.

---

## Files Status

Current implementation has these files:
- `main.cpp` - Implicit solver with normalization (UNSTABLE - energy goes negative)
- `visualize.py` - Visualization (works fine)
- `diagnostics.dat` - Shows the energy instability
- `wavefunction.dat` - The unstable evolution

**Next steps:** Revert to a simpler, more stable approach without normalization.
