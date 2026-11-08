# From a SMILES string to the Gibbs matrix

This is the exact math the code runs, stage by stage, from the characters you type to the Newton system that minimizes Gibbs free energy. Each stage names the source file that implements it.

```
"CCO"  ──tokenize──▶  G = (V, E, ord)  ──valence──▶  h ∈ ℕ^|V|  ──count──▶  a ∈ ℕ^E
                                                                              │
             μ°/RT ◀── NASA-7 or Joback(G) ◀──────────────────────────────────┤
                │                                                             ▼
                └──────────▶  min G(n) = Σ nᵢ μᵢ   s.t.  A n = b   ◀──  A = [a₁ … a_S]
```

---

## 1. String → molecular graph (`src/smiles.cpp`)

The parser is a single-pass state machine with state `(prev, pending_bond, branch_stack, ring_table[100])`. Each character either creates an atom vertex `v`, sets the pending bond, pushes or pops a branch, or opens or closes a ring slot. It produces:

- vertices `V`, with per-atom `(Z, q, isotope, chirality, bracket?, aromatic?, h_bracket)`
- edges `E ⊆ V×V`, with `ord(e) ∈ {1,2,3,4, ar}`, where `ar` is aromatic

An implicit bond (no symbol) between two aromatic atoms is `ar`; anywhere else it is single.

### 1.1 Implicit hydrogens

Let `σ(v) = Σ_{e∋v} w(e)`, where `w(e) = ord(e)` for ordinary bonds and `w(e) = 1` for `ar` bonds. Let `𝒱(Z)` be the OpenSMILES default valence set (C: {4}, N: {3,5}, S: {2,4,6}, …).

- **Bracket atom:** `h(v)` is whatever was written inside the brackets (0 if no H was given).
- **Aliphatic organic atom:** `h(v) = min{ν ∈ 𝒱(Z) : ν ≥ σ(v)} − σ(v)`, or 0 if no such ν exists.
- **Aromatic organic atom:** with `ν₀ = min 𝒱(Z)` and `f = ν₀ − σ(v)`: if `f ≥ 1`, the atom *needs a π bond* and `h(v) = f − 1`; otherwise `h(v) = 0`.

### 1.2 Kekulization as a perfect matching

Let `N = {v aromatic : v needs a π bond}`, and let `E_ar|_N` be the aromatic edges whose two endpoints are both in `N`. The kekulized structure is a **perfect matching** `M ⊆ E_ar|_N`:

```
ord(e) = 2  if e ∈ M,     ord(e) = 1  for e ∈ E_ar \ M
```

The matching is found by depth-first search with backtracking, always branching on the most-constrained vertex first. If no perfect matching exists, the parser emits `KekulizationFailed`. The classic case is `c1ccnc1`: five atoms need π bonds, and an odd count can never be perfectly matched.

### 1.3 Valence validation

After kekulization, let `τ(v) = Σ ord(e) + h(v)`. For a charged main-group atom, the bound comes from its isoelectronic neighbour `Z′ = Z − q` (so N⁺ is bounded like C, O⁻ like F, C⁻ like N). The rule is:

```
τ(v) > ν_max(Z′)   ⇒   ValenceExceeded(v, position)
```

Every violating atom gets its own entry in the error vector.

---

## 2. Graph → composition vector (`src/formula.cpp`)

For species `i`, let `𝓔` be the set of elements, ordered by atomic number. The composition is:

```
a_ji = |{v ∈ Vᵢ : Z(v) = j}|  +  [j = H] · Σ_{v ∈ Vᵢ} h(v)
```

Stacking the composition vectors as columns gives the **element–species matrix**:

```
A = [ a₁ a₂ … a_S ] ∈ ℤ^{E×S}
```

Example: water-gas shift. Species are `[C-]#[O+]`, `O`, `O=C=O`, `[H][H]`, and the element rows are C, H, O:

```
          CO  H2O  CO2  H2
     C  [  1   0    1    0 ]
A =  H  [  0   2    0    2 ]
     O  [  1   1    2    0 ]
```

Reaction stoichiometry is the null space of `A`. Every stoichiometric vector ν satisfies `A ν = 0`; here `null(A) = span{(−1, −1, 1, 1)ᵀ}`, which is CO + H₂O → CO₂ + H₂.

`check_balance` (`src/reaction.cpp`) evaluates `A ν` for a user-written reaction. Each nonzero component becomes its own `ElementImbalance{element j, delta (Aν)_j}`, and a nonzero `qᵀν` (where `q` is the net-charge vector) becomes a `ChargeImbalance`.

---

## 3. Composition → standard chemical potentials (`src/thermo.cpp`, `src/joback.cpp`)

Each species gets one `gᵢ(T) = μᵢ°(T)/RT`. Every species uses the same convention: `H°(298.15 K) = Δ_fH°` and `S°` is absolute.

**NASA-7** (database species, GRI-Mech 3.0 coefficients):

```
H°/RT = a₀ + a₁T/2 + a₂T²/3 + a₃T³/4 + a₄T⁴/5 + a₅/T
S°/R  = a₀ ln T + a₁T + a₂T²/2 + a₃T³/3 + a₄T⁴/4 + a₆
g     = H°/RT − S°/R
```

**Joback** (fallback, computed from the graph). Let `k` index groups and `c_k(G)` be the number of times group `k` occurs in the graph. Then:

```
Δ_fH° = 68.29 + Σ_k c_k ΔH_k          Δ_fG° = 53.88 + Σ_k c_k ΔG_k        [kJ/mol]
S°    = (Δ_fH° − Δ_fG°)/298.15 + Σ_j a_j S°_ref(j)
g(T)  ≈ Δ_fH°/RT − S°/R                    (ΔH, S frozen at 298 K; flagged ThermoApproximation)
```

`S°_ref(j)` is the per-atom entropy of element `j` in its reference state, for example S°(H₂)/2 for hydrogen.

This step is where the SMILES string reaches the Gibbs vector without a database lookup. Changing the SMILES changes `c_k(G)`, which changes `g`.

---

## 4. Gibbs minimisation (`src/gibbs.cpp`)

For an ideal-gas phase at pressure `P` with `N = Σ nᵢ`:

```
G(n)/RT = Σᵢ nᵢ ( cᵢ + ln(nᵢ/N) ),      cᵢ = gᵢ(T) + ln(P/P°)
minimise  G(n)   subject to   A n = b,   n ≥ 0,      where  b = A n⁰
```

### 4.1 Pre-processing
1. Drop species that contain an element with `b_j = 0`; they are fixed at `nᵢ = 0`.
2. Keep a linearly independent set of rows of `A`, using modified Gram–Schmidt with a 1e-10 relative cutoff. Example: for pure H₂O the H row is 2 × the O row, so one of the two is dropped.

### 4.2 Element-potential (dual) form
At the optimum, the Lagrange (KKT) conditions give `cᵢ + ln xᵢ = Σ_j a_ji π_j`, where `π_j = λ_j/RT` is the element potential and `xᵢ = nᵢ/N`. So every mole fraction is an explicit function of the duals:

```
xᵢ(π) = exp( aᵢᵀ π − cᵢ ),        nᵢ = N xᵢ
```

The unknowns are `u = (π₁…π_E, ν = ln N)`, so there are **E + 1** of them no matter how many species exist. The residual is:

```
F_j = ( N Σᵢ a_ji xᵢ − b_j ) / b_j        j = 1..E        (element balance, relative)
F_E =   Σᵢ xᵢ − 1                                         (mole fractions sum to 1)
```

### 4.3 The Newton matrix
With `R_jk = Σᵢ a_ji a_ki xᵢ` and `β_j = Σᵢ a_ji xᵢ`:

```
            ⎡ N R₁₁/b₁  …  N R₁E/b₁   N β₁/b₁ ⎤
            ⎢    ⋮            ⋮          ⋮    ⎥
J = ∂F/∂u = ⎢ N R_E1/b_E … N R_EE/b_E N β_E/b_E⎥ ,      J δ = −F
            ⎣   β₁     …    β_E          0    ⎦
```

`R = A diag(x) Aᵀ` is symmetric positive definite once dependent rows are removed. The bordered system is solved by Gaussian elimination with partial pivoting.

The step `δ` is clamped to at most 5 in log units. A backtracking line search then requires `‖F(u+λδ)‖² ≤ (1 − 10⁻⁴λ)‖F(u)‖²`. Convergence is declared when `‖F‖_∞ < 10⁻¹²`.

The dual form is used instead of a primal RAND step for one reason: in a primal method, a single trace species approaching zero pulls the fraction-to-boundary rule down for the entire step. In the dual form every `nᵢ = N·exp(·)` is positive by construction, so trace species such as CH₄ at 2500 K (8.7×10⁻¹⁵ mol) converge as quickly as major ones. That run converges in 10 iterations.

The initial guess is a least-squares fit, `(A Aᵀ) π⁰ = A (c + ln x⁰)` with `ν⁰ = ln Σ n⁰`.

### 4.4 Outputs
`nᵢ`, `xᵢ`, `π` (the equilibrium chemical potential of element `j` per RT), `G/RT = Σ nᵢ(cᵢ + ln xᵢ)`, and the achieved `max_j |A n − b|_j / b_j`.

---

## 5. How the matrix updates during a simulated reaction

What gets recomputed depends on what changed in `ChemicalSystemState`:

| Event | `A` | `g` | `b` | Newton system |
|---|---|---|---|---|
| New SMILES species added | append column `a_{S+1}` from §2 | append `g_{S+1}` from §3 | unchanged (n⁰ = 0) | same size E+1; the new species enters only through `xᵢ(π)` |
| Species introduces a new element | append row | — | append `b_j` | grows to E+2 |
| Temperature step `T → T′` | — | recompute every `gᵢ(T′)` | — | warm start from previous `π` |
| Pressure step | — | `cᵢ += ln(P′/P)` | — | warm start |
| Feed / extent `n⁰ → n⁰ + ν ξ` | — | — | `b` unchanged if `Aν = 0`; otherwise `ElementImbalance` | — |

Because the solver's dimension is `E + 1`, not `S`, adding species to a reacting system is cheap. A mechanism with hundreds of C/H/O/N species is still a 5×5 solve per Newton step, plus O(E²·S) work to assemble `R`.

---

## 6. Worked example (reproducible)

```bash
./build/omnichem equilibrate -T 1000 -P 1 '[C-]#[O+]' 1 O 1 O=C=O 0 '[H][H]' 0
```

Because the shift reaction conserves moles, `x_CO₂ x_H₂ / (x_CO x_H₂O) = K = exp(−Σ νᵢ gᵢ)`. The test `water_gas_shift_matches_equilibrium_constant` asserts that the solver's composition reproduces `K` to 10⁻⁸, which verifies §4 end to end. At 1000 K the result is n_CO₂ = n_H₂ = 0.5451 mol, so K ≈ 1.44, matching the literature value.

---

## 7. Kinetics: from a mechanism to the stiff ODE (`src/kinetics.cpp`)

Each species in a `.ocm` mechanism is bound to a SMILES string. Element balance is therefore checked against the same composition matrix `A` as §2, and thermochemistry comes from §3.

**Rate of progress** of reaction `r`, with concentrations `C` in mol/m³:

```
q_r = [M]_r · ( k_f,r Π_j C_j^ν'_jr  −  k_r,r Π_j C_j^ν''_jr )
k_f = A T^b exp(−E_a/RT)                         [M]_r = Σ_j ε_jr C_j   (third-body reactions, else 1)
k_r = k_f / K_c,     ln K_c = −Σ_j (ν''_jr − ν'_jr) g_j(T)  +  Δν_r ln(P°/RT)
ω_i = Σ_r (ν''_ir − ν'_ir) q_r
```

Because every `k_r` comes from detailed balance with the same `g_j(T)` used in §4, every reversible mechanism's long-time limit is exactly the Gibbs-minimiser equilibrium. `kinetics_relaxes_to_gibbs_equilibrium` and the ignition test check this to 10⁻⁶ and 10⁻⁴ respectively.

**Pressure-dependent rate constants** (`k_f` above):

```
Fall-off:              k = k_∞ · Pr/(1+Pr) · F,     Pr = k_0 [M] / k_∞
Chemically activated:  k = k_0 · 1/(1+Pr) · F
Lindemann F = 1.  Troe:  log₁₀F = log₁₀F_cent / (1 + f²),   f = (log₁₀Pr + c)/(n − 0.14(log₁₀Pr + c))
  F_cent = (1−a) e^(−T/T3) + a e^(−T/T1) + e^(−T2/T),   c = −0.4 − 0.67 log₁₀F_cent,   n = 0.75 − 1.27 log₁₀F_cent
PLOG:  ln k = ln k_i + (ln k_{i+1} − ln k_i) · (ln P − ln P_i)/(ln P_{i+1} − ln P_i), clamped outside [P_1, P_n]
```

For fall-off and PLOG reactions, `k_r = k_f/K_c` scales with the same factor, so `q = k_f(fwd − rev/K_c)`. Their Jacobians gain `ε_j · (dk_f/d[M]) · (fwd − rev/K_c)` and `RT · (dk_f/dP) · (fwd − rev/K_c)` respectively, where `P = Σ C_i RT`. Both derivatives are analytic, including `dF/dPr` for Troe.

**Jacobian** (analytic, exact):

```
∂q_r/∂C_j = [M] ( k_f ∂Π'/∂C_j − k_r ∂Π''/∂C_j ) + ε_jr ( k_f Π' − k_r Π'' )
∂ω_i/∂C_j = Σ_r (ν''_ir − ν'_ir) ∂q_r/∂C_j
```

**Reactor models.** The state is `y = (n_1..n_S [, T])` in moles. For constant volume, `V = N₀RT₀/P₀`; for constant pressure, `V = NRT/P₀`. Then:

```
dn_i/dt = V ω_i(n/V)
adiabatic, constant V:  Σ n_i c_v,i · dT/dt = −Σ u_i dn_i/dt      (u = h − RT, c_v = c_p − R)
adiabatic, constant P:  Σ n_i c_p,i · dT/dt = −Σ h_i dn_i/dt
```

The reactor Jacobian combines three parts:
- the species block `∂ω/∂C`, which is exact at constant volume because the V factors cancel;
- one forward difference in T, giving the T column (k(T), K_c(T), V(T));
- the T row, computed analytically:

```
∂Ṫ/∂n_j = −(Σ_i e_i J_ij)/D + Ṫ_num · c_j / D²,    D = Σ n_i c_i
```

## 8. Variable-order BDF (`src/ode.cpp`)

Let the corrector polynomial `p` of order `k` interpolate `y_{n+1}, y_n, …, y_{n+1−k}` at their actual, unequal times. BDF requires `p'(t_{n+1}) = f(t_{n+1}, y_{n+1})`. With Lagrange-basis derivatives `c_j = ℓ_j'(t_{n+1})`, this is:

```
G(y) = c_0 y + Σ_{j≥1} c_j y_{n+1−j} − f(t_{n+1}, y) = 0
c_0 = Σ_{i≥1} 1/(t_0 − t_i),     c_j = Π_{i≠0,j}(t_0 − t_i) / Π_{i≠j}(t_j − t_i)
```

It is solved by modified Newton with `(c_0 I − J)` factored once per attempt. The predictor extrapolates through `k+1` past points. The local error is estimated by Milne's device: with `|C_BDF| = 1/(k+1)` and `|C_pred| = 1`, `LTE ≈ ‖y_corr − y_pred‖_w / (k+2)`, in the weighted RMS norm `w_i = rtol·|y_i| + atol_i`.

After each accepted step, the solver estimates the error the neighbouring orders `k±1` would have had from the same data. It switches order when another order allows a step at least 1.3× (lower) or 1.4× (higher) larger. Dense output reuses the corrector polynomial, so requested output times never shorten steps.

Benchmarks in the test suite:

| Problem | Result |
|---|---|
| Prothero–Robinson, λ = −10⁶ | max error < 10⁻⁶ in < 3000 steps (an explicit method needs ~10⁷) |
| Robertson to t = 4·10⁵ (stiffness ratio ~10¹¹) | matches Hairer–Wanner reference values at t = 40 to 10⁻⁶; < 2000 steps |
| Adiabatic H₂/O₂/Ar ignition, 9 species | about 600 steps, no Newton failures; U and elements conserved to 10⁻⁶ and 10⁻⁹ |
