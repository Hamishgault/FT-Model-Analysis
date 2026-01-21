# Brubach et al. (2022) Paper Analysis and Implementation Comparison

## Paper Details
**Title:** Detailed Kinetic Modeling of CO2-Based Fischer–Tropsch Synthesis  
**Authors:** Lucas Brübach, Daniel Hodonj, Linus Biffar, Peter Pfeifer  
**Journal:** Catalysts 2022, 12(6), 630  
**DOI:** https://doi.org/10.3390/catal12060630  

---

## Paper's Model Architecture (Full 82-Equation System)

### Surface Species (10 types in model)
From Table 1 and mechanism description:
1. **Free sites (∗)**
2. **θ_H** (H* from dissociative H2 adsorption)
3. **θ_CO2** (CO2* from associative CO2 adsorption)
4. **θ_CO** (CO* from CO adsorption or CO2 dissociation)
5. **θ_OH** (OH* pseudo-component for oxygen-containing species)
6. **θ_CH2** (CH2* monomer for chain growth)
7. **θ_R** (R_n* - alkyl species, n=1 to n_max)
8. **θ_IR** (IR_n* - iso-alkyl species for branched chains)
9. **(Implicit: O* mentioned in two-step RWGS but merged with OH* via quasi-equilibrium)**

### Key Reactions from Paper (Table 1, 2, 3)

#### 1. H2 Adsorption (Equilibrium, Step 1)
```
H2 + 2* ⇌ 2H*
K1 = θ_H^2 / (θ_* ^2 * pH2)
```
- K1 (bar^-1) at T_ref = 300°C

#### 2. CO2 Adsorption (Equilibrium, Step 2)
```
CO2 + * ⇌ CO2*
K2 = θ_CO2 / (θ_* * pCO2)
```
- K2 (bar^-1) at T_ref = 300°C

#### 3. CO Adsorption/Desorption (Kinetic, Steps 3)
```
CO + * → CO*    (k3+)
CO* → CO + *    (k3-)
r_CO_ads_des = k3+ * pCO * θ_* - k3- * θ_CO
```
- k3+ (mol g^-1 h^-1 bar^-1), k3- (mol g^-1 h^-1)

#### 4. H2O Adsorption/Desorption (Kinetic, Step 4)
```
H2O + 2* ⇌ OH* + H*
r4 = k4+ * pH2O * θ_*^2 - k4- * θ_OH * θ_H
```
- k4+ (mol g^-1 h^-1 bar^-1), k4- (mol g^-1 h^-1)
- Temperature dependent: EA4+ = 151 kJ/mol, EA4- = 199 kJ/mol

#### 5. RWGS - Direct CO2 Dissociation (Step 5, Table 2)
Paper uses **direct CO2 dissociation mechanism** from DFT studies:
```
Step 5a: CO2* + * → CO* + O*      (quasi-equilibrium, K5a)
Step 5b: O* + H* → OH*            (quasi-equilibrium, K5b)
Overall: CO2* + H* → CO* + OH*
```
Expression (Table 1, step 5):
```
r5 = k5 * θ_CO2 * θ_H / θ_OH
```
- Paper states: "O* expressed in terms of θ_OH via step 5b quasi-equilibrium"
- k5 (mol g^-1 h^-1) at T_ref

**CRITICAL:** Paper uses **θ_H / θ_OH** term, NOT separate O* tracking!

#### 6. H-Assisted CO Dissociation (Step 6, Table 3)
```
Step 6a: CO* + H* ⇌ HCO*       (quasi-equilibrium, K6a)
Step 6b: HCO* + H* → CH2* + O* (RDS, k6b)
Overall: CO* + 2H* → CH2* + O*
```
Expression (Table 1, step 6):
```
r6 = k6 * K6a * θ_CO * θ_H^2 / θ_OH^phi
```
- k6 = k6b (rate-determining step constant, mol g^-1 h^-1)
- K6a (equilibrium constant for 6a)
- **phi exponent on θ_OH:** Paper uses θ_OH^phi damping factor

**CRITICAL:** Paper's final expression has **K6a * θ_H^2 / θ_OH^phi**, NOT separate HCO tracking!

#### 7. Chain Initiation (Step 7)
```
CH2* + H* → R1*
r7 = k7 * θ_CH2 * θ_H
```
- k7 (mol g^-1 h^-1)

#### 8. Chain Growth (Step 8)
```
R_n* + CH2* → R_(n+1)*
r8(n) = k8 * θ_R_n * θ_CH2
```
- k8 (mol g^-1 h^-1) assumed constant for all n

#### 9. Associative Desorption to n-alkane (Step 9)
```
R_n* + H* → P_n + *
r9(n) = k9 * exp(-Γ9*n/(RT)) * θ_R_n * θ_H
```
- k9 (mol g^-1 h^-1)
- Γ9 chain-length dependence parameter (kJ/mol)

#### 10. Dissociative Desorption to 1-alkene (Step 10, β-hydride elimination)
```
R_n* → Ol_n + H* + *
r10(n) = k10 * exp(-Γ10*n/(RT)) * θ_R_n / θ_OH^phi
```
- k10 (mol g^-1 h^-1)
- Γ10 chain-length dependence (kJ/mol)
- EA10 activation energy (kJ/mol) for temperature dependence
- **θ_OH^(-1) empirical damping factor** (see Section 3.1)

#### 11. Branching (Step 11)
```
R_n* → IR_n*
r11(n) = k11 * exp(-Γ11*n/(RT)) * θ_R_n * θ_CH2
```
- k11 (mol g^-1 h^-1)
- Γ11 chain-length dependence

#### 12. Iso-alkyl Growth (Step 12)
```
IR_n* + CH2* → IR_(n+1)*
r12(n) = k12 * θ_IR_n * θ_CH2
```

#### 13. Iso-termination to iso-alkene (Step 13)
```
IR_n* → I_n + H* + *
r13(n) = k13 * θ_IR_n / θ_OH^phi
```

### Temperature Dependencies (Equation 2)
```
k_j,n = k_j,Tref * exp(-Γ_j * n / (RT)) * exp(-EA_j/R * (1/T - 1/T_ref))
```
- T_ref = 300°C = 573 K
- Γ_j: chain-length dependence (for steps 9, 10, 11)
- EA_j: activation energy (only for steps 4+, 4-, 10)

### Site Balance (Equation 5)
```
1 = θ_* + θ_CO2 + θ_CO + θ_H + θ_OH + θ_CH2 + Σ_n θ_R_n + Σ_n θ_IR_n
```

### Steady-State Surface Coverage (Equation 4)
For each adsorbed species k:
```
dθ_k/dt = Σ_j ν_k,j * r_j = 0
```
Results in 82 algebraic equations (for all R_n and IR_n up to n_max).

### Parameter Values (Table 4)
Paper provides 23 adaptable parameters from regression:

| Parameter | Value | Units | Description |
|-----------|-------|-------|-------------|
| K1 | 5.82×10^-1 | bar^-1 | H2 adsorption equilibrium |
| K2 | 1.37 | bar^-1 | CO2 adsorption equilibrium |
| k3+ | 1.43×10^3 | mol g^-1 h^-1 bar^-1 | CO adsorption |
| k3- | 9.21×10^3 | mol g^-1 h^-1 | CO desorption |
| k4+ | 3.49×10^2 | mol g^-1 h^-1 bar^-1 | H2O adsorption |
| k4- | 1.45×10^1 | mol g^-1 h^-1 | H2O desorption |
| k5 | 2.21×10^1 | mol g^-1 h^-1 | RWGS effective |
| K6a | 4.10×10^1 | - | HCO equilibrium |
| k6 (k6b) | 1.64×10^-1 | mol g^-1 h^-1 | CO hydrogenation RDS |
| k7 | 9.57×10^2 | mol g^-1 h^-1 | Initiation |
| k8 | 1.79×10^3 | mol g^-1 h^-1 | Chain growth |
| k9 | 2.18×10^2 | mol g^-1 h^-1 | n-alkane termination |
| k10 | 2.49×10^3 | mol g^-1 h^-1 | 1-alkene termination |
| k11 | 7.67×10^2 | mol g^-1 h^-1 | Branching |
| k12 | (not listed) | mol g^-1 h^-1 | Iso-growth |
| k13 | (not listed) | mol g^-1 h^-1 | Iso-termination |
| Γ9 | (not listed) | kJ/mol | n-alkane chain-length dep |
| Γ10 | 4.52×10^2 | kJ/mol | 1-alkene chain-length dep |
| Γ11 | 2.89×10^3 | kJ/mol | Branching chain-length dep |
| phi | (not listed, likely ~1) | - | OH exponent damping |
| EA4+ | 1.51×10^2 | kJ/mol | H2O adsorption activation |
| EA4- | 1.99×10^2 | kJ/mol | H2O desorption activation |
| EA10 | 5.97 | kJ/mol | 1-alkene termination activation |

**Note:** Paper uses k9, k10, k11 without subscript 'a' or 'b' distinction for base rates.

---

## Current Implementation Analysis

### ✅ What's Correct

1. **Surface species structure (10 unknowns):** Current implementation tracks θ_*, θ_H, θ_CO2, θ_CO, θ_OH, θ_CH2, θ_R, θ_IR, θ_O, θ_HCO ✓

2. **H2 and CO2 adsorption equilibria (steps 1-2):** Implemented correctly with K1, K2 ✓

3. **CO adsorption/desorption kinetics (step 3):** k3p, k3m correctly applied ✓

4. **H2O adsorption/desorption (step 4):** k4p, k4m correctly applied ✓

5. **O* formation and hydrogenation:** Explicit tracking of O* formation and O* + H* → OH* ✓

6. **HCO formation/decomposition:** Explicit tracking with K6a quasi-equilibrium ✓

7. **Termination damping with Gamma10:** Applied exp(-Gamma10/(RT)) to termination rates ✓

8. **Rate unit conversion:** Converts mol g^-1 h^-1 to mol m^-3 s^-1 using cat_loading ✓

9. **Non-negative clamping:** Prevents NaN in fractional power operations ✓

10. **Robust solver strategy:** Bounded least_squares → root() → continuation → fallback ✓

### ❌ Major Discrepancies

#### 1. **RWGS Mechanism Expression (CRITICAL)**
**Paper (Table 1, step 5):**
```
r5 = k5 * θ_CO2 * θ_H / θ_OH
```
Uses θ_H / θ_OH term (O* implicitly via quasi-equilibrium with OH*)

**Current Implementation:**
```python
r_O_form = prm.k5 * thCO2 * ths          # CO2* + * → CO* + O*
r_O_to_OH = prm.k5m * thO * thH          # O* + H* → OH*
r_rwgs = r_O_to_OH
```
**Issue:** We're explicitly tracking O* and using k5m (not in paper). Paper's expression embeds O* balance via θ_H/θ_OH.

**Required Fix:** Eliminate explicit O* tracking; use paper's expression `r5 = k5 * θ_CO2 * θ_H / max(θ_OH, 1e-12)`.

---

#### 2. **CO Hydrogenation (RDS) Expression (CRITICAL)**
**Paper (Table 1, step 6):**
```
r6 = k6 * K6a * θ_CO * θ_H^2 / θ_OH^phi
```
Combines steps 6a (quasi-equilibrium) and 6b (RDS) into single expression using K6a.

**Current Implementation:**
```python
r_HCO_form = k6f * thCO * thH           # CO* + H* ⇌ HCO*
r_HCO_decomp = k6r * thHCO
r2 = prm.k6 * thHCO * thH               # HCO* + H* → CH2*
```
**Issue:** We're explicitly tracking HCO* and decomposition kinetics. Paper's expression uses **K6a * θ_H^2 / θ_OH^phi**, not separate HCO balance.

**Required Fix:** Eliminate explicit HCO* tracking; use paper's expression `r6 = k6 * K6a * θ_CO * θ_H^2 / max(θ_OH^phi, 1e-12)`.

---

#### 3. **Parameter Values Don't Match Table 4**
**Paper values vs Current defaults:**

| Parameter | Paper | Current | Match? |
|-----------|-------|---------|--------|
| K1 | 5.82×10^-1 | 5.82×10^-1 | ✓ |
| K2 | 1.37 | 1.37 | ✓ |
| k3p | 1.43×10^3 | 1.43×10^3 | ✓ |
| k3m | 9.21×10^3 | 9.21×10^3 | ✓ |
| k4p | 3.49×10^2 | 3.49×10^2 | ✓ |
| k4m | 1.45×10^1 | 1.45×10^1 | ✓ |
| k5 | 2.21×10^1 | 2.21×10^1 | ✓ |
| **K6a** | **4.10×10^1** | **1.0** | ❌ |
| **k6** | **1.64×10^-1** | **4.10×10^1** | ❌ |
| **k7** | **9.57×10^2** | **1.64×10^-1** | ❌ |
| **k8** | **1.79×10^3** | **1.0×10^4** | ❌ |
| **k9** | **2.18×10^2** | **k9a=1.0×10^3** | ❌ |
| **k10** | **2.49×10^3** | **k10a=9.57×10^2** | ❌ |
| k11 | 7.67×10^2 | k11a=2.18×10^2 | ❌ |
| Γ10 | 4.52×10^2 | 0.452×10^3 | ✓ |
| Γ11 | 2.89×10^3 | 2.89×10^3 | ✓ |

**Major issues:**
- K6a = 1.0 instead of 41.0
- k6 = 41.0 instead of 0.164 (values swapped!)
- k7 = 0.164 instead of 957
- k8 = 10,000 instead of 1790
- k9a = 1000 instead of 218 (CH4 termination)
- k10a = 957 instead of 2490 (alkene termination)

---

#### 4. **Chain-Length Tracking Missing**
**Paper:** Solves for individual θ_R_n and θ_IR_n for n=1 to n_max (82 equations).

**Current:** Lumps all R and IR species into single θ_R, θ_IR coverages.

**Impact:** Cannot reproduce chain-length-dependent termination rates (exp(-Γ_n/(RT))) or individual product distributions.

**Required Fix:** Expand algebraic system to track θ_R_1, θ_R_2, ..., θ_R_nmax separately.

---

#### 5. **Chain-Length-Dependent Rate Expressions Missing**
**Paper (steps 9, 10, 11):**
```
r9(n) = k9 * exp(-Γ9*n/(RT)) * θ_R_n * θ_H
r10(n) = k10 * exp(-Γ10*n/(RT)) * θ_R_n / θ_OH^phi
r11(n) = k11 * exp(-Γ11*n/(RT)) * θ_R_n * θ_CH2
```

**Current:** Uses constant k9a, k9b, k10a, k10b without chain-length exponential decay.

**Required Fix:** Implement chain-length-dependent rate expressions with Γ parameters.

---

#### 6. **Temperature Dependencies Incomplete**
**Paper (Equation 2):**
```
k_j,n = k_j,Tref * exp(-Γ_j*n/(RT)) * exp(-EA_j/R * (1/T - 1/T_ref))
```
Only steps 4+, 4-, 10 have activation energies.

**Current:** Only applies Gamma10 damping to terminations; no EA4+, EA4-, EA10 implementation.

**Required Fix:** Add Arrhenius temperature dependencies for k4+, k4-, k10.

---

#### 7. **Branching and Iso-chain Implementation**
**Paper:** Separate R_n → IR_n branching step (k11) with chain-length dependence.

**Current:** Uses single k11a without chain-length dependence; iso-growth (k12) and iso-termination (k13) exist but with incorrect parameter values.

**Required Fix:** Implement full chain-length-dependent branching and iso-chain tracking.

---

#### 8. **Product Lumping**
**Paper:** Resolves individual n-alkanes (P_n), 1-alkenes (Ol_n), iso-alkenes (I_n) up to C15+.

**Current:** Lumps products into CH4, C2_4, C5+ with fixed stoichiometry in nu matrix.

**Impact:** Cannot compare individual hydrocarbon selectivities or ASF distribution plots.

**Required Fix:** Export individual product rates for n=1 to n_max.

---

#### 9. **Site Balance**
**Paper (Equation 5):**
```
1 = θ_* + θ_CO2 + θ_CO + θ_H + θ_OH + θ_CH2 + Σ_n θ_R_n + Σ_n θ_IR_n
```

**Current:**
```python
bal_site = 1.0 - (ths + thH + thCO2 + thCO + thOH + thCH2 + thR + thIR + thO + thHCO)
```
**Issue:** Includes θ_O and θ_HCO, which should NOT be in site balance if using paper's quasi-equilibrium expressions.

**Required Fix:** Remove θ_O and θ_HCO from site balance; use paper's reduced set.

---

## Summary of Required Changes

### High Priority (Mechanism Correctness)

1. **Remove explicit O* and HCO* tracking from surface unknowns**
   - Reduce unknowns from 10 to 8: [θ_*, θ_H, θ_CO2, θ_CO, θ_OH, θ_CH2, Σθ_R_n, Σθ_IR_n]
   - Use paper's quasi-equilibrium expressions for RWGS and CO hydrogenation

2. **Fix RWGS expression to match paper:**
   ```python
   r5 = prm.k5 * thCO2 * thH / max(thOH, 1e-12)
   ```

3. **Fix CO hydrogenation (step 6) to match paper:**
   ```python
   r6 = prm.k6 * prm.K6a * thCO * thH**2 / max(thOH**prm.phi, 1e-12)
   ```

4. **Correct parameter values to match Table 4:**
   - K6a: 1.0 → 41.0
   - k6: 41.0 → 0.164
   - k7: 0.164 → 957
   - k8: 10,000 → 1,790
   - k9: 1,000 → 218 (base n-alkane termination)
   - k10: 957 → 2,490 (base 1-alkene termination)
   - k11: 218 → 767 (base branching)

5. **Expand to full chain-length tracking:**
   - Track θ_R_1, θ_R_2, ..., θ_R_nmax individually
   - Track θ_IR_1, θ_IR_2, ..., θ_IR_nmax individually
   - Results in ~82 algebraic equations for n_max=20

6. **Implement chain-length-dependent rate expressions:**
   ```python
   r9_n = prm.k9 * np.exp(-prm.Gamma9 * n / (R_gas * T)) * theta_R_n * thH
   r10_n = prm.k10 * np.exp(-prm.Gamma10 * n / (R_gas * T)) * theta_R_n / max(thOH**prm.phi, 1e-12)
   r11_n = prm.k11 * np.exp(-prm.Gamma11 * n / (R_gas * T)) * theta_R_n * thCH2
   ```

7. **Add temperature dependencies:**
   ```python
   k4p_T = prm.k4p * np.exp(-prm.EA4p / R_gas * (1/T - 1/573.15))
   k4m_T = prm.k4m * np.exp(-prm.EA4m / R_gas * (1/T - 1/573.15))
   k10_T = prm.k10 * np.exp(-prm.EA10 / R_gas * (1/T - 1/573.15))
   ```

### Medium Priority (Quantitative Accuracy)

8. **Correct site balance to paper's specification:**
   - Remove θ_O and θ_HCO from site balance

9. **Export individual product distributions:**
   - Return P_n (n-alkanes), Ol_n (1-alkenes), I_n (iso-alkenes) rates for each n

10. **Add missing Γ9 parameter** (chain-length dependence for n-alkane termination)

11. **Implement special handling for C3 branching** (paper mentions extra k11 value for C3)

### Low Priority (Model Completeness)

12. **Add secondary reactions** (optional, paper discarded due to computational cost):
    - Secondary hydrogenation of 1-alkenes
    - VLE for heavy products

13. **Implement 2-site model** (paper mentions as future improvement):
    - Separate RWGS and FTS sites to avoid empirical damping factor

---

## Recommendations

### Option A: Full 82-Equation Implementation (Most Accurate)
Pros:
- Matches paper exactly
- Can reproduce all figures and parity plots
- Individual hydrocarbon selectivities
- Proper chain-length effects

Cons:
- Large algebraic system (82 unknowns)
- Slower solver convergence
- More complex implementation

### Option B: Corrected Reduced Model (Balanced Approach)
Pros:
- Fixes mechanism expressions (RWGS, CO hydrogenation)
- Corrects parameter values
- Faster than full model
- Still captures key physics

Cons:
- Lumped chain tracking
- Cannot reproduce detailed ASF plots
- Missing chain-length dependencies

### Option C: Hybrid Approach (Recommended)
1. **Phase 1:** Fix expressions and parameters in current reduced model
   - Correct RWGS and CO hydrogenation expressions
   - Update parameters to Table 4 values
   - Remove O* and HCO* explicit tracking
   - Validate against paper's conversion/selectivity results

2. **Phase 2:** Expand to full chain-length tracking
   - Implement θ_R_n, θ_IR_n individual tracking
   - Add chain-length-dependent rates
   - Export individual product distributions
   - Validate against paper's ASF plots and hydrocarbon fractions

---

## Next Steps

1. **Implement Phase 1 fixes** (corrected reduced model)
2. **Validate against paper's experimental results:**
   - CO2 conversion (Figure 4a)
   - H2 conversion (Figure 4b)
   - CO selectivity (Figure 5a)
   - CH4 selectivity (Figure 5b)
   - Short-chain alkene selectivities (Figure 6)
3. **Assess need for Phase 2** based on user's modeling goals
4. **Add experimental data from paper's SI** for direct comparison

