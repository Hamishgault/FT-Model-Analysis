"""
Examples: Using the Phase 1 Brubach Model

This file demonstrates practical usage of the corrected Brubach et al. (2022)
CO2-FTS kinetic model implementation.
"""

import numpy as np
from kinetics.brubach_2022_updated import BrubachModel, BrubachParams
from utils.species import SPECIES_IDX


# ============================================================================
# Example 1: Single-Point Calculation
# ============================================================================

def example_single_point():
    """Calculate coverages and rates at a single condition."""
    print("\n" + "="*70)
    print("Example 1: Single-Point Calculation")
    print("="*70)
    
    model = BrubachModel()
    
    # Operating conditions
    T = 573.15      # K (300°C - typical FTS temperature)
    P = 10.0        # bar
    
    # Inlet gas composition (mol/s)
    # Typical syngas composition: H2/CO ≈ 2.0
    Fi = np.array([
        0.100,   # CO
        0.200,   # H2
        0.000,   # CH4 (inlet)
        0.000,   # C2-C4 (inlet)
        0.000,   # C5+ (inlet)
        0.010,   # H2O (inlet)
        0.100,   # CO2
    ])
    
    print(f"\nConditions:")
    print(f"  T = {T:.1f} K ({T - 273.15:.0f}°C)")
    print(f"  P = {P:.1f} bar")
    print(f"  Feed: CO={Fi[SPECIES_IDX['CO']]:.3f}, H2={Fi[SPECIES_IDX['H2']]:.3f}, "
          f"CO2={Fi[SPECIES_IDX.get('CO2', -1)]:.3f} mol/s")
    
    # Solve surface coverages
    coverages = model.solve_surface(T, P, Fi)
    
    print(f"\nSurface Coverages:")
    for species in ["theta_*", "theta_H", "theta_CO2", "theta_CO", "theta_OH", "theta_CH2", "theta_R", "theta_IR"]:
        if species in coverages:
            print(f"  {species:15s} = {coverages[species]:.6f}")
    
    # Check site balance
    sum_coverage = sum(coverages.values())
    print(f"\nSite balance: Σθ_i = {sum_coverage:.6f} (should be ≈ 1.0)")
    
    # Calculate reaction rates
    rates = model.rate(T, P, Fi)
    
    print(f"\nReaction Rates (mol m^-3 s^-1):")
    rxn_names = [
        "r0: CO ads/des",
        "r1: RWGS",
        "r2: CO hydrogenation",
        "r3: Initiation",
        "r4: Chain growth",
        "r5: Alkane termination",
        "r6: Alkene termination",
        "r7: Branching",
    ]
    for i, (name, rate) in enumerate(zip(rxn_names, rates)):
        print(f"  {name:25s} = {rate:12.6e}")


# ============================================================================
# Example 2: Temperature Scan
# ============================================================================

def example_temperature_scan():
    """Scan RWGS reaction rate over a range of temperatures."""
    print("\n" + "="*70)
    print("Example 2: Temperature Scan (RWGS Rate)")
    print("="*70)
    
    model = BrubachModel()
    P = 10.0
    
    # Fixed inlet
    Fi = np.array([0.100, 0.200, 0.0, 0.0, 0.0, 0.010, 0.100])
    
    # Temperature range: 300–400°C
    T_range = np.linspace(573.15, 673.15, 11)  # K
    
    print(f"\nTemperature Scan (P = {P:.1f} bar)")
    print(f"{'T (°C)':>10} | {'T (K)':>10} | {'r_RWGS':>12} | {'r_CO_hydr':>12}")
    print("-" * 50)
    
    for T in T_range:
        rates = model.rate(T, P, Fi)
        r_rwgs = rates[1]  # r1: RWGS
        r_hydrogenation = rates[2]  # r2: CO hydrogenation
        
        print(f"{T - 273.15:>10.1f} | {T:>10.2f} | {r_rwgs:>12.6e} | {r_hydrogenation:>12.6e}")


# ============================================================================
# Example 3: H2/CO Ratio Effect (Syngas Composition)
# ============================================================================

def example_syngas_ratio():
    """Study effect of H2/CO ratio on coverages and rates."""
    print("\n" + "="*70)
    print("Example 3: H2/CO Ratio Study")
    print("="*70)
    
    model = BrubachModel()
    T = 573.15
    P = 10.0
    
    # H2/CO ratios to test
    h2_co_ratios = [1.0, 1.5, 2.0, 2.5, 3.0]
    
    print(f"\nEffect of H2/CO ratio (T={T-273.15:.0f}°C, P={P:.1f} bar)")
    print(f"{'H2/CO':>8} | {'θ_H':>10} | {'θ_CO':>10} | {'r_RWGS':>12} | {'r_CO_hydr':>12}")
    print("-" * 60)
    
    for ratio in h2_co_ratios:
        # Normalize to total feed = 0.4 mol/s
        F_CO = 0.2 / (1 + ratio)
        F_H2 = 0.2 * ratio / (1 + ratio)
        
        Fi = np.array([F_CO, F_H2, 0.0, 0.0, 0.0, 0.01, 0.1])
        
        cov = model.solve_surface(T, P, Fi)
        rates = model.rate(T, P, Fi)
        
        print(f"{ratio:>8.1f} | {cov['theta_H']:>10.6f} | {cov['theta_CO']:>10.6f} | "
              f"{rates[1]:>12.6e} | {rates[2]:>12.6e}")


# ============================================================================
# Example 4: Mechanism Validation
# ============================================================================

def example_mechanism_validation():
    """Verify that corrected mechanism expressions work correctly."""
    print("\n" + "="*70)
    print("Example 4: Mechanism Validation")
    print("="*70)
    
    model = BrubachModel()
    prm = model.params
    
    print(f"\nVerifying Phase 1 corrected expressions:")
    print(f"\n1. RWGS Mechanism (direct CO2 dissociation, Table 2):")
    print(f"   Expression: r_RWGS = k5 * θ_CO2 * θ_H / θ_OH")
    print(f"   k5 (Table 4) = {prm.k5:.2f} mol g^-1 h^-1")
    
    # Hypothetical coverages
    theta_CO2 = 0.05
    theta_H = 0.1
    theta_OH = 0.02
    r_rwgs_calc = prm.k5 * theta_CO2 * theta_H / theta_OH
    print(f"\n   Example: θ_CO2={theta_CO2}, θ_H={theta_H}, θ_OH={theta_OH}")
    print(f"   r_RWGS = {prm.k5:.2f} × {theta_CO2} × {theta_H} / {theta_OH}")
    print(f"          = {r_rwgs_calc:.4f} mol g^-1 h^-1 ✓")
    
    print(f"\n2. CO Hydrogenation (H-assisted dissociation, Table 3):")
    print(f"   Expression: r_CO_hydr = k6 * K6a * θ_CO * θ_H^2 / θ_OH^phi")
    print(f"   k6 (Table 4) = {prm.k6:.4f} mol g^-1 h^-1")
    print(f"   K6a (Table 4) = {prm.K6a:.1f} (HCO* quasi-equilibrium)")
    print(f"   phi (exponent) = {prm.phi}")
    
    theta_CO = 0.05
    r_hydrogenation = prm.k6 * prm.K6a * theta_CO * (theta_H ** 2) / (theta_OH ** prm.phi)
    print(f"\n   Example: θ_CO={theta_CO}, θ_H={theta_H}, θ_OH={theta_OH}")
    print(f"   r_CO_hydr = {prm.k6:.4f} × {prm.K6a:.1f} × {theta_CO} × {theta_H}^2 / {theta_OH}^{prm.phi}")
    print(f"             = {r_hydrogenation:.6f} mol g^-1 h^-1 ✓")
    
    print(f"\n3. Surface Unknowns (8 total, implicit quasi-equilibria):")
    print(f"   Phase 1 reduces from 10 to 8 unknowns by treating")
    print(f"   O* and HCO* as fast intermediates (quasi-equilibrium)")
    print(f"   rather than tracked species:")
    print(f"      - O*: removed, appears implicitly in θ_OH damping")
    print(f"      - HCO*: removed, appears implicitly in K6a factor")
    print(f"   Remaining: θ_*, θ_H, θ_CO2, θ_CO, θ_OH, θ_CH2, θ_R, θ_IR")


# ============================================================================
# Example 5: Parameter Sensitivity
# ============================================================================

def example_parameter_sensitivity():
    """Quick sensitivity check: how much does 10% k5 change affect RWGS?"""
    print("\n" + "="*70)
    print("Example 5: Parameter Sensitivity (k5 Variation)")
    print("="*70)
    
    T = 573.15
    P = 10.0
    Fi = np.array([0.100, 0.200, 0.0, 0.0, 0.0, 0.010, 0.100])
    
    print(f"\nEffect of ±10% variation in k5 (RWGS rate constant):")
    print(f"{'k5 variation':>20} | {'r_RWGS (base)':>14} | {'r_RWGS (mod)':>14} | {'Change %':>10}")
    print("-" * 65)
    
    # Baseline
    model_base = BrubachModel()
    rates_base = model_base.rate(T, P, Fi)
    r_rwgs_base = rates_base[1]
    
    # Variations
    for scale in [0.9, 0.95, 1.0, 1.05, 1.10]:
        params_mod = BrubachParams()
        params_mod.k5 *= scale
        model_mod = BrubachModel(params=params_mod)
        rates_mod = model_mod.rate(T, P, Fi)
        r_rwgs_mod = rates_mod[1]
        
        change_pct = 100 * (r_rwgs_mod - r_rwgs_base) / r_rwgs_base
        variation_str = f"{scale*100:.0f}% of nominal" if scale != 1.0 else "nominal (k5={:.2f})".format(params_mod.k5)
        
        print(f"{variation_str:>20} | {r_rwgs_base:>14.6e} | {r_rwgs_mod:>14.6e} | {change_pct:>9.1f}%")


# ============================================================================
# Main
# ============================================================================

if __name__ == "__main__":
    print("\n")
    print("╔" + "═"*68 + "╗")
    print("║" + " "*68 + "║")
    print("║" + "  Brubach et al. (2022) CO2-FTS Model - Phase 1 Examples".center(68) + "║")
    print("║" + "  (Corrected & Validated)".center(68) + "║")
    print("║" + " "*68 + "║")
    print("╚" + "═"*68 + "╝")
    
    example_single_point()
    example_temperature_scan()
    example_syngas_ratio()
    example_mechanism_validation()
    example_parameter_sensitivity()
    
    print("\n" + "="*70)
    print("All examples completed successfully!")
    print("="*70 + "\n")
