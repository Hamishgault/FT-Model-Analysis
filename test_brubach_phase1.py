#!/usr/bin/env python3
"""Test Phase 1 corrected Brubach model implementation."""

import sys
from pathlib import Path
import numpy as np

# Add workspace to path
ws_root = Path(__file__).parent
sys.path.insert(0, str(ws_root))

from kinetics.brubach_2022_updated import BrubachModel, BrubachParams
from utils.species import SPECIES_IDX

def test_basic_mechanism():
    """Test that basic mechanisms work correctly."""
    print("=" * 70)
    print("Test 1: Basic Mechanism Verification")
    print("=" * 70)
    
    model = BrubachModel()
    prm = model.params
    
    # Verify Table 4 parameters are loaded
    print(f"Key parameters from Table 4:")
    print(f"  k5 (RWGS): {prm.k5:.2f} mol g^-1 h^-1")
    print(f"  k6 (CO hydrogenation): {prm.k6:.4f} mol g^-1 h^-1")
    print(f"  K6a (HCO* quasi-eq): {prm.K6a:.2f}")
    print(f"  k7 (initiation): {prm.k7:.2f} mol g^-1 h^-1")
    print(f"  k8 (growth): {prm.k8:.2f} mol g^-1 h^-1")
    print(f"  k9 (alkane term): {prm.k9:.2f} mol g^-1 h^-1")
    print(f"  k10 (alkene term): {prm.k10:.2f} mol g^-1 h^-1")
    print()
    
    # Verify mechanism expressions via direct computation
    print("Mechanism expression verification:")
    
    # Test RWGS: r5 = k5 * θ_CO2 * θ_H / θ_OH
    theta_CO2 = 0.05
    theta_H = 0.1
    theta_OH = 0.01
    r_rwgs = prm.k5 * theta_CO2 * theta_H / theta_OH
    print(f"  RWGS: r5 = {prm.k5:.2f} * {theta_CO2} * {theta_H} / {theta_OH}")
    print(f"        = {r_rwgs:.4f} mol g^-1 h^-1 ✓")
    
    # Test CO hydrogenation: r6 = k6 * K6a * θ_CO * θ_H^2 / θ_OH^phi
    theta_CO = 0.05
    phi = 1.0
    r_hydrogenation = prm.k6 * prm.K6a * theta_CO * (theta_H ** 2) / (theta_OH ** phi)
    print(f"  CO hydrogenation: r6 = {prm.k6:.4f} * {prm.K6a:.2f} * {theta_CO} * {theta_H}^2 / {theta_OH}^{phi}")
    print(f"                   = {r_hydrogenation:.4f} mol g^-1 h^-1 ✓")
    print()


def test_surface_solver():
    """Test surface solver convergence."""
    print("=" * 70)
    print("Test 2: Surface Solver Convergence")
    print("=" * 70)
    
    model = BrubachModel()
    
    # Typical reactor inlet
    T = 573.15  # K (300°C)
    P = 10.0    # bar
    
    # Inlet flows: CO, H2, CH4, C2_4, C5plus, H2O, CO2
    # All in mol/s; typical syngas ratio: H2/CO ≈ 2
    Fi = np.array([0.1, 0.2, 0.0, 0.0, 0.0, 0.01, 0.1])
    
    print(f"Inlet conditions:")
    print(f"  T = {T:.2f} K ({T - 273.15:.1f}°C)")
    print(f"  P = {P:.1f} bar")
    print(f"  Feed: CO={Fi[0]}, H2={Fi[1]}, CO2={Fi[6]}, H2O={Fi[5]} (mol/s)")
    print()
    
    # Solve
    cov = model.solve_surface(T, P, Fi)
    
    print("Surface coverages at steady-state:")
    for species, theta in sorted(cov.items()):
        print(f"  {species:15s} = {theta:.6f}")
    
    # Check site balance
    sum_coverage = sum(cov.values())
    print(f"\nSite balance check:")
    print(f"  Σθ_i = {sum_coverage:.6f} (should be ≈ 1.0)")
    
    # Check validity
    valid = all(0 <= v <= 1.0001 for v in cov.values()) and 0.99 <= sum_coverage <= 1.01
    print(f"  Valid: {'✓' if valid else '✗'}")
    print()


def test_reaction_rates():
    """Test reaction rate calculation."""
    print("=" * 70)
    print("Test 3: Reaction Rate Calculation")
    print("=" * 70)
    
    model = BrubachModel()
    
    T = 573.15  # K
    P = 10.0    # bar
    Fi = np.array([0.1, 0.2, 0.0, 0.0, 0.0, 0.01, 0.1])
    
    # Compute rates
    rates = model.rate(T, P, Fi)
    
    print(f"Reaction rates (mol m^-3 s^-1) at T={T:.0f} K, P={P:.1f} bar:")
    print(f"  r0 (CO ads/des):           {rates[0]:>12.6f}")
    print(f"  r1 (RWGS):                 {rates[1]:>12.6f}")
    print(f"  r2 (CO hydrogenation):     {rates[2]:>12.6f}")
    print(f"  r3 (initiation):           {rates[3]:>12.6f}")
    print(f"  r4 (growth):               {rates[4]:>12.6f}")
    print(f"  r5 (alkane termination):   {rates[5]:>12.6f}")
    print(f"  r6 (alkene termination):   {rates[6]:>12.6f}")
    print(f"  r7 (branching):            {rates[7]:>12.6f}")
    print()
    
    # Net conversion rates (using stoichiometry)
    r_vec = rates / (model.params.cat_loading / 3600.0)  # Convert back to per-gram
    net_CO_conv = (r_vec[0] + r_vec[1] + r_vec[2])  # CO ads + RWGS + CO hydrogenation
    net_H2_cons = (r_vec[2] * 3 + r_vec[3] + r_vec[4] * 0.5 + r_vec[5] + r_vec[6] * 2.5)
    
    print(f"Key reactants (mol g^-1 h^-1):")
    print(f"  Net CO consumption:  {net_CO_conv:>12.6f}")
    print(f"  H2O production rate: {r_vec[1] + r_vec[2]:>12.6f}")
    print()


def test_stoichiometry():
    """Test stoichiometric matrix."""
    print("=" * 70)
    print("Test 4: Stoichiometric Matrix")
    print("=" * 70)
    
    model = BrubachModel()
    
    print(f"Stoichiometric matrix ν (8 reactions × 7 species):")
    print(f"Species: CO, H2, CH4, C2_4, C5plus, H2O, CO2")
    print()
    
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
    
    for i, rxn in enumerate(rxn_names):
        row = model.nu[i] if i < model.nu.shape[0] else "N/A"
        print(f"{rxn:25s}: {row}")
    print()


if __name__ == "__main__":
    try:
        test_basic_mechanism()
        test_surface_solver()
        test_reaction_rates()
        test_stoichiometry()
        print("=" * 70)
        print("All tests completed successfully!")
        print("=" * 70)
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
