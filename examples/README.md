# Examples Directory

This directory contains runnable demonstration scripts for the FT-Model-Analysis framework.

## Available Examples

### 1. Bifunctional Packed-Bed Reactor Demo (`run_bifunctional_demo.py`)

A comprehensive demonstration of the bifunctional reactor model combining RWGS, Fischer-Tropsch, and zeolite upgrading reactions.

**To run:**
```bash
python examples/run_bifunctional_demo.py
```

**What it demonstrates:**
- Loading the component registry from YAML database
- Building an IDAES flowsheet with time discretization
- Creating a simplified property package
- Setting up inlet feed conditions (CO2/H2 mixture)
- Building the bifunctional reactor model structure
- Displaying reaction network stoichiometry
- Model initialization approach
- Inlet composition analysis

**Key Features:**
- **19-component system**: RWGS feeds + FT products + zeolite products
- **1D spatial discretization**: ContinuousSet along catalyst weight W ∈ [0, 1]
- **Three coupled reaction networks**:
  - RWGS: 1 reversible reaction (CO2 + H2 ↔ CO + H2O)
  - Fischer-Tropsch: 7 product distribution reactions
  - Zeolite upgrading: 3 cracking/aromatization pathways
- **State variables**: flow_mol_comp[t, W, component], temperature[t, W], pressure[t, W]

**Inlet Conditions (Default):**
- Feed: 70% CO2 + 30% H2
- Total molar flow: 1.0 kmol/s
- Temperature: 523.15 K (250°C)
- Pressure: 20 bar (2 MPa)

**Output Example:**
```
================================================================================
BIFUNCTIONAL PACKED-BED REACTOR DEMO
================================================================================

Step 1: Loading component registry...
✓ Loaded 19 components
  Available components: CO2, H2, CO, H2O, CH4...

Step 2: Building Pyomo model and IDAES flowsheet...
✓ Created ConcreteModel and FlowsheetBlock

...

DEMO SUMMARY
================================================================================

Bifunctional Reactor Demo Completed Successfully!

Model Components:
  - 19-species component system (RWGS + FT + Zeolite products)
  - 1D spatial discretization along catalyst weight (W)
  - Three coupled reaction networks:
    • RWGS: 1 reversible reaction
    • FT: 7 product distribution reactions
    • Zeolite: 3 upgrading pathways

Inlet Conditions:
  - Feed: 70% CO2 + 30% H2
  - Flow: 1.0 kmol/s
  - Temperature: 523.1 K
  - Pressure: 20.0 bar
```

## Modifying the Examples

To customize the demo script:

1. **Change inlet composition**: Modify the `setup_feed_inlet_conditions()` function
2. **Adjust temperature/pressure**: Update values in the inlet_conditions dict
3. **Add solver call**: Uncomment solver section in Step 10 and implement IPOPT call
4. **Display additional outputs**: Extend Step 9 to report more species or states

## Next Steps

For production use, consider:

1. **Full reactor instantiation**: Use the actual `BifunctionalPackedBedReactor` class
2. **Solver integration**: Add IPOPT solver with appropriate tolerances
3. **Discretization**: Configure finite difference or collocation methods
4. **Energy balance**: Include non-isothermal operation
5. **Pressure drop**: Add momentum balance equations
6. **Flowsheet integration**: Connect reactor to downstream separation units

## Related Files

- [bifunctional_reactor.py](../src/ft_model/bifunctional_reactor.py) - Full reactor implementation
- [rwgs_reactor.py](../src/ft_model/rwgs_reactor.py) - RWGS-only reactor model
- [test_bifunctional_reactor.py](../tests/test_bifunctional_reactor.py) - Structural tests (56 tests)
- [test_bifunctional_reactor_integration.py](../tests/test_bifunctional_reactor_integration.py) - Integration tests (27 tests)
- [README.md](../README.md) - Main project documentation
