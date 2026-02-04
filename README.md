# FT-Model-Analysis: IDAES-based Process Modelling Framework

## Project Purpose

This project implements mechanistic Fischer-Tropsch (FT) catalyst modelling combined with lumped zeolite post-processing model. The framework models hydrocarbon synthesis and subsequent conversion using state-of-the-art process modelling tools. Includes custom unit models for reverse water-gas shift (RWGS) reactors with spatial discretization.

## Tools & Technologies

- **Python 3.8+**: Core programming language
- **IDAES-PSE**: Integrated DAE Environment and Solver with Process Systems Engineering models
- **Pyomo**: Python-based optimization modelling language
- **NumPy/SciPy**: Numerical computing and scientific algorithms
- **Pandas**: Data manipulation and analysis
- **Matplotlib**: Visualization

## Project Structure

```
src/
├── components/         # Component database with thermodynamic properties
│   ├── component.py   # Component class definition
│   ├── registry.py    # ComponentRegistry for managing component definitions
│   └── components.yaml # YAML definitions for 39 chemical species
├── ft_model/           # Fischer-Tropsch and specialized reactor models
│   ├── rwgs_reactor.py # RWGS packed-bed reactor with 1D spatial discretization
│   ├── bifunctional_reactor.py  # Bifunctional reactor combining RWGS + FT + zeolite
│   ├── kinetics.py    # FT rate expressions and reaction networks
│   └── reactor.py     # IDAES-based FT reactor unit models
├── zeolite_model/      # Zeolite conversion kinetics and models
│   ├── kinetics.py    # Zeolite reaction kinetics
│   └── reactor.py     # IDAES-based zeolite reactor unit models
├── flowsheet/          # Process flowsheet assembly
│   └── flowsheet.py   # Complete process flowsheet construction
└── utils/             # Utility modules
    ├── lumping.py     # FT-to-zeolite product lumping and mapping
    └── parameters.py  # Parameter loading and configuration management
```

## Installation

1. Create and activate the virtual environment:
   ```bash
   python -m venv .venv
   .venv\Scripts\activate  # Windows
   # or
   source .venv/bin/activate  # Linux/Mac
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Configure environment variables:
   - Copy `.env` and update paths as needed:
   ```bash
   IDAES_CONFIG_PATH="/path/to/idaes/config"
   DATA_PATH="./data"
   ```

## Key Features

### RWGS Reactor Model
- **Reverse Water-Gas Shift (RWGS) Packed-Bed Reactor**: Custom IDAES unit model with full 1D spatial discretization
  - Reaction: CO2 + H2 ↔ CO + H2O (endothermic, equilibrium-limited)
  - Spatial domain: Catalyst weight (W) discretization using Pyomo.DAE `ContinuousSet`
  - Material balance equations: Differential-algebraic with `DerivativeVar` formulation
  - Rate expression: Forward/reverse kinetics with equilibrium constant dependency
  - Partial pressure calculations from mole fractions
  - Isothermal operation with optional temperature constraint
  - Performance reporting with inlet/outlet conditions

### Bifunctional Reactor Model
- **Combined RWGS + Fischer-Tropsch + Zeolite Upgrading**: Advanced packed-bed reactor with coupled multi-stage reactions
  - **Stage 1 - RWGS**: CO2 + H2 ↔ CO + H2O (produces syngas for FT)
  - **Stage 2 - Fischer-Tropsch Synthesis**: CO + H2 → CH4, C2H4, C2H6, C3H6, C3H8, C5+, wax
  - **Stage 3 - Zeolite Upgrading**: Wax → Distillate → Naphtha, Olefins → Aromatics → Coke
  - 1D spatial discretization along catalyst bed (catalyst weight dimension)
  - 19-component system including RWGS feeds, FT products, and zeolite products
  - Coupled material balance equations with contributions from all three reaction networks
  - Pressure-dependent FT kinetics and first-order zeolite cracking rates
  - Partial pressure calculations and isothermal operation
  - Full performance reporting with inlet/outlet product yields

### Component Database System
- Comprehensive database of 39+ chemical species (CO2, H2, CO, H2O, CH4, C2H4, etc.)
- YAML-based component definitions with thermodynamic properties
- Component registry for managing definitions across models
- Full test coverage (39 component definition tests)

## Usage

### Basic Workflow

```python
from src.ft_model import kinetics, reactor, rwgs_reactor, bifunctional_reactor
from src.zeolite_model import kinetics as zeo_kinetics
from src.flowsheet import flowsheet
from src.utils import lumping, parameters
from src.components import registry

# Load component database
comp_registry = registry.ComponentRegistry()
comp_registry.load_from_yaml('src/components/components.yaml')

# Load parameters from environment
params = parameters.load_parameters()

# Build RWGS reactor (single reaction system)
rwgs = rwgs_reactor.PackedBedRWGSReactorData()

# Build bifunctional reactor (RWGS + FT + zeolite)
inlet_conditions = {
    'CO2': 1.0,
    'H2': 2.0,
    'CO': 0.1,
    'H2O': 0.05,
    # ... other components
}
bifunctional = bifunctional_reactor.build_bifunctional_reactor(
    flowsheet=None,  # parent flowsheet
    property_package=None,  # property package
    inlet_flow=inlet_conditions,
    inlet_temperature=523.15,  # 250°C
    inlet_pressure=30 * 101325.0,  # 30 bar
    W_total=5.0  # 5 kg catalyst
)

# Build FT reactor
ft_reactor = reactor.build_ft_reactor(params)

# Build zeolite reactor
zeo_reactor = zeolite_reactor.build_zeolite_reactor(params)

# Build complete flowsheet
fs = flowsheet.build_flowsheet(rwgs, bifunctional, ft_reactor, zeo_reactor)
```

## Testing

Run all tests with pytest:
```bash
pytest tests/          # Run all 161 tests
pytest tests/test_rwgs_reactor.py -v  # Run RWGS reactor tests (23 tests)
pytest tests/test_bifunctional_reactor.py -v  # Run bifunctional reactor structural tests (56 tests)
pytest tests/test_bifunctional_reactor_integration.py -v  # Run bifunctional integration tests (27 tests)
```

Test coverage includes:
- Component database validation (39 tests)
- RWGS reactor physics and structure (23 tests)
- Bifunctional reactor structural tests (56 tests)
- Bifunctional reactor integration tests (27 tests):
  * Spatial domain creation and bounds
  * State variable initialization
  * Partial pressure calculations
  * Total flow balance constraints
  * RWGS, FT, and zeolite rate expressions
  * Material balance with DerivativeVar
  * Inlet boundary conditions
  * Kinetic parameter settings
  * Non-negative flow constraints
  * Temperature and pressure bounds
  * Isothermal and constant pressure modes
  * Helper function availability
- Stoichiometric balance validation
- Reaction kinetics verification
- Material balance equations with DerivativeVar
- Partial pressure calculations
- Coupled reaction networks
- Zeolite cracking stoichiometry
- Legacy tests (16 tests)

## License

[Add license information]

## Authors

[Add author information]
