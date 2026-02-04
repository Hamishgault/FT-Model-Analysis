# FT-Model-Analysis: IDAES-based Process Modelling Framework

## Project Purpose

This project implements mechanistic Fischer-Tropsch (FT) catalyst modelling combined with lumped zeolite post-processing model. The framework models hydrocarbon synthesis and subsequent conversion using state-of-the-art process modelling tools.

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
├── ft_model/           # Fischer-Tropsch reactor kinetics and models
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

## Usage

### Basic Workflow

```python
from src.ft_model import kinetics, reactor
from src.zeolite_model import kinetics as zeo_kinetics
from src.flowsheet import flowsheet
from src.utils import lumping, parameters

# Load parameters from environment
params = parameters.load_parameters()

# Build FT reactor
ft_reactor = reactor.build_ft_reactor(params)

# Build zeolite reactor
zeo_reactor = zeolite_reactor.build_zeolite_reactor(params)

# Build complete flowsheet
fs = flowsheet.build_flowsheet(ft_reactor, zeo_reactor)
```

## Testing

Run tests with pytest:
```bash
pytest tests/
```

## License

[Add license information]

## Authors

[Add author information]
