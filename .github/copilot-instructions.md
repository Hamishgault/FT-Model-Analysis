# Copilot instructions (FT-Model-Analysis)

## Architecture (what matters first)
- Core model is the IDAES unit model `FTRWGSReactor` in src/ft_model/ft_rwgs_zeolite_reactor.py. It owns all stoichiometry, parameters, and constraints for RWGS + FT + optional zeolite upgrading.
- The reactor is a 1D packed-bed DAE: `ContinuousSet W ∈ [0,1]` with `DerivativeVar` on `flow_mol_comp` (dF/dW). Discretize with `discretize_reactor()` before solving.
- Stoichiometry is validated on module import via `check_atom_balance()` and `ATOMIC_COMPOSITION`; `verify_atom_conservation()` is the post-solve check.
- Feature toggles are config flags on `FTRWGSReactor` (e.g., `include_zeolite_reactions`, `energy_balance`, `pressure_drop`, `ergun_pressure_drop`, `heat_transfer`, `mass_transfer`). Kinetics switch on `kinetics_model` (`lumped_simple`, `marvast_2005`, `rwgs_2017`).
- The component system uses YAML in src/components/data loaded by `ComponentRegistry` (src/components/registry.py). This is the source of truth for species definitions in demos.

## Developer workflows (actual, not aspirational)
- Self-test run: `python src/ft_model/ft_rwgs_zeolite_reactor.py` (prints stoichiometry verification and atom balance).
- CO2 feed sweep: `python examples/compare_co2_feed_analysis.py` (uses staged solve with `rate_multiplier` and `pressure_drop_multiplier`).
- Demo flowsheet: `python examples/run_bifunctional_demo.py` (builds a simplified structure; no full solve).
- Dependencies: IDAES + Pyomo + IPOPT (see requirements.txt and README.md for setup/env vars).

## Project-specific patterns to follow
- Initialization fixes inlet at `W=0`, then sets profiles; when `energy_balance`/`pressure_drop` are off, temperature/pressure are fixed along W (see `FTRWGSReactorData.initialize`).
- Keep new reaction stoichiometries atom-balanced by updating `ATOMIC_COMPOSITION` and the relevant reaction dictionaries in src/ft_model/ft_rwgs_zeolite_reactor.py.
- When adding constraints/vars, follow the existing naming scheme: `flow_mol_comp`, `flow_mol_total`, `partial_pressure`, `rate_*`, `reaction_rate`, `dF_dW`.

## Known gaps / consistency notes
- tests/ reference modules that do not exist in this repo (e.g., src/ft_model/bifunctional_reactor.py, rwgs_reactor.py, kinetics.py). Expect tests to fail unless updated to the current model file.
- src/flowsheet/flowsheet.py, src/zeolite_model/kinetics.py, src/zeolite_model/reactor.py, and src/utils/lumping.py are stubs (`pass`). Avoid assuming they are implemented without adding code.
