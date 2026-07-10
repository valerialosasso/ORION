# ORION
Fast computational workflow that generates orientation-dependent volume fraction and SLD profiles directly from atomic structures without requiring full MD simulations.
This workflow combines small stochastic torsional perturbations (“wriggling”) with systematic sampling of protein orientation. Profiles are calculated from van der Waals representations and converted into SLD profiles using Gaussian smoothing. The resulting profiles are then assembled into orientation lookup tables.
The wriggling procedure produces smoother and more experimentally realistic profiles than purely rigid-body models, while preserving the overall protein structure.
The resulting profiles capture the expected orientation dependence observed in neutron reflectometry experiments and provide a practical tool for rapid orientational screening.

