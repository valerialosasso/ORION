# ORION
Fast computational workflow that generates orientation-dependent volume fraction and SLD profiles directly from atomic structures without requiring full MD simulations.
This workflow combines small stochastic torsional perturbations (“wriggling”) with systematic sampling of protein orientation. Profiles are calculated from van der Waals representations and converted into SLD profiles using Gaussian smoothing. The resulting profiles are then assembled into orientation lookup tables.
The wriggling procedure produces smoother and more experimentally realistic profiles than purely rigid-body models, while preserving the overall protein structure.
The resulting profiles capture the expected orientation dependence observed in neutron reflectometry experiments and provide a practical tool for rapid orientational screening.


## Optional hydrogen accessibility split

Install the optional dependency with `python -m pip install -r requirements-hydrogen.txt`.
MDTraj supplies the published Shrake–Rupley SASA implementation; ORION supplies
the covalent-parent inference and binary grouping. No PSF, force-field topology,
or MD trajectory is required. An explicit-hydrogen, solute-only PDB is required
(the bundled Colicin example lacks hydrogens). Add hydrogens/protonate your
structure at the intended conditions before running ORION. Resolve alternate
locations and use a single model. Do not include explicit water.

```bash
python run_pipeline.py --pdb protein_with_H.pdb \
  --theta-angles 0:90:15 --phi-angles 0:90:15 \
  --outdir orientations_H --box-size 150 --nstruc 20 \
  --split-hydrogens --sasa-cutoff 0.0 --probe-radius 1.4 --sasa-points 960
```

The Python `main_api` accepts `split_hydrogens=True`, `sasa_cutoff_A2=0.0`,
`probe_radius_A=1.4`, and `sasa_points=960`.

- `density_H_group1.dat`: C-bound H plus buried N/O/S-bound H.
- `density_H_group2.dat`: N/O/S-bound H with SASA strictly above the cutoff.
- `density_H.dat`: total H, retained for existing consumers.
- Lookup JSON gains `nH_group1` and `nH_group2`.
- `hydrogen_classification.json`: per-conformation atom names/indices, inferred
  parents, H SASA in A^2, and group assignments. Run metadata records parameters.

Density files contain z (A), mean number density (A^-3), and sample standard
deviation across conformations (zero for a single conformation). Their mean
profiles satisfy `nH_group1 + nH_group2 = nH` bin by bin. Standard deviations
are not additive.

Parent inference requires exactly one heavy atom within a covalent-distance
window: lower bound 0.60 A and upper bounds C 1.25, N 1.20, O 1.15, S 1.50 A.
These are conservative geometric heuristics, not force-field bond assignments.
Missing or ambiguous parents raise an error rather than silently assigning a
group. Inspect the audit for unusual/proton-transferring structures. All H in
the supplied solute are classified; for protein-only profiles supply protein
coordinates only. All supplied heavy atoms, including sulfur, occlude SASA.

SASA uses MDTraj elemental radii and a 1.4 A probe by default. The default zero
cutoff means any sampled accessible surface, **not a calibrated exchange
threshold**. Check sensitivity to the cutoff and sphere-point resolution when
using the results quantitatively. Classification is computed separately for
each wriggled conformation, once before rigid rotation, and reused for all its
orientations. This avoids orientation-dependent surface-sampling noise. No
membrane/substrate shielding is inferred if those components are absent. Local
SASA does not establish connectivity of internal cavities to bulk solvent.

The split models solvent exposure, not H/D-exchange kinetics or history.
Buried labile H remain chemically exchangeable; conformational opening can
allow exchange over time. The two densities can be used to apply separate H/D
scattering weights downstream. This option does not itself change the existing
protiated `sld_total` or introduce a D2O fraction. Existing total-SLD code still
sums C/H/N/O only; it is not a complete sulfur/phosphorus SLD calculation.

The density counter now counts individual atom centres rather than occupied
voxels, including atoms on grid planes. Consequently regenerated densities may
differ from earlier ORION outputs even without hydrogen splitting. Volume
voxelisation is unchanged. With no specified box, ensemble profiles use common
bounds so corresponding z bins align.

Validation: `python -m unittest discover -s tests -p 'test_*.py' -v`.

References:
- McGibbon et al. (2015), MDTraj, *Biophysical Journal* 109, 1528–1532,
  https://doi.org/10.1016/j.bpj.2015.08.015
- https://mdtraj.readthedocs.io/en/stable/api/generated/mdtraj.shrake_rupley.html
