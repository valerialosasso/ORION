python3 ../run_pipeline.py --pdb colicin_h.pdb \
  --theta-angles 0:90:15 --phi-angles 0:90:15 \
  --outdir orientations_H --box-size 150 --nstruc 5 \
  --bb-window 0 \
  --split-hydrogens --sasa-cutoff 0.0 --probe-radius 1.4 --sasa-points 960
