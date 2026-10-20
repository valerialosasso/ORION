python3 ../run_pipeline.py --pdb colicin_h.pdb \
  --theta-angles 0:90:30 --phi-angles 0:180:30 \
  --outdir orientations_Colicin_split_short --box-size 120 --nstruc 2 \
  --bb-window 0 \
  --split-hydrogens --sasa-cutoff 0.0 --probe-radius 1.4 --sasa-points 960
