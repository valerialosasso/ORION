python3 ../run_pipeline.py --pdb colicin_h.pdb \
  --theta-angles 0:90:30 --phi-angles 0:180:30 \
  --outdir test_rawoutput --box-size 120 --nstruc 2 \
  --bb-window 0 \
  --split-hydrogens --sasa-cutoff 0.0 --probe-radius 1.4 --sasa-points 960
