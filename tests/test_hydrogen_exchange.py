import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np
from hydrogen_exchange import classify_hydrogens, infer_hydrogen_parents
from density_volume_profiles import count_cubes, calcAll
from io_utils import write_xyz_trajectory
from run_pipeline import main_api


class HydrogenTests(unittest.TestCase):
    def test_chemistry_and_accessibility(self):
        atoms = [('C', -6, 0, 0), ('H', -5, 0, 0),
                 ('N', 0, 0, 0), ('H', 1, 0, 0),
                 ('O', 6, 0, 0), ('H', 7, 0, 0),
                 ('S', 12, 0, 0), ('H', 13.34, 0, 0)]
        mask, records = classify_hydrogens(atoms)
        np.testing.assert_array_equal(mask, [False, False, False, True, False, True, False, True])
        self.assertGreater(records[0]['sasa_A2'], 0)
        mask, _ = classify_hydrogens(atoms, sasa_cutoff_A2=1000)
        self.assertFalse(mask.any())
        buried = [('N', -1, 0, 0), ('H', 0, 0, 0)]
        for axis in range(3):
            for sign in (-1, 1):
                xyz = [0., 0., 0.]
                xyz[axis] = sign * 2.5
                buried.append(('C', *xyz))
        mask, records = classify_hydrogens(buried)
        self.assertFalse(mask.any())
        self.assertEqual(records[0]['sasa_A2'], 0)

    def test_bad_geometry_fails(self):
        for atoms in ([('C', 0, 0, 0)], [('H', 0, 0, 0), ('N', 3, 0, 0)],
                      [('H', 0, 0, 0), ('N', 1, 0, 0), ('C', -1, 0, 0)]):
            with self.assertRaises(ValueError):
                infer_hydrogen_parents(atoms)

    def test_atoms_not_occupied_cells(self):
        atoms = [('H', 0, 0, 0, 1.2), ('H', .1, .1, .1, 1.2), ('H', 2, 2, 2, 1.2)]
        _, counts, _ = count_cubes(atoms, ((-2, 2),)*3, 1)
        self.assertEqual(sum(c['H'] for c in counts), 3)

    def test_partition_common_grid_and_single_frame_std(self):
        frames = [[('N', 0, 0, 0), ('H', 1, 0, 0)],
                  [('N', 0, 0, 1), ('H', 1, 0, 1)]]
        with tempfile.TemporaryDirectory() as folder:
            prev = os.getcwd()
            try:
                os.chdir(folder)
                write_xyz_trajectory('tempFiles/trajectory.xyz', frames)
                with contextlib.redirect_stdout(io.StringIO()):
                    calcAll(None, .5, [[False, True], [False, False]])
                total = np.loadtxt('density_H.dat')
                g1 = np.loadtxt('density_H_group1.dat')
                g2 = np.loadtxt('density_H_group2.dat')
                np.testing.assert_allclose(g1[:,1]+g2[:,1], total[:,1])
                self.assertTrue(np.isfinite(total).all())
                self.assertGreater(g1[:,1].sum(), 0)
                self.assertGreater(g2[:,1].sum(), 0)
            finally:
                os.chdir(prev)

    def test_pipeline_and_rotation_reuse(self):
        atom_data = [('N','N',0,0,0), ('H','H',-1,0,0), ('CA','C',1.45,0,0),
                     ('HA2','H',1.45,1.05,0), ('C','C',2.1,-1.35,0),
                     ('O','O',1.6,-2.4,0)]
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            pdb = root/'input.pdb'
            pdb.write_text(''.join(
                f'ATOM  {i:5d} {name:^4s} GLY A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n'
                for i,(name,el,x,y,z) in enumerate(atom_data,1))+'END\n')
            structures=root/'out'/'wriggles_base'/'structures'
            structures.mkdir(parents=True)
            (structures/'conf_001.pdb').write_text(pdb.read_text())
            prev = os.getcwd()
            try:
                os.chdir(root)
                with contextlib.redirect_stdout(io.StringIO()):
                    result=main_api(pdb, outdir=str(root/'out'), skip_wriggle=True,
                                    theta_angles='0:90:90', phi_angles='0:0:1',
                                    box_size=20, bin_size=1, split_hydrogens=True)
                for row in result:
                    entry=row[0]
                    np.testing.assert_allclose(np.array(entry['nH_group1'])+entry['nH_group2'], entry['nH'])
                    self.assertAlmostEqual(sum(entry['nH'])*400, 2)
                    self.assertAlmostEqual(sum(entry['nH_group2'])*400, 1)
                audit=json.loads((root/'out'/'hydrogen_classification.json').read_text())
                self.assertEqual(len(audit['conformations']),1)
                self.assertTrue(np.isfinite(np.loadtxt(root/'out'/'theta0_phi0'/'density_H_group1.dat')).all())
            finally:
                os.chdir(prev)


if __name__ == '__main__':
    unittest.main()
