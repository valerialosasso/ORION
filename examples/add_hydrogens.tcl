# 1. Load the psfgen package
package require psfgen

# 2. Read the topology file (make sure you point to a valid CHARMM topology on your system)
# VMD usually includes these in its installation directory
topology top_all27_prot_lipid.inp

# 3. Create a segment from your original PDB file (replace 'colicin.pdb' with your file name)
segment U {pdb colicin.pdb}

# 4. Read the coordinates from your PDB into the segment
coordpdb colicin.pdb U

# 5. Tell psfgen to calculate coordinates for all missing atoms (including hydrogens)
guesscoord

# 6. Write out the new files with hydrogens included
writepdb protein_with_hydrogens.pdb
writepsf protein_with_hydrogens.psf
