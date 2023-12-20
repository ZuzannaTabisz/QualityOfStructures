# QualityOfStructures

clashscore.py:

Calculates clash scores for a given PDB file based on atomic clashes.
Identifies clashes between atoms and provides details such as atom types, overlap, and distances.

rmsd.py:

Computes the Root Mean Square Deviation (RMSD) between protein structures.
Superimposes structures using the SVDSuperimposer from Biopython.
Outputs RMSD values for comparing a reference structure to models.

recalc.py:

Recalculates and writes a new PDB file with aligned atom numbering based on a reference structure.
Outputs the recalculated and renumbered strucure named original name + _rec.pdb

superimpose.py:

Superimposes protein structures onto a reference structure.
Outputs transformed structures with aligned coordinates named original name + _sup.pdb


Requirements:
Biopython library (install using pip install biopython)
