from Bio import PDB
from Bio.PDB import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

def readfile_pdb(file_path):

    structure = PDB.PDBParser().get_structure('structure', file_path)
    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atoms.append(atom)
    return atoms

def atoms_to_coordinates(atoms):

    return np.array([atom.coord for atom in atoms])


file_paths = [
    "R1107_reference.pdb",
    "R1107TS081-model1_rec.pdb",
    "R1107TS081-model2_rec.pdb",
    "R1107TS081-model3_rec.pdb",
    "R1107TS081-model4_rec.pdb"
]

reference_atoms = readfile_pdb(file_paths[0])


for file_path in file_paths[1:]:
    model_atoms = readfile_pdb(file_path)


    common_atoms = set(reference_atoms) & set(model_atoms)

    ref_coords = atoms_to_coordinates([atom for atom in reference_atoms if atom in common_atoms])
    mod_coords = atoms_to_coordinates([atom for atom in model_atoms if atom in common_atoms])

    #superimposer object - svd -  singular value decomposition
    super_imposer = SVDSuperimposer()
    #setting common positions in superimpose objects
    super_imposer.set(ref_coords, mod_coords)


    super_imposer.run()
    #tran - translation vector is a three-dimensional vector that describes the displacement of
    #the geometry of one structure relative to another in space
    #rot - rotation matrix describes the rotational transformation applied to the coordinates of
    #one structure to align it with another structure

    rot, tran = super_imposer.get_rotran()


    structure = PDB.PDBParser().get_structure('structure', file_path)
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for i, atom in enumerate(residue):
                    atom.transform(rot, tran)
    
    output_file_path = file_path.replace('.pdb', '_sup.pdb')
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_file_path)


