from Bio import PDB
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

def readfile_pdb(file_path):
    atom_data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('ATOM'):
                atom_info = {
                    'name': line[12:16].strip(),
                    'residue_number': int(line[22:26]),
                    'line': line,
                }
                atom_data.append(atom_info)
    return atom_data


def recalc_and_write_pdb(file_path1, file_path2, output_file_path):
    pdb_data1 = readfile_pdb(file_path1)
    pdb_data2 = readfile_pdb(file_path2)

    aligned_data = []

    for atom1 in pdb_data1:
        for atom2 in pdb_data2:
            if atom1['name'] == atom2['name'] and atom1['residue_number'] == atom2['residue_number']:
                aligned_data.append(atom1['line'])
                break
    
    renumbered_data = []
    atom_counter = 1
    for line in aligned_data:
        renumbered_line = line[:6] + f"{atom_counter:5d}" + line[11:]
        renumbered_data.append(renumbered_line)
        atom_counter += 1



    with open(output_file_path, 'w') as output_file:
        for line in renumbered_data:
            output_file.write(line)


model_files = ['R1107TS081-model1.pdb', 'R1107TS081-model2.pdb', 'R1107TS081-model3.pdb', 'R1107TS081-model4.pdb']


reference_file = 'R1107_reference.pdb'


for model_file in model_files:

    file_path1 = model_file
    file_path2 = reference_file
    output_file_path = model_file.replace('.pdb', '_rec.pdb')


    recalc_and_write_pdb(file_path1, file_path2, output_file_path)
