def calculate_clash_score(pdb_file):
    
    with open(pdb_file, 'r') as file:
        pdb_lines = file.readlines()

    
    clashes = 0
    clash_details = []

    
    for i, line in enumerate(pdb_lines):
        
        if line.startswith('ATOM') or line.startswith('HETATM'):
            
            atom_type = line[13:16].strip()
            x_coord = float(line[30:38].strip())
            y_coord = float(line[38:46].strip())
            z_coord = float(line[46:54].strip())
            residue_num = int(line[22:26].strip())  

            
            atom_radii = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'H': 1.2, 'P': 1.8}  
            atom_radius = atom_radii.get(atom_type, 1.52)  
            for key in atom_radii:
                if key in atom_type:
                    atom_radius = atom_radii[key]
                    break
            
            for j, other_line in enumerate(pdb_lines):
                if j != i and (other_line.startswith('ATOM') or other_line.startswith('HETATM')):
                    
                    other_x = float(other_line[30:38].strip())
                    other_y = float(other_line[38:46].strip())
                    other_z = float(other_line[46:54].strip())
                    other_residue_num = int(other_line[22:26].strip())

                    
                    # comparing only residues that are not neighbours of the residue (so that no atoms have bonds with each other and thus are calculated as a clash)
                    if residue_num + 1 < other_residue_num:
                        distance = ((x_coord - other_x)**2 + (y_coord - other_y)**2 + (z_coord - other_z)**2)**0.5
                        other_atom_type = other_line[13:16].strip()
                        other_atom_radius = atom_radii.get(other_atom_type, 1.55)
                        for key in atom_radii:
                            if key in other_atom_type:
                                other_atom_radius = atom_radii[key]
                                break
                        

                        overlap = atom_radius + other_atom_radius-distance
                        #x = ((atom_radius + other_atom_radius) * 2 - 0.4)
                        if overlap >=  0.4:
                            clashes += 1
                            #print(f"{atom_type} {i} :{overlap:.3f} with {other_atom_type} {j}")
                            clash_details.append((i, j, overlap,distance))

   


    for clash in clash_details:
        i, j,overlap,distance = clash
        line_i = pdb_lines[i].strip()
        line_j = pdb_lines[j].strip()
        print(f"{line_i} :{overlap:.3f},{distance:.3f} with {line_j}")

    #print(f"\nclashscore = {clashes / (len(pdb_lines) / 1000):.2f}")


    return clashes / (len(pdb_lines) / 1000)



pdb_file_path = 'R1107_reference.pdb'
#pdb_file_path = 'R1107TS081-model4_rec.pdb'

score = calculate_clash_score(pdb_file_path)
print(f'Clash Score: {score}')
