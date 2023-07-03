import pymol
import argparse
import os
import requests

parser = argparse.ArgumentParser(
                    prog='PyMOL surface D/E residue identification',
                    description='For a given list of pdb files, generate the surface D/E residues. This script uses the `get_sasa_relative` command from PyMOL with a default cutoff of 0.2',
                    )
parser.add_argument('-l', '--list', type=str, required=True, help='Path to the input csv file with `UniProtID, file path, chain ID` ')
parser.add_argument('-d', '--directory', type=str, required=True, help='Path to the output directory.')
parser.add_argument('-c', '--cutoff', type=float, default=0.2, help='cutoff for relative per residue solvent accessible area' )

args = parser.parse_args()
pdb_list = args.list
output_path = args.directory
cutoff = args.cutoff

def fetch_pdb_from_github(file_path):
    raw_url = f"https://raw.githubusercontent.com/{file_path}"
    response = requests.get(raw_url)
    if response.status_code == requests.codes.ok:
        pdb_content = response.text
        return pdb_content
    else:
        print("Failed to fetch PDB file:", response.status_code)
        return None
     

def pymol_surf_DE(pdb_file, chain_id, cutoff, protein_id):
    pymol.cmd.load(pdb_file, 'structure')
    sasa_dict = pymol.cmd.get_sasa_relative()
    rel_sa_cutoff = cutoff
    chain_id = 'A'
    surface_residues = [int(x[-1]) for x, val in sasa_dict.items() if val > rel_sa_cutoff and x[2] == chain_id]

    # Dictionary to map three-letter amino acid codes to single-letter codes
    aa_dict = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    }

    # Containers for residues and indices
    residues = []
    indices = []

    # Iterate over each residue in the structure
    cmd.iterate('structure and chain {} and name CA'.format(chain_id), 'residues.append(resn); indices.append(resi)', space=locals())

    # Convert three-letter codes to single-letter codes
    residues = [aa_dict[res] for res in residues]

    # Convert indices to numbers
    indices = [int(x) for x in indices]
    
    # Join residues into a string
    prot_seq = ''.join(residues)

    # look up the indices in `indices` to check for residue identity to D/E
    surf_DE_list = [int(x) for x in surface_residues if prot_seq[indices.index(x)] in ['D', 'E']]
    with open(os.path.join(output_path, '{}_SurfDE.txt'.format(protein_id)), 'w') as f:
        f.write(','.join(map(str, surf_DE_list)))

pymol.finish_launching(['pymol', '-qc']) 

with open(pdb_list, 'r') as f:
    for line in f:
        print(line)
        prot_id, pdb_file_path, chain_id = line.split(',')
        prot_id, pdb_file_path, chain_id = prot_id.strip(), pdb_file_path.strip(), chain_id.strip()
        pdb_file_content = fetch_pdb_from_github(pdb_file_path)
        pdb_file_name = './{}.pdb'.format(prot_id)
        with open(pdb_file_name,'w') as f:
            f.write(pdb_file_content)
        pymol_surf_DE(pdb_file_name, chain_id, cutoff, prot_id)
        os.remove(pdb_file_name)

pymol.cmd.quit()
    

