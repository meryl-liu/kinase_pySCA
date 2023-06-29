import pymol
import argparse
import os
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

def pymol_surf_DE(pdb_file, chain_id, cutoff, protein_id):
    pymol.cmd.load(pdb_file, 'structure')
    sasa_dict = pymol.cmd.get_sasa_relative()
    rel_sa_cutoff = 0.2
    chain_id = 'A'
    surface_residues = [int(x[-1]) for x, val in sasa_dict.items() if val > rel_sa_cutoff and x[2] == chain_id]
    prot_seq = pymol.cmd.get_fastastr('structure and chain {}'.format(chain_id))
    # prot_seq = pymol.cmd.get_fastastr(key='model_{}'.format(chain_id))
    prot_seq = ''.join(prot_seq.split('\n')[1:])
    print(prot_seq)
    surf_DE_list = [int(x) for x in surface_residues if prot_seq[x-1] in ['D', 'E']]
    with open(os.path.join(output_path, '{}_SurfDE.txt'.format(protein_id)), 'w') as f:
        f.write(','.join(map(str, surf_DE_list)))

pymol.finish_launching(['pymol', '-qc']) 

with open(pdb_list, 'r') as f:
    for line in f:
        prot_id, pdb_file, chain_id = line.split(',')
        prot_id, pdb_file, chain_id = prot_id.strip(), pdb_file.strip(), chain_id.strip()
        pymol_surf_DE(pdb_file, chain_id, cutoff, prot_id)

pymol.cmd.quit()
    

