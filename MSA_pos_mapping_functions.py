import requests as r

def get_sequence(uniprot_id):
    """get sequence from UniProt for a given UniProt ID, return as a string"""
    #print(uniprot_id)
    baseUrl="http://www.uniprot.org/uniprot/"
    currentUrl=baseUrl+uniprot_id+".fasta"
    response = r.post(currentUrl)
    seq = response.text
    seq = ''.join(seq.split('\n')[1:])
    
    return seq

def msa_pos_mapping(full_seq, sub_seq, cutoff=4):
    """
    This function takes in the full protein sequence `full_seq` and the sequence as given in
    the MSA entry `sub_seq`, and returns the positions of all residues appeared in the full
    length sequence.
    """
    sub_seq = sub_seq.replace("-", "").replace("\n","").upper()
    full_seq = full_seq.replace("\n","").upper()

    full_seq_pos = 0  # Initialize position in full sequence (1-based)
    #print(full_seq_pos)
    sub_seq_pos = 0  # Initialize position in subsequence (1-based)
    mappings = []  # List to hold position mappings
    
    # set a cutoff to require at least certain number of continuous residues. 
    continuity_cutoff = 4
    while full_seq_pos < len(full_seq) and sub_seq_pos < len(sub_seq):
        # find initial matches
        while  (full_seq_pos < len(full_seq) and sub_seq[sub_seq_pos] != full_seq[full_seq_pos]):
            #print(full_seq_pos)
            full_seq_pos += 1
            continue
        # mappings.append(full_seq_pos)
        tentative_mappings = [] # List to hold the most recent continuous pos mappings.
        #tentative_mappings.append(full_seq_pos)
        while (
            full_seq_pos < len(full_seq)
            and sub_seq_pos < len(sub_seq)
            and sub_seq[sub_seq_pos] == full_seq[full_seq_pos]
        ):
            tentative_mappings.append(full_seq_pos)
            sub_seq_pos += 1
            full_seq_pos += 1
        if len(tentative_mappings) >= continuity_cutoff:
            for res_id in tentative_mappings:
                mappings.append(res_id)
    mappings = [x+1 for x in mappings]
    mappings_string = '\n'.join(str(element) for element in mappings)
    return mappings_string

def get_MSA_seq(msa_file, target_id):
    """
    find the MSA entry for a given ID, e.g. `>467818` and return it
    """
    #print(target_id)
    with open(msa_file, 'r') as file:
        lines = file.readlines()
        for i in range(len(lines)):
            if target_id in lines[i]:
                return lines[i + 1]
