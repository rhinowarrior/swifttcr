"""
For now this works with the 2iam.pdb file. It will identify the chains in the PDB file and replace them with A, B, C, D, E. But not sure if this is possible cause i don't know how other PDB files are structured that can be used for this because it stops working if their are more or less than 5 chains in the PDB file. (and i'm not sure if the chains are always in the same order in the PDB file) (pandora has a method but I do not know how to implement it for my own use)
"""
from Bio import PDB

def get_chains(pdb_file):
    parser = PDB.PDBParser()
    structure = parser.get_structure('structure', pdb_file)
    chains = [chain.id for chain in structure.get_chains()]
    return chains

def replace_chains(pdb_data, replacements):
    modified_lines = []
    for line in pdb_data.splitlines():
        if line.startswith('ATOM'):
            chain_id = line[21]
            for old_chain, new_chain in replacements:
                if chain_id == old_chain:
                    line = line[:21] + new_chain + line[22:]
                    break
        modified_lines.append(line)
    return '\n'.join(modified_lines)

if __name__ == "__main__":
    pdb_file = '/home/nils/swifttcr/example/input/2iam.pdb'
    
    chains = get_chains(pdb_file)
    
    # Read the input PDB file
    with open(pdb_file, 'r') as file:
        pdb_data = file.read()

    # Define the chain replacements
    replacements = [
        (str(chains[0]), 'A'),
        (str(chains[1]), 'B'),
        (str(chains[2]), 'C'),
        (str(chains[3]), 'D'),
        (str(chains[4]), 'E')
    ]

    # Perform the chain replacements
    new_pdb_data = replace_chains(pdb_data, replacements)

    # Write the modified PDB data to a new file
    with open('modified_file.pdb', 'w') as file:
        file.write(new_pdb_data)

    print('Chain replacements complete.')