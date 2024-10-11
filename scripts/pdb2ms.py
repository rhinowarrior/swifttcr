"""Script to convert pdb files to ms files with selected attractive residues.
Author: Jan Aarts

ms file format
    COLUMNS        DATA  TYPE    FIELD        DEFINITION
     1 -  6        Record name   "ATOM  "
     7 - 11        Integer       serial       Atom  serial number.
    13 - 16        Atom          name         Atom name.
    17             Character     altLoc       Alternate location indicator.
    18 - 20        Residue name  resName      Residue name.
    21             Residue suffix             Libmol extension adding a suffix to residue name.
    22             Character     chainID      Chain identifier.
    23 - 26        Integer       resSeq       Residue sequence number.
    27             AChar         iCode        Code for insertion of residues.
    31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)     attraction   Attraction
    61 - 66        Real(6.2)     surface      Surface accessibility
"""
from pathlib import Path


def open_file(filename):
    """Open a file and read the lines into a list.

    Args:
        filename (str): The path to the file to open.
        
    Returns:
        list: A list of lines from the file.
    """
    f = open(filename)
    lines = f.readlines()
    f.close()
    return lines


def pdb2ms(pdb):
    """Convert a pdb file to an ms file.
    
    Args:
        pdb (list): A list of lines from a pdb file.
    
    Returns:
        list: A list of lines for the ms file.
    """
    ms = []
    for line in pdb:
        ms_line = line[:54]
        attraction = "   0.0"
        ms_line = ms_line + attraction
        ms.append(ms_line)
    return ms


def add_attraction_tcr(ms, attractive_res):
    """Add attraction to the residues in the attractive_res dictionary.
    
    Args:
        ms (list): A list of lines from an ms file.
        attractive_res (dict): A dictionary of attractive residues.
    
    Returns:
        list: A list of lines for the ms file with attractions added.
    """
    new_ms = []
    for line in ms:
        if line.startswith("ATOM"):
            if line[21] in attractive_res.keys():
                resID = int(line[22:26])
                resdict = attractive_res[line[21]]
                start = resdict['start']
                end = resdict['end']
                hit = False
                
                for (s, e) in zip(start,end):
                    if resID >= s and resID <= e:
                        new_line = line[:54] + "   1.0"
                        hit = True
                        new_ms.append(new_line)
                if hit == False:
                    new_ms.append(line)
            else:
                new_ms.append(line)
        else:
            new_ms.append(line)
    return new_ms


def write_ms(ms, filename):
    """Write the ms lines to a file.
    
    Args:
        ms (list): A list of lines for the ms file.
        filename (str): The path to the file to write.
    """
    f = open(filename, "w")
    for line in ms:
        f.write(line + "\n")
    f.close()


def isPrepared(modelname):
    """Checks if the model is prepared.
    
    Args:
        modelname (str): The name of the model.
        
    Returns:
        bool: True if the model is prepared, False otherwise.
    """
    if modelname.endswith("pnon"):
        return True
    else:
        return False


def pdb2ms_main(file_1, file_2):
    """main function to convert pdb files to ms files with selected attractive residues.
    
    Args:
        file_1 (str): path to the first pdb file
        file_2 (str): path to the second pdb file
    """
    #can now use this because renumbering to IMGT numbering is done
    attractive_res = {'D':{'start' : [26, 55,104], 'end': [39,66,118]}, 'E': {'start' : [26, 55,104], 'end': [39,66,118]}, 'C':{'start' : [-1], 'end': [1000]}}
    
    # make a list of files to process
    file_list = [file_1, file_2]
    
    # loop over the files
    for file in file_list:
        p = Path(file)
        case_id = p.name[:3]
        model_type = p.name[4:6]
        prepared = isPrepared(p.stem)
        
        if p.suffix == ".pdb":
            pdb = open_file(p)  # Assuming open_file is a pre-defined function
            ms = pdb2ms(pdb)    # Assuming pdb2ms is a pre-defined function
            ms = add_attraction_tcr(ms, attractive_res)  # Assuming add_attraction_tcr is a pre-defined function
            output_file = str(Path(p.parent, p.stem)) + ".ms"
            print(f"Wrote to: {output_file}")
            write_ms(ms, output_file)  # Assuming write_ms is a pre-defined function
        else:
            print(f"Skipped {p.stem}")
