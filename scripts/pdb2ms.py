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
from sys import argv
from pathlib import Path


def open_file(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()
    return lines

def pdb2ms(pdb):
    ms = []
    for line in pdb:
        ms_line = line[:54]
        attraction = "   0.0"
        ms_line = ms_line + attraction
        ms.append(ms_line)
    return ms

def add_attraction_tcr(ms, attractive_res):
    new_ms = []
    for line in ms:
        #print("ATOM   4163  CH2 TRP E 253      31.771  -0.145  97.650   0.0".startswith("ATOM"))
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
    f = open(filename, "w")
    for line in ms:
        f.write(line + "\n")
    f.close()

def isPrepared(modelname):
    if modelname.endswith("pnon"):
        return True
    else:
        return False

def pdb2ms_main(outdir, attractive_res):
    #benchmark
    # attractive_res = {'D':{'start' : [26, 55,104], 'end': [39,66,118]}, 'E': {'start' : [26, 55,104], 'end': [39,66,118]}, 'C':{'start' : [-1], 'end': [1000]}}

    #lyra
    #attractive_res = {'D':{'start' : [26, 57, 106], 'end': [40, 69, 128]}, 'E': {'start' : [25, 53, 103], 'end': [37, 64, 129]}, 'C':{'start' : [1], 'end': [8]}}

    #tcr_model2
    #attractive_res = {'D':{'start' : [26, 57, 106], 'end': [41, 69, 139]}, 'E': {'start' : [26, 57, 106], 'end': [41, 69, 139]}, 'C':{'start' : [1], 'end': [8]}}

    #capri
    #attractive_res = {'D':{'start' : [27, 50, 92], 'end': [35, 60, 103]}, 'E': {'start' : [27, 49, 92], 'end': [33, 56, 108]}, 'C':{'start' : [-1], 'end': [1000]}}
    #attractive_res = {'D':{'start' : [27, 59, 107], 'end': [40, 72, 139]}, 'E': {'start' : [27, 59, 107], 'end': [40, 68, 139]}, 'C':{'start' : [-1], 'end': [1000]}}
    p = Path(outdir)
    for model in p.iterdir():
        case_id = model.name[:3]
        model_type = model.name[4:6]
        prepared = True#isPrepared(model.stem)
        if model.suffix == ".pdb":
            pdb = open_file(model)
            ms = pdb2ms(pdb)
            ms = add_attraction_tcr(ms, attractive_res)
            print("Wrote to: ", str(Path(p, model.stem)) + ".ms")
            write_ms(ms , str(Path(p, model.stem)) + ".ms")
        else:
            print("Skipped {}".format(model.stem))




if __name__ == "__main__":
     pdb2ms_main("/home/nils/swifttcr/example/output/")