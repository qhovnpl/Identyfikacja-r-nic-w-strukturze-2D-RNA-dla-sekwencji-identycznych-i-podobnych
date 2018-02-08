#!/usr/bin/python

# This program is a simple python code to clean PDB structures and automatically adds the structures to the FREAD database.
# I separately design this code in python since I understand that everyone has struggled to parse the PDB format.
# A bad PDB structure can cause "Segmentation Fault".
# There are plenty of ways to clean PDB structures. You can freely modify this code for your own sake.
# The file name must start with 4 letters of PDB code.
# This program will cut the query structure chain by chain and produce "[code][chain].pdb" files.
# The original file will be removed.
# The database file should be located in the same folder as the folder containing PDB structures.
# This programme only works for X-ray determined structures.

import sys, glob, getopt

REARRANGE = 0
# !!!!! Important: If you set this variable '1', it will rearrange all the first residue of chains starting with 1. !!!!!!!!!!!!
# You can turn this option off, simply by setting it '0'

global mainchain

mainchain = ["N", "CA", "C", "O"]

def usage():
    print ("./FREAD-Database.py [Expression] [optional: -c[chain identifier]")
    print( "\nFor example,")
    print ("$ ./FREAD-Database.py \"*.pdb\"")
    print ("$ ./FREAD-Database.py 2v9t.pdb")
    print ("     - This will parse all chains from the \"2v9t.pdb\" file")
    print( "$ ./FREAD-Database.py 2v9t.pdb -cA")
    print ("     - This will parse only A chain from the \"2v9t.pdb\" file")
    sys.exit()

def Parsing(file):
    chains = {}
    open_file = open(file).readlines()
    for line in open_file:
        if line[:4] == "ATOM":
            Chain = line[21]
            atom = line[13:15].strip()
            Residue_Number = int(line[22:26].strip())
            if not chains.has_key(Chain):
                chains[Chain] = {Residue_Number:{atom:line}}
            else:
                if not chains[Chain].has_key(Residue_Number):
                    chains[Chain][Residue_Number] = {atom:line}
                else:
                    if mainchain.count(atom):
                        if chains[Chain][Residue_Number].has_key(atom):
                            pass
                        else:
                            chains[Chain][Residue_Number][atom] = line
                    else:
                        if not chains[Chain][Residue_Number].has_key("Sidechain"):
                            chains[Chain][Residue_Number]["Sidechain"] = [line]
                        else:
                            chains[Chain][Residue_Number]["Sidechain"].append(line)
        elif line[:6] == "ENDMDL":
            break
    return chains

def MakePDB(file, Residue, coordinate, shift):
    write_file = open(file, "w")
    for res_num in Residue:
        Coordinate_key = set(coordinate[res_num].keys())
        New_Resnum = "%4d"%(res_num-shift)
        if coordinate[res_num].has_key("Sidechain"):
            all_atom = mainchain + ["Sidechain"]
        else:
            all_atom = mainchain
        if not Coordinate_key == set(all_atom):
            pass
        else:
            for atoms in all_atom:
                if mainchain.count(atoms):
                    line = coordinate[res_num][atoms]
                    New_line = line[:22] + New_Resnum + line[26:]
                    write_file.write(New_line)
                else:
                    try:
                        for sidechain in coordinate[res_num][atoms]:
                            New_line = sidechain[:22] + New_Resnum + sidechain[26:]
                            write_file.write(New_line)
                    except:
                        pass
    return 0

if __name__ == "__main__":
    # if len(sys.argv) < 2:
    #     usage()
    #else:
        FILES = glob.glob('/1ESY.pdb')
        optlist = getopt.getopt('A', "c:")
        
        if not optlist[0]:
            chain_wanted = []
        else:
            chain_wanted = optlist[0][0][1]

        for file in FILES:
            PDB_code = file[:4]
            Coordinate = Parsing(file)
            if chain_wanted:
                for i in Coordinate.keys():
                    if i != chain_wanted:
                        del Coordinate[i]
                    else: pass
            else: pass
            for chain in Coordinate:
                Residue_No = Coordinate[chain].keys()
                Residue_No.sort()
                First_ResNum = Residue_No[0]
                FileName = "%s%s.pdb"%(PDB_code, chain)
                if not REARRANGE:
                    MakePDB(FileName, Residue_No, Coordinate[chain], 0)
                else:
                    MakePDB(FileName, Residue_No, Coordinate[chain], First_ResNum - 1)
