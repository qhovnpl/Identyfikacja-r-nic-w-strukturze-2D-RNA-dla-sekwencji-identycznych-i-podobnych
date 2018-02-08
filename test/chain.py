from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import *

parser = PDBParser()
ppb = PPBuilder()
structure = parser.get_structure('tmp', '1ESY.pdb')
polypeptides = ppb.build_peptides(structure[0]['A'])
sequence = str(polypeptides[0].get_sequence())

print( sequence)
