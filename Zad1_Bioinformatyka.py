import os
import re
from filecmp import dircmp
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser

import unittest
import sys
sys.path.append("..")
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, SingleLetterAlphabet
import forgi.graph.bulge_graph as fgb
import forgi.utilities.debug as fud
from forgi.utilities.exceptions import GraphConstructionError
import forgi.utilities.stuff as fus
import forgi.threedee.model.coarse_grain as ftmc

punktacja = {'dopasowanie':1, 'niedopasowanie':-1, 'przerwa':-1}
def pobranie_PDB():
    print("PDB")
def funkcja(string, tmplista):
    parser = PDBParser()
    structure = parser.get_structure('X', string)
    model=structure[0]
    for chain in model:
        tmplista.append([string[-8:-4],chain])
def sprawdzDopasowanie(x, y):
    if x == y:
        return punktacja['dopasowanie']
    elif x == "-" or y == "-":
        return punktacja['przerwa']
    else:
        return punktacja['niedopasowanie']

def NeedlemanWunsch(seq1, seq2, wyniki):
    m, n = len(seq1), len(seq2)
    punkty = np.zeros((m+1, n+1))  
 # Faza inicjalizacji macierzy---------------------------------------------------------
    for i in range(m+1):
        punkty[i][0] = punktacja['przerwa'] * i
    for j in range(n+1):
        punkty[0][j] = punktacja['przerwa'] * j
 # Wypelnienie macierzy punktacji-------------------------------------------------------
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            przek = punkty[i-1][j-1] + sprawdzDopasowanie(seq1[i-1], seq2[j-1])
            D = punkty[i-1][j] + punktacja['przerwa']
            I = punkty[i][j-1] + punktacja['przerwa']
            punkty[i][j] = max(przek, D, I)
    i = m
    j = n
    dopasowanie_1 = "" 
    dopasowanie_2 = ""
 # Przejscie przez macierz-----------------------------------------------------------
    while (i>0 and j>0):
        punkty_teraz = punkty[i][j]
        punkty_przek = punkty[i-1][j-1]
        punkty_L = punkty[i][j-1]
        punkty_gora = punkty[i-1][j]
       
        if punkty_teraz == punkty_przek + sprawdzDopasowanie(seq1[i-1], seq2[j-1]):
            dop_1 = seq1[i-1]
            dop_2 = seq2[j-1]
            i = i-1
            j = j-1
        elif punkty_teraz == punkty_gora + punktacja['przerwa']:
            dop_1 = seq1[i-1]
            dop_2 = "-"
            i -= 1
        elif punkty_teraz == punkty_L + punktacja['przerwa']:
            dop_1 = "-"
            dop_2 = seq2[j-1]
            j -= 1
        dopasowanie_1+=dop_1
        dopasowanie_2+=dop_2
            
    while i>0:
        dop_1 = seq1[i-1]
        dop_2 = "-"
        dopasowanie_1+=dop_1
        dopasowanie_2+=dop_2
        i -= 1
    while j>0:
        dop_1 = "-"
        dop_2 = seq2[j-1]
        dopasowanie_1+=dop_1
        dopasowanie_2+=dop_2
        j -= 1
    
    dopasowanie_1 = dopasowanie_1[::-1]
    dopasowanie_2 = dopasowanie_2[::-1]
    sekwencja = len(dopasowanie_1)
    punktacjaSekwencji = 0
    identycznosc = 0
    for i in range(sekwencja):
        dop_1 = dopasowanie_1[i]
        dop_2 = dopasowanie_2[i]
        if dop_1 == dop_2:
            identycznosc += 1
            punktacjaSekwencji += sprawdzDopasowanie(dop_1, dop_2)
        else: 
            punktacjaSekwencji += sprawdzDopasowanie(dop_1, dop_2)
        
    identycznosc = identycznosc/sekwencja * 100
    wy = ' '.join(['Po dopasowaniu', '\nsekwencja 1:',  str(dopasowanie_1), '\nsekwencja 2:', str(dopasowanie_2), '\nProcent identycznosci:', str(identycznosc), '\nPunktacja:', str(punktacjaSekwencji), '\n'])
    wyniki.write(wy)

    print("Po dopasowaniu")
    print("sekwencja 1:",dopasowanie_1)
    print("sekwencja 2:",dopasowanie_2)
    print("Procent identycznosci: %2.1f" % identycznosc)
    print("Punktacja:", punktacjaSekwencji)
    
    return identycznosc

    

def metrykaGory(S1, S2):
    vS1 = []
    vS2 = []
    S1c = 0
    S2c = 0

    for j in S1:
            if j == '(': 
                S1c+=1
                vS1.append(S1c)
            elif j == ')':
                S1c-=1
                vS1.append(S1c)
            else:
                vS1.append(S1c)

    for j in S2:
            if j == '(': 
                S2c+=1
                vS2.append(S2c)
            elif j == ')':
                S2c-=1
                vS2.append(S2c)
            else:
                vS2.append(S2c)

    return abs(sum(vS1)-sum(vS2))

def metryka(s1, s2): #procent niezgodnych pozycji
    if len(s1)!=len(s2):
        return 100
    d=0
    for i in range(len(s1)):
        if s1[i]!=s2[i]:
            d=d+1
    return 100*d/len(s1)
def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
def get_pdb_file(PDBlist2,path):
    pdbl = PDBList()
    for i in PDBlist2:
        #pdbl.retrieve_pdb_file(i)
        pdbl.retrieve_pdb_file(i, file_format="pdb",pdir=path)
        #pdbl.update_pdb()
def convert_pdb_to_fasta_tring(pdb_file):
     cg = ftmc.from_pdb(pdb_file)
     return cg.to_fasta_string()
    
def check_lancuchow(pdb_file):
    amino_code = {
	'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
	'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
	'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M',
	'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T',
	'TRP':'W', 'TYR':'Y', 'VAL':'V', 'HIS':'H',
	'ASX':'B', 'GLX':'Z', 'UNK':'K'
    }
    fa = {}
    with open(pdb_file) as fh:
        for buff in fh:
            if (buff[0:4] != 'ATOM'):
                continue
            chain_name = buff[21:22]
            res_number = int(buff[22:26])
            amino_acid = buff[17:20]
            if not (chain_name in fa):
                fa[chain_name] = []
            aa = 'X'
            if (amino_acid in amino_code):
                aa = amino_code[amino_acid]
            if (len(fa[chain_name]) != res_number):
                fa[chain_name] += ['X'] * (res_number - len(fa[chain_name]))
            fa[chain_name][res_number - 1] = aa
    for k, v in sorted(fa.items()):
	#print (len(k))
        if not(len(fa) >=2 and set(''.join(v))=={'X'}): 
            return False
        return True

def drzewoDFS(elem):
    print("wejscie: ", elem)

    #return wyjscie
def utworzZbiorZPliku(sciezka):
    f = open(sciezka, "r")
    print(f.name)
    struktury = []
   
    for lin in f:
        struktury.append(lin)
    return struktury
if __name__ == "__main__":
     path = Path('/pdb')
     lista = list(path.glob('**/*.ent'))

     sciezka='/pdb'
     parser = PDBParser()
     
     import Bio
     from Bio.PDB import PDBList
     listaLancuchow=[]
     PDBlist1=[]
     
     PDBlist=['1EVV','1EHZ','1ESY','1DDY','2G1W','1ALK','2lbk','3g78','4jkw']
     #pdbl = PDBList()
     get_pdb_file(PDBlist,path)
     for i in lista:
         if check_lancuchow(i) is True:
             PDBlist1.append(i)
             print(i)
     lista_struktury = 'pdb/Struktury.txt'
     listy_struktury = utworzZbiorZPliku(lista_struktury)
     bg1 = fgb.BulgeGraph()
     bg2 = fgb.BulgeGraph()
      
     _1evv_ ='(((((((..((((.....[..)))).((((.........)))).....(((((..]....))))))))))))....'
     _1ehz_ ='(((((((..((((.....[..)))).((((.........)))).....(((((..]....))))))))))))....'
     _1esy_ = '((((.((....))..))))'
     _1ddy_ ='......[[[.{((....((]]]...).).}.))..'
     _2g1w_ = '((((.[[..))))......]]'
     _2lbk_ ='.(((((.....))))).'
     _1clq_ ='(((((......(((....(((....)))....)))...)))))'
     _1f7g_ ='((((((.(((((....)))))))))))'
     _1f79_ ='.(((((.(((((....)))))))))).'
     listy = list([_1evv_,_1ehz_,_1esy_,_1ddy_,_2g1w_,_2lbk_,_1f7g_,_1f79_])
     bg = fgb.BulgeGraph()
     path_wynik =r"wynik/wyniki.txt"
     
     wyniki = open(path_wynik, "w")
     #wyj=''
     for ind, z in enumerate(listy):
         bg.from_dotbracket(z)
         elem_str = bg.to_element_string()
         print("struktury",ind+1,":", z)
         print("            ",elem_str)
         wyj = ' '.join(['struktury',str(ind+1),':', str(z) , '\n','            ',str(elem_str)])
     #wyniki.write(wyj)
     #wynik1 =''
     for ind, z in enumerate(listy):
         for ind2, z2 in enumerate(listy):
             if ind == ind2:
                 continue
             bg1.from_dotbracket(z)
             elem_str1 = bg1.to_element_string()
             bg2.from_dotbracket(z2)
             elem_str2 = bg2.to_element_string()
             print("=============================================================================================\n")
             print("odleglosc struktury", ind+1, "od struktury", ind2+1, "wynosi")
             NeedlemanWunsch(elem_str1, elem_str2, wyniki)
             print("metryka gory:", metrykaGory(z, z2))
             print("metryka procent niezgodnych pozycji:",metryka(z, z2))
             print("=============================================================================================\n")
     ##wyniki.write(wynik1)
