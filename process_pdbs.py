import sys
import numpy as np
import pandas as pd
from itertools import zip_longest
import csv

import json
import requests
from time import sleep

import Bio
from Bio import PDB
import Bio.PDB
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio import SeqIO
from Bio.PDB import PDBParser, MMCIFIO
from Bio.PDB import PDBList, PDBIO, FastMMCIFParser, Select
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.Polypeptide import PPBuilder

class ChainSelect(Select):
    def __init__(self, chain):
        self.chain = chain
    def accept_chain(self, chain):
        if chain.get_id() == self.chain:
            return 1
        else:          
            return 0


file = open('Database2.csv')
csvreader = csv.reader(file)

rows = []
for row in csvreader:
    rows.append(row)
#print(rows) # start from rows[1] b/c first row is simply header
file.close()

parser = FastMMCIFParser(QUIET = True)
#structure = parser.get_structure('4n9f', '4n9f.cif')

io = MMCIFIO()
io.set_structure(structure)
io.save("raw_pdbs_xray/"+rows[i][0].lower()+'.cif', ChainSelect('U'))

   
# how to find which sequence is host and which is pathogen
# mkdir for polypeptide id into the correct folders


# extra folder - unknown, host, pathogen, coordinates     

    
