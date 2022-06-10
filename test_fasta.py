import Bio
import re
from Bio import PDB
import Bio.PDB
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio import SeqIO
from Bio.PDB import PDBParser, MMCIFIO


f = open('rcsb_pdb_1A14.fasta','r')
lines = f.readlines()
print(lines)

for i in range(len(lines)):
    print(str(i) + " " + lines[i])

for i in range(0, len(lines), 2):
    directory = ""
    if ("VIRUS" in lines[i]) or ("virus" in lines[i]):
        directory = "pathogen_chains/"
    else:
        directory = "host_chains/"

#hre=re.compile('>(\S+)')
#outh = hre.search(lines[0])
#if outh:
    #print(outh.group(1))

#from Bio import SeqIO

#for record in SeqIO.parse("rcsb_pdb_1A14.fasta", "fasta"):
    #print(record.id)



sequences = [i for i in SeqIO.parse("rcsb_pdb_1A14.fasta", "fasta")]
sequence_a = sequences[0]
print(sequence_a.name)
