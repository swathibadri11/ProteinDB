# -*- coding: utf-8 -*-

import sys
import numpy as np
import pandas as pd
from itertools import zip_longest
import itertools
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


URL="https://search.rcsb.org/rcsbsearch/v1/query"

class queryPDB:
    def __init__(self,query):
        self.rows = 100
        self.return_type = "entry"
        self.query = query
     
    def make_params(self, start=0):
        return dict(
            query=self.query,
            return_type=self.return_type,
            request_options=dict(pager=dict(start=start, rows=self.rows)),
        )
        
    def single_query(self, start=0):
        parms = self.make_params(start)
        response=requests.get(url = URL,params={'json': json.dumps(parms, separators=(',', ':'))})
        response.raise_for_status()

        if response.status_code == requests.codes.OK:
            return response.json()
    
        elif response.status_code == requests.codes.NO_CONTENT:
            return None
        else:
            raise Exception(f"Unexpected status: {response.status_code}") 
    
    def get_identifiers(self, response):
        identifiers = [result["identifier"] for result in response["result_set"]]
        return identifiers
    
    def __iter__(self):
        
        start = 0
        response = self.single_query(start=start)
        if response is None:
            return
        identifiers = self.get_identifiers(response)
        start += self.rows
        total = response["total_count"]
        
        if len(identifiers) > 0:
            yield from identifiers
        
        while start < total:
            response = self.single_query(start)
            identifiers = self.get_identifiers(response)
            start += self.rows
            yield from identifiers

#query only for proteins 
query_xray = {
    "type": "group",
    "logical_operator": "and",
    "nodes": [
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "exact_match",
          "value": "Viruses",
          "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name"
        }
      },
        {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "exact_match",
          "value": "Eukaryota",
          "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name"
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "exact_match",
          "value": "X-RAY DIFFRACTION",
          "attribute": "exptl.method"
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "less_or_equal",
          "value": 3.5,
          "attribute": "rcsb_entry_info.resolution_combined"
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "greater",
          "attribute": "rcsb_entry_info.nonpolymer_entity_count",
          "value": 1
        }
      }
    ]
  }

query_em = {
    "type": "group",
    "logical_operator": "and",
    "nodes": [
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "exact_match",
          "value": "Viruses",
          "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name"
        }
      },
        {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "exact_match",
          "value": "Eukaryota",
          "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name"
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "exact_match",
          "value": "ELECTRON MICROSCOPY",
          "attribute": "exptl.method"
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_entry_info.resolution_combined",
          "operator": "range",
          "value": {
            "from": 0,
            "include_lower": True,
            "to": 6,
            "include_upper": True
          }
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "greater",
          "attribute": "rcsb_entry_info.nonpolymer_entity_count",
          "value": 1
        }
      }
    ]
  }

query_nmr = {
    "type": "group",
    "logical_operator": "and",
    "nodes": [
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "exact_match",
          "value": "Viruses",
          "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name"
        }
      },
        {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "exact_match",
          "value": "Eukaryota",
          "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name"
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "exact_match",
          "value": "pdbx_nmr_spectrometer.model",
          "attribute": "exptl.method"
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "less_or_equal",
          "value": 3.5,
          "attribute": "rcsb_entry_info.resolution_combined"
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "greater",
          "attribute": "rcsb_entry_info.nonpolymer_entity_count",
          "value": 1
        }
      }
    ]
  }

ids_xray = queryPDB(query_xray)

ids_em = queryPDB(query_em)

ids_nmr = queryPDB(query_nmr)

all_ids_xray = [i for i in ids_xray]

#all_ids = all_ids_xray + all_ids_em + all_ids_nmr
all_ids = [i for i in all_ids_xray]


no_pdb = len(all_ids)
#no_pdb = 50

pdbl=PDB.PDBList()
for id in all_ids_xray[0:no_pdb]:
    # download the specified pdbs
    pdbl.retrieve_pdb_file(id,pdir='raw_pdbs_xray/')
    # download the fasta file also - easier to find species 
    
# module 1 organism
organism_xray = [None]*no_pdb
virus_xray = [None]*no_pdb

def pdb_organism(i, dict_xray):
    if "_entity_src_nat.pdbx_organism_scientific" in dict_xray:
        for k in dict_xray["_entity_src_nat.pdbx_organism_scientific"]:
            if ("virus" in k) or ("VIRUS" in k):
                virus_xray[i] = k
            else:
                organism_xray[i] = k
    return organism_xray

# module 1.5 virus - add _pdbx_entity_src_syn.organism_scientific
def pdb_virus(i, organism_xray, dict_xray):
    if "_entity_src_gen.pdbx_gene_src_scientific_name" in dict_xray:
        for k in dict_xray["_entity_src_gen.pdbx_gene_src_scientific_name"]:
            if (not(k in organism_xray) and ("virus" in k) and ("VIRUS" in k)):
                virus_xray[i] = k
            else:
                organism_xray[i]= k
                
    #else if "_pdbx_entity_src_syn.organism_scientific" in dict_xray:
        #for k in dict_xray["_pdbx_entity_src_syn.organism_scientific"]:
            #if (not(k in organism_xray) and ("virus" in k) and ("VIRUS" in k)):
                #virus_xray[i] = k
            #else:
                #organism_xray[i] = k
            

# module 2
pubid_xray = []

def pdb_pubid(dict_xray):
    for k in dict_xray["_citation.pdbx_database_id_PubMed"]:
        if (k != "?"):
            pubid_xray.append(k)
    return pubid_xray

# module 3 - helper module
max = 0

def polymer_no(dict_xray):
    polymer = 0
    for k in range(len(dict_xray["_entity.type"])):
        if (dict_xray["_entity.type"][k] == "polymer"):
            polymer += 1
    return polymer

def k(dict_xray):
    k_list = []
    for k in range(len(dict_xray["_entity.type"])):
        if (dict_xray["_entity.type"][k] == "polymer"):
            k_list.append(k)
    return k_list


# module 4
def polymer_name(polymer_names, k_list, dict_xray):
    for k in range(len(k_list)):
        polymer_names[k].append(dict_xray["_entity.pdbx_description"][k_list[k]])
    return polymer_names

# module 5 - helper module
def ids(dict_xray):
    return len(dict_xray["_pdbx_struct_assembly_gen.assembly_id"])

# module 6

# maybe do a dictionary inside a dictionary 
def assembly(dict_xray):
    return dict_xray["_pdbx_struct_assembly_gen.assembly_id"]

# module 7 - does not need to be shown in the table, write and read this dictionary, important to test that you are correctly writing and reading the informtaion, structure of nested dict should not be lost
def asym(dict_xray):
    return dict_xray["_pdbx_struct_assembly_gen.asym_id_list"]    
              
# main method to call all sub-modules

polymer_no_list = []

for i in range(len(all_ids_xray[0:no_pdb])):
    dict_xray = MMCIF2Dict("raw_pdbs_xray/"+str(all_ids_xray[i]).lower()+".cif")
    
    org_ids = pdb_organism(i, dict_xray)
    pdb_virus(i, org_ids, dict_xray)

    pdb_pubid(dict_xray)

    #print(organism_xray)
    #print(virus_xray)

    polymer_no_list.append(polymer_no(dict_xray))

    if (polymer_no(dict_xray) > max):
        max = polymer_no(dict_xray)

# module 8 - geometry

#def geometry(structure):
    # geometry description
        # cartesian - x, y, z in space position (angstroms)
        # internal coordinates - translation invariant
            # translate molecules but coordinates will not change
            # describing them w respect to molecule itself
            # 2 distances, 1 angle --> define geometry of molecule independent of reference
    #p = PDBParser()
    #structure = p.get_structure(structure)
    #for model in structure:
        #for chain in model:
            #for residue in chain:
                #return residue.get_coord()

    #with coordinates, we want them in independent files, so we can take different chains and compare them
            # host pathogen interactions of a certain kind of virus, compare 2 different virus binding receptoprs to different hosts
            # libraries to save coordinates in their own file
                # chain select and extract 
    

polymer_names = [[] for i in range(max)]


max_id = 0
for i in range(len(all_ids_xray[0:no_pdb])):
    dict_xray = MMCIF2Dict("raw_pdbs_xray/"+str(all_ids_xray[i]).lower()+".cif")
    polymer_name(polymer_names, k(dict_xray), dict_xray)

    if (ids(dict_xray) > max_id):
        max_id = ids(dict_xray)

    #print(polymer_names)

    #polymer_names[i] = polymer_names[i][0:polymer_no(dict_xray)]

   
assembly_ids = []
asym_ids = []

for i in range(len(all_ids_xray[0:no_pdb])):
    dict_xray = MMCIF2Dict("raw_pdbs_xray/"+str(all_ids_xray[i]).lower()+".cif")
    # save chains as mmcif files themselves 
    assembly_ids.append(assembly(dict_xray))
    asym_ids.append(asym(dict_xray))

# how to initialize table
# make empty table, initalize table at the bginning, add info as you are getting it
# work with dictionaries instead

# nested dictionaries would give us more flexiblity to pickle into file
# submodule to take relevant info into a table
# currently: id, pubid, organism, virus, polymer number, polymer names, map id

# calling it: read everything, organize into nested dict, write file, open empty python script w file, access all info, use pandas data frame to write a samll table that would contain some of the info we would need about the files
print()

#transpose + crop + re-transpose + if element = none, delete
polymer_names_copy = polymer_names
polymer_names = list(itertools.zip_longest(*polymer_names))

for i in range(len(all_ids_xray[0:no_pdb])):
    dict_xray = MMCIF2Dict("raw_pdbs_xray/"+str(all_ids_xray[i]).lower()+".cif")
    polymer_names[i] = polymer_names[i][0:polymer_no(dict_xray)]

#for i in range(len(polymer_names)): # max number of columns 
    #for j in range(len(polymer_names[i]): # rows on each one
        #if (polymer_names[i][j] == "None"):
            #polymer_names.remove(polymer_names[i][j])

polymer_names = list(itertools.zip_longest(*polymer_names))
 
d = [all_ids_xray, pubid_xray, organism_xray, virus_xray, polymer_no_list]
for i in range(len(polymer_names)):
    print(polymer_names[i])
    d.append(polymer_names[i])

export_data = zip_longest(*d, fillvalue = '')
with open('Database.csv', 'w', encoding="ISO-8859-1", newline='') as myfile:
      wr = csv.writer(myfile)
      wr.writerow(("X-Ray IDs", "Pub ID", "X-Ray Organisms", "X Ray Virus", "No. of Chains", "Name of Chains"))
      wr.writerows(export_data)

myfile.close()

# unique chains that form the complex - check
# want to extract chains from the righ sub complex 
# sub directory called structures/ coordinates, with host pathogen directories
    # host and pathogen mmcif files saved
# in pathogen file, one chain for b, o, etc
# reconstruct history files
# can have table as pickle or csv file

# extract which chain where
# get sequences as well
# ppBuilder
# geometry of complexes, sequences will be imp to annotate
# no need to add to the table 

  


# do the table, save cif file for each entry
# new python script - process structures
    # input: table
    # iterate over table and create 2 folders: host proteins and pathogen proteins
          # each row - extract the chain associated with host to a host folder, same for pathogen
          # name each file by the pdb code and chain identifier

# directory
    # scripts, complex mmcif file, host mmcif files, pathogen mmcif file, host coord, pathogen coord 
    # get the coordinates from them individually
          # extract file tells u which ones are chains
          # chain select
          # which one of the chains is host and which one is pathogen 
                

# Notes - case sensitive, remove all commas, format the same, ribosome thing in the first lines of the table
# Estimate looking through table, how many have live chain, heavy chain, or antibody 
        # first pass of counting antibodies by counting heavy and light chain
        # some labelelled as fab or heavy or light or antibody (case sensitive)
        # also count how many have DNA/RNA


