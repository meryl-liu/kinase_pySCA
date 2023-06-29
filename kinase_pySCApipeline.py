#import the necessary modules

%matplotlib inline
from __future__ import division
import os
import time
import matplotlib.pyplot as plt
import matplotlib
import scipy.cluster.hierarchy as sch
import numpy as np
import copy
import math
from pysca import scaTools as sca
from scipy.stats import t
from scipy.stats import scoreatpercentile 
from scipy.stats import fisher_exact
import pickle as pickle
from Bio import SeqIO
from Bio import Entrez
from Bio import PDB
from bioservices import UniProt
import requests
import subprocess
import time
import pandas as pd

#### Function Definitions

class KinaseEntry:
    def __init__(self, GI_number, family, species, MSA_index, UniProt_ID):
        self.GI_number = GI_number
        self.family = family
        self.species = species
        self.MSA_index = MSA_index
        self.UniProt_ID = UniProt_ID
    
 
def parseAlignment(alignment, selection, type):
    ''' create a function that, using a original, annotated but unprocessed multiple sequence alignment as its input, 
    returns all the desired proteins of a certain family or species, the sequence index in the alignment, and the 
    protein name/PDB ID (if it exists?) alignment, selection, and type are strings. The index must be multiplied by 2, +1
    to reflect the actual index in the alignment fasta file. '''  
    columnNum = 0;
    if type == 'family':
        columnNum = 1;
        printNum = 2
    else:
        columnNum = 2
        printNum = 1
        
    # read in the .fasta or .an file containing the alignment
    fasta_file = '../sca/'+ alignment
    records = SeqIO.parse(fasta_file, 'fasta')

    # define the output path and text name
    output_file = '../output/'+alignment+'_'+selection+'_index.txt'
    
    # initialize the empty dictionary of kinase MSA entries, each dictionary entry is a KinaseEntry object
    kinase_dict = {};
    
    # Open the output file in write mode
    with open(output_file, 'w') as f:
    # Loop through the records
        for MSA_index, record in enumerate(records, start=0):
            #print(index)
            header = record.description.split('|');
            header = [part.strip(' ') for part in header] # strip out all whitepsace
            #print(header)
            if header[columnNum] == selection:
                
                # Extract the protein GI number
                gi_number = header[0]
                #uniprot_ID = get_uniprotACC(gi_number)
                #print(uniprot_ID)
                # Extract the protein sequence family and species
                family = header[1]
                species = header[2]
                
                # Write the protein GI numbers to a file, for copy-pasting into the UniProt ID mapping tool
                f.write(f"{gi_number}\n")
                kinase_dict[gi_number] = KinaseEntry(gi_number, family, species, MSA_index, None)

    return kinase_dict


def appendUniProtIDs(dict, tableFile):
    ''' function appendUniProtIDs takes the kinase dictionary constructed from parseAlignment and the .tsv file
    from the UniProt ID Mapping, and adds the UniProt ID to the attribute of the Kinase Entry in each dictionary item'''
    counter = 0;
    df = pd.read_table(tableFile)
    for index, row in df.iterrows():
        counter+=1;
        dict[str(row['From'])].UniProt_ID = row['Entry']
        #print(dict[str(row['From'])].UniProt_ID, counter)
    #print(df)
    #return df
''' 
parseAlignment function creates a python dictionary with GI numbers as keys to protein atrributions (family, species, index), but also returns a .txt file in ../output that contains a list of all the GI numbers corresponding to the species/family of interest. This list is to be copy and pasted into the UniProt ID mapping tool, to download a tsv file of all the correspomding successfully mapped UniProt Primary Accession Numbers.'''

# main program
kinase_dict = parseAlignment('masterAln.an', 'Homo sapiens', 'species')
appendUniProtIDs(kinase_dict, 'HumanKinase_UniProtID.tsv') # .tsv comes from takign the .txt output of parseAlignment and mapping

# After running appendUniProtIDs, kinase_dict is now a dictionary that contains keys in the form of GI numbers, leading to KinaseEntry 
# objects that have specific attributes about each protein



