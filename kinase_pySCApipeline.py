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
import requests as r
import subprocess
import time
import pandas as pd
import MSA_pos_mapping_functions
import csv
#### Function Definitions

# Construct an alignment-to-structure mapping (ATS)
# This is used to convert between alignment numbering schemes and structure positions in a
# given PDB
def buildATS(alignment, indexFile, refposFileArg = ''):
    f = open(indexFile, 'r')
    protLines = f.readlines()
    ATSmap = {}
    seqIx = {}
    for k in protLines:
        prot,index =k.split()
        if refposFileArg == '':
            refposFile = prot+'.pos'
        else: 
            refposFile = refposFileArg
        refpos = 'Refpos/'+refposFile
        outputfile= 'Outputs/ATS_'+prot
        cmd1='../sca/scaMakeATS.py Inputs/'+alignment+ ' -i '+index+' -o '+refpos+' --output '+outputfile+' > '+outputfile+'.log'
        print(cmd1)
        os.system(cmd1)
        dbtmp= pickle.load(open(outputfile+'.db', 'rb'))
        ATSmap[prot] = (dbtmp['sequence']['ats'])
        seqIx[prot] = (int(index))
    if len(seqIx) == 1:
        return ATSmap[ATSmap.keys()[0]],seqIx[ATSmap.keys()[0]]
    else:
        return ATSmap,seqIx

#Identify all residues contacting another group of residues
#Used to define sector contacting positions.
def findContacts(distMat,distATS,resGroup):
    conn = []
    for i,pos in enumerate(resGroup):
        if pos in distATS:
            contacts = np.where(distMat[distATS.index(pos),:] < 4 )[0]
            contactPos = [distATS[k] for k in contacts]
            conn = conn + contactPos
    return set(conn)

# Read in a distance matrix (assumes csv format) - used for finding structural contacts
def readDistMat(matrixFile):
    f=open(matrixFile, 'r')
    dist = f.readlines()
    distMat = np.zeros([len(dist),len(dist)])
    for i,line in enumerate(dist):
        distMat[i,:] = line.split(',')
    return distMat

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
    output_file = '../sca/output/'+alignment+'_'+selection+'_index.txt'
    
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
                #f.write(f"{gi_number}\n")
                f.write(f"{gi_number} {family} {species} {MSA_index} \n")
                kinase_dict[gi_number] = KinaseEntry(gi_number, family, species, MSA_index, None)

    return kinase_dict


def appendUniProtIDs(dict, tableFile):
    ''' function appendUniProtIDs takes the kinase dictionary constructed from parseAlignment and the .tsv file
    from the UniProt ID Mapping, and adds the UniProt ID to the attribute of the Kinase Entry in each dictionary item'''
    #counter = 0;
    df = pd.read_table(tableFile)
    for index, row in df.iterrows():
        #counter+=1;
        dict[str(row['From'])].UniProt_ID = row['Entry']
        #print("GI #: " + dict[str(row['From'])].GI_number, "Uniprot ID: " + dict[str(row['From'])].UniProt_ID, counter)
    #print(df)
    #return df

def write_rows_to_csv(file_path, rows):
    with open(file_path, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        for row in rows:
            writer.writerow(row)

def cleanUpUniProtMappings(tsv_file):
    df = pd.read_table(tsv_file)
    df = df.drop_duplicates(subset="Entry", keep=False)  # Remove rows with duplicates (including the original row)
    df = df.drop_duplicates(subset="From", keep='first')  # Remove duplicate rows based on the specified column
    df.to_csv("HumanKinase_UniProtID_cleaned.tsv", sep='\t', index=False, quoting=csv.QUOTE_NONE)  # Save the modified DataFrame back to the CSV file


def generateSurfaceDECSV(dict):
    csv = [];
    #counter = 0;
    for key, value in dict.items():
        if value.UniProt_ID != None:
            csv_row = [value.UniProt_ID, "meryl-liu/kinase_pySCA/main/pdb/" + value.UniProt_ID + ".pdb", "A"]
            csv.append(csv_row)
        
    write_rows_to_csv("surfaceDE_input.csv", csv)
    

''' 
parseAlignment function creates a python dictionary with GI numbers as keys to protein atrributions (family, species, index), but also returns a .txt file in ../output that contains a list of all the GI numbers corresponding to the species/family of interest. This list is to be copy and pasted into the UniProt ID mapping tool, to download a tsv file of all the correspomding successfully mapped UniProt Primary Accession Numbers.'''

# main program
kinase_dict = parseAlignment('masterAln.an', 'Homo sapiens', 'species')
cleanUpUniProtMappings('HumanKinase_UniProtID.tsv') # cleans up .tsv file for duplicates, and returns a cleaned up tsv file
appendUniProtIDs(kinase_dict, 'HumanKinase_UniProtID_cleaned.tsv') # .tsv comes from takign the .txt output of parseAlignment and mapping

# After running appendUniProtIDs, kinase_dict is now a dictionary that contains keys in the form of GI numbers, leading to KinaseEntry 
# objects that have specific attributes about each protein

generateSurfaceDECSV(kinase_dict)

