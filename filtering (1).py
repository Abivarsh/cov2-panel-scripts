#!/usr/bin/env python
# coding: utf-8

# # Basic Workflow
# 
# ### DATA FILTERING
# 
# This Jupyter notebook takes the output of a the cellranger pipeline and returns structures that are ready for synthesis
# 
# This assumes that the user is ordering their vectors in a monosystronic format and it provides a couple of output files which should aid in the diagnosis of problems.
# 
# Versions:
#     Python 3.7.6
# 
# Requirements:
# 
#     Bio Python
#     
#     igblast installed with a version of germline databases set up
#     
#     
#     Either:
#         1.
#             Data that has been processed with Cellranger
#         2.
#             cellranger installed and available at the path 
# 
#             cellranger germline references installed 
#     
#     
#     

# In[1]:


from Bio import SeqIO
from shutil import copyfile
import csv
import os
import functools
import re
import tarfile
import zipfile
import subprocess


# In[2]:


#Constants - Do Not Change
chainTypes = ['IGL','TRA', 'IGK', 'TRD', 'Multi', 'None',  'TRB', 'TRG', 'IGH'] # Chain Types called by cellranger's pipeline
allIsotypes = ["IGHA1", "IGHA2", "IGHD", "IGHE", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM", "IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC7", "TRAC", "TRBC1", "TRBC2", "TRDC","TRGC1","TRGC2"]


# In[3]:


# Intended to be changed by the user to impact the execution of the pipeline.


folders = ["4658-AT-1","4658-AT-2","4658-AT-3","4658-AT-4"] # sample names from the cellranger pipeline

outputs = ['outs/consensus.fastq', 'outs/clonotypes.csv', 'outs/filtered_contig.fastq',"outs/all_contig.fastq", "outs/web_summary.html", "outs/vloupe.vloupe", "outs/consensus_annotations.csv", "outs/all_contig_annotations.csv", "outs/filtered_contig_annotations.csv","outs/metrics_summary.csv", "outs/cell_barcodes.json"]# Files that you want to provide as part of the tar. 
exportsDir = "pipelineOutputs" #The name of the directory where the output files will be stored
numberOfTopSequences = 999999999 # Number of sequences that you plan to order -- set high to disable filtering
outputPrefix = "COV2-" #Prefix that proceeds the clone names
outputIdStart = 1 # Where output numbering should be started
keepChains = ['IGL', 'IGH', 'IGK'] # The values of chain type that are going to be kept 
isotypesToKeep = ['IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', "IGHA1", "IGHA2"] # Heavy Chain Isotypes that must be in a cell for it to be valid

vhLeader = "ACCGCCGCCACCATGGAGTTCGGTCTTAGCTGGGTGTTTCTTGTCGCCCTGTTCAGAGGGGTACAATGC"
vconstantLinkerVlLeader='GCCTCCACAAAGGGTCCAAGTGTCTTCCCACTGGCTCCCAGCAGCAAGAGTACTTCAGGTGGGACTGCAGCTCTCGGGTGCCTGGTCAAGGACTACTTTCCCGAGCCCGTAACAGTATCTTGGAACTCCGGTGCTCTGACAAGTGGAGTGCATACTTTCCCAGCTGTGTTGCAGTCAAGCGGGTTGTACTCCCTCAGTAGTGTAGTTACTGTCCCTTCATCTTCACTGGGGACTCAAACCTACATTTGTAACGTGAATCACAAACCAAGCAATACTAAAGTAGATAAGAAGGTGGAGCCAAAAAGTTGTGATAAAACTCATACTTGTCCCCCCTGTCCTGCACCAGAGCTGTTGGGCGGTCCCAGTGTATTTTTGTTTCCCCCTAAACCCAAAGACACACTGATGATTTCTCGAACTCCCGAAGTGACCTGCGTCGTCGTTGATGTAAGTCACGAAGACCCCGAGGTAAAATTCAATTGGTACGTAGACGGCGTAGAGGTGCATAACGCTAAAACCAAACCAAGGGAGGAACAATACAACAGCACTTACAGAGTCGTATCTGTACTGACTGTTCTCCACCAAGACTGGCTTAATGGGAAGGAGTATAAATGCAAGGTGTCAAACAAGGCTCTGCCTGCCCCCATAGAGAAAACCATAAGTAAAGCTAAAGGACAACCTAGAGAGCCCCAGGTTTATACTCTTCCCCCCTCCCGAGACGAGCTGACCAAGAACCAAGTTTCCTTGACCTGTCTGGTTAAGGGTTTTTATCCAAGCGATATAGCCGTAGAATGGGAGAGCAACGGACAACCTGAGAATAACTACAAAACAACTCCCCCAGTGCTGGACTCAGATGGCTCATTTTTCCTGTATTCAAAGCTCACCGTGGACAAATCTCGGTGGCAGCAAGGGAATGTATTCTCCTGTTCCGTCATGCACGAGGCGCTGCATAATCACTATACCCAGAAATCTCTGTCCCTTTCTCCTGGAAAGCGGAGAAAACGGGGATCAGGGGAAGGTAGGGGGTCATTGCTGACTTGCGGTGACGTCGAGGAGAATCCTGGACCTATGGACATGCGAGTCCCTGCCCAACTGCTGGGTCTTCTCCTCCTTTGGCTGTCAGGAGCCCGCTGC'

humanVGeneReference = "/home/taylor/software/ncbi-igblast-1.16.0/Ig/human/human_gl_V"
humanDGeneReference = "/home/taylor/software/ncbi-igblast-1.16.0/Ig/human/human_gl_D"
humanJGeneReference = "/home/taylor/software/ncbi-igblast-1.16.0/Ig/human/human_gl_J"
auxData             = "/home/taylor/software/ncbi-igblast-1.16.0/aux_data/aux/human_gl.aux"
igblastIGDATA     = "/home/taylor/software/ncbi-igblast-1.16.0"

# A dictionary of controls that should be included in the order
# {"clone_name":[List of chains - must have one heavy and one light - order does not matter]}
controlsDict ={"CR3022": ["cagatgcaactggtgcagagcgggaccgaggtgaaaaaacccggcgagtccttgaagattagttgcaaaggttcaggttatggattcattacctactggataggatgggtaagacaaatgccaggtaaaggtttggagtggatggggataatataccctggggattctgagacaagatactcccccagcttccagggtcaggtaactatatctgcagataaaagtataaacactgcctatcttcagtggtcctcactgaaagcatctgacactgctatctactactgtgccggaggtagtggcatatcaactccaatggacgtatggggtcaagggacaactgtcaccgtctccagcg", "gacatacagttgacacaatctccagattctctggcagtatcacttggggaacgtgcaacaataaattgcaaatcatcccagtctgtactctacagcagcatcaataaaaattacctcgcatggtatcagcaaaagccaggacaacccccaaaacttctcatatattgggcaagtacacgggagtctggcgtccccgataggttttctggtagtgggagtggtactgactttacccttaccatatcatcacttcaagccgaggatgtcgctgtgtactactgccaacaatattattcaaccccttatacattcggccagggcacaaaggtagaaattaaa"]}

"""
filterOneToOne: Cells must be one to one
removeAllWithExtraChains removes cells that have transcprits that you would not expect to find in B Cells - only if you have more data than you can order
minUmiFilter: Requires all kept consensuses to have at least X UMIS to be considered valid
isotypeFiltering: Require there to be at least one chain in the kept chains array for a cell to stay in the collection
"""
filters = {
    "filterOneToOne": True,
    "removeAllWithExtraChains": False,
    "minUmiFilter": True,
    "isotypeFiltering": True
}
minUmisToBeValid = 2
intermediateLogging = True
numThreads = 12


# In[4]:


# Computed from your intputs
dropChains = (set(chainTypes) - set(keepChains)) #Chains that will not be included
isotypesToDiscard = list(set(allIsotypes) - set(isotypesToKeep)) # Isotypes that do not want to be kept
os.environ["IGDATA"] = igblastIGDATA


# In[5]:


# Function definitions

# Lambda functions used for filtering
getUmis = lambda x : int(x['umis'])
getReads = lambda x : x['reads']
getConsensusName = lambda x : x['raw_consensus_id']
getProductive = lambda x : x['productive'] == "True"
getFullLength = lambda x : x['full_length'] == "True"
getIsCell = lambda x : x['is_cell'] == "True"
getHighConfidence = lambda x : x['high_confidence'] == "True"
flattenLists = lambda lists: [val for sublist in lists for val in sublist]
getCellAllGood = lambda x: set(map(getChainAllGood, flattenLists(list(x[1].values())))) == {True}
getChainAllGood = lambda x: x['high_confidence'] == "True" and x['productive'] == "True" and x['full_length'] == "True" and x['is_cell'] == "True"
getFilterOneToOne = lambda x: len(x[1]['IGH']) == 1 and ((len(x[1]['IGK']) + len(x[1]['IGL'])) == 1)
getCellContainsNoDropableChain = lambda x: set(map(getDropChain, flattenLists(list(x[1].values())))) == {False}
getDropChain = lambda x: x['chain'] in dropChains
getAtLeastXUmis = lambda x: int(x['umis']) >= minUmisToBeValid
getCellLimitedUmiChain = lambda x: set(map(getAtLeastXUmis, flattenLists(list(x[1].values())))) == {True}

# Take a dictionary and apply a given filter to the key and the value elements. Return a dictionary of the elements that pass the filter.
def filterToDict(fil, oldDict):
    replacement = {}
    for key, value in oldDict:
        if fil((key,value)):
            replacement[key] = value
    return replacement

#Translate a nucleotide sequence into an amino acid sequence.
def translate(seq): 
    
    #codon table mapping codons to Amino Acids
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    }
    
    protein ="" 
    
    if len(seq) % 3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= table[codon] 
            
    return protein 

# This did not end up getting used
def lookupClonotype(element, keyedOnClonotypeId):
    
    heavy = pairing['IGH']

    if pairing.get('IGL',False):
        light = pairing['IGL']
    if pairing.get('IGK',False):
        light = pairing['IGK']

    return keyedOnClonotypeId[heavy['sequence_id']]

# Return the percent identity for of both chains 
def getPercentIdentityBothChains(heavy, light):
    
    bases = 0
    correct = 0
    
    for thing in [heavy, light]:
        stuff = getBasesCorrect(thing)
        correct += stuff[0]
        bases += stuff[1]
        
    return correct / bases

# Seaprate a cells chains out into heavy and light and return them
def getHeavyLight(cell):

    heavy = cell['IGH']
    try:
        light = cell['IGK'] if cell['IGK'] else ['IGL']
    except:
        try:
            light = cell['IGK']
        except:
            light = cell['IGL']
            
    return heavy, light

# Take in a cell and return the percent identity of the chains.
def getPercentIdentityBothChainsCell(cell):
    heavy, light = getHeavyLight(cell)
    return getPercentIdentityBothChains(heavy, light)

# Get the percent identity of a single chain
def getPercentIdentity(element):
    vBases = (int(element['v_germline_end']) - int(element['v_germline_start']))
    vCorrect = vBases * float(element['v_identity'])
    jBases = (int(element['j_germline_end']) - int(element['j_germline_start']))
    jCorrect = jBases * float(element['j_identity'])
    
    if element['d_call']!='':
        dBases = (int(element['d_germline_end']) - int(element['d_germline_start']))
        dCorrect = dBases * float(element['d_identity'])
    else:
        dBases = 0
        dCorrect =0
    
    return ((vCorrect+dCorrect+jCorrect) / (vBases+dBases+jBases))

# Hand this function a single chain in the format listed above and it returns the number of bases and the number of bases mapped from germline
def getBasesCorrect(element):
    
    vBases = (int(element['v_germline_end']) - int(element['v_germline_start']))
    vCorrect = vBases * float(element['v_identity'])
    jBases = (int(element['j_germline_end']) - int(element['j_germline_start']))
    jCorrect = jBases * float(element['j_identity'])
    
    if element['d_call']!='':
        dBases = (int(element['d_germline_end']) - int(element['d_germline_start']))
        dCorrect = dBases * float(element['d_identity'])
        
    else:
        dBases = 0
        dCorrect =0
        
    return ((vCorrect+dCorrect+jCorrect), (vBases+dBases+jBases))

# From the igblast returned TSV file this will return the needed sequence
def getRawSequence(element):
    return element['sequence']

# This will grab the variable region from the igblast result formatted in TSV format
def getVariableRegion(element):
    """
    There is a situation where the alignment does not extend all the way to the end of the genes.
        If the first or last codon does not match germline then the ideal alignment will involve simply clipping the first or last codon
    This function uses the cigar strings for the gene matches in order to determine if the alignment needs to be extended
    """
    
    v_cigar_split = re.findall("\d+[NSM=XDI]", element['v_cigar'])
    vShift = 0
    
    if v_cigar_split[1][-1] =="N":
        vShift = int(v_cigar_split[1][:-1])    

    j_cigar_split = re.findall("\d+[NSM=XDI]", element['j_cigar'])    
    jExtension = 0
    
    if j_cigar_split[-1][-1] == "N":
        jExtension = int(j_cigar_split[-1][:-1])
    
    sequence = element['sequence'][int(element['v_sequence_start'])-vShift-1:int(element['j_sequence_end'])+jExtension]
    
    endSequence = len(sequence)%3
    
    if endSequence != 0:
        sequence = sequence[:-1*endSequence]
    
    return sequence
    
# extract the cdr3
def getCdr3(element):
    return element['cdr3']

# get cdr3AA
def getCdr3AA(element):
    return element['cdr3_aa']
    
# Count the number of element which are present in 
def getCleanPairingClonotypeCount(element):
    
    clonotypeName = element['IGH']['sequence_id']
    heavy = keyedOnClonotypeId[clonotypeName]['IGH']
            
    if keyedOnClonotypeId[clonotypeName].get('IGL',False):
        light = keyedOnClonotypeId[clonotypeName]['IGL']
    if keyedOnClonotypeId[clonotypeName].get('IGK',False):
        light = keyedOnClonotypeId[clonotypeName]['IGK']

    return len(v3jPairDict[(heavy[0]['v_gene'], heavy[0]['cdr3'], heavy[0]['j_gene'], light[0]['v_gene'], light[0]['cdr3'], light[0]['j_gene'])])

# Format a row in the CSV from the TSV igblast output
def createRow(pairing):
    results = {}
    
    heavy = pairing['IGH']
    if pairing.get('IGL',False):
        light = pairing['IGL']
    if pairing.get('IGK',False):
        light = pairing['IGK']
        
    results["10x_unique_identifier"] = heavy['sequence_id']
    try:
        results["clean_cells_count_of_clonotype"] = getCleanPairingClonotypeCount(pairing)
    except:
        results["clean_cells_count_of_clonotype"] = ""

    try:
        results["heavy_isotype"] = keyedOnClonotypeId[results["10x_unique_identifier"]]['IGH'][0]['c_gene']
    except:
        results['heavy_isotype'] = ""
    results["light_isotype"] = light['locus']
    
    results['heavy_untrimmed'] = getRawSequence(heavy)
    results['light_untrimmed'] =  getRawSequence(light)
    
    heavyVariableRegion = getVariableRegion(heavy)
    if len(heavyVariableRegion) % 3 ==0:
        results["heavy_nt_alignment"] = heavyVariableRegion
    else:
        results["heavy_nt_alignment"] = heavyVariableRegion[0:-1*(len(heavyVariableRegion) % 3)]
    
    results["heavy_nt_alignment_length"] = len(results["heavy_nt_alignment"])
    
    lightVariableRegion = getVariableRegion(light)
    if len(lightVariableRegion) % 3 ==0:
        results["light_nt_alignment"] = lightVariableRegion
    else:
        results["light_nt_alignment"] = lightVariableRegion[0:-1*(len(lightVariableRegion) % 3)]
    
    results["light_nt_alignment_length"] = len(results["light_nt_alignment"])
    results["heavy_aa"] = translate(results["heavy_nt_alignment"])
    results["heavy_aa_length"] = len(translate(results["heavy_nt_alignment"]))
    results["light_aa"] = translate(results["light_nt_alignment"])
    results["light_aa_length"] = len(translate(results["light_nt_alignment"]))
    results['heavy_v'] = heavy['v_call']
    results['heavy_d'] = heavy['d_call']
    results['heavy_j'] = heavy['j_call']
    results['heavy_cdr3_aa'] = getCdr3AA(heavy)
    results['heavy_cdr3_aa_length'] = len(getCdr3AA(heavy))
    results['heavy_percent_identity'] = getPercentIdentity(heavy)
    results['light_v'] = light['v_call']
    results['light_j'] = light['j_call']
    results['light_cdr3'] = getCdr3(light)
    results['light_cdr3_aa'] = getCdr3AA(light)
    results['light_cdr3_aa_length'] = len(getCdr3AA(light))
    results['light_percent_identity'] = getPercentIdentity(light)
    results['heavy_minus_light_percent_identity']  = getPercentIdentity(heavy) - getPercentIdentity(light)
    results['percent_identity_both'] = getPercentIdentityBothChains(heavy, light)
    try:
        results['heavy_umi_count'] = umisPerConsensus[heavy['sequence_id']]['IGH']
    except:
        results['heavy_umi_count'] = ""
    try:
        results['light_umi_count'] = umisPerConsensus[heavy['sequence_id']].get('IGK',False) if  umisPerConsensus[heavy['sequence_id']].get('IGK',False) else umisPerConsensus[heavy['sequence_id']].get('IGL',False)
    except:
        results['light_umi_count'] = ""
    results['VH leader'] = vhLeader
    results['V constant linker light leader']= vconstantLinkerVlLeader
    results['leader vh linker leader vl'] = results['VH leader'] + results["heavy_nt_alignment"] + results['V constant linker light leader'] + results["light_nt_alignment"]    
    
    return results


# In[6]:


# Create the directory for the output files of the pipeline. If it already exists don't stop the execution of the pipeline
try:
    os.mkdir(exportsDir)
except:
    pass


# In[7]:


# For each of the folders above grab the relevant files and include them in the outputs directory
for folder in folders:
    os.makedirs(os.path.join(exportsDir,"inputFiles",folder, "outs"), exist_ok=True)
    for output in outputs:
        copyfile(os.path.join(folder, output),os.path.join(exportsDir,'inputFiles', folder,output))


# In[8]:


# Go through all of the outs directories and read in their data 
# Change the names to ensure that their definitions will all be unique
allChains = []

for folder in folders:
    file = os.path.join(folder, "outs", "filtered_contig_annotations.csv")
    with open(file, 'r') as f:
        dictreader = csv.DictReader(f)
        for element in dictreader:
            element["barcode"] = folder + "_" + element["barcode"]
            element["contig_id"] = folder + "_" + element["contig_id"]
            element["raw_clonotype_id"] = folder + "_" + element["raw_clonotype_id"]
            element["raw_consensus_id"] = folder + "_" + element["raw_consensus_id"]
            allChains.append(element)

if intermediateLogging:
    print("Chains Observed:", len(allChains))


# In[9]:


# Iterate all of the transcripts and group them into cells
pairingsChainGroup = {} 

for element in allChains:
    try:
        pairingsChainGroup[element['barcode']][element['chain']].append(element)
    except:
        pairingsChainGroup[element['barcode']] = {'IGL':[],'TRA': [], 'IGK': [], 'TRD': [], 'Multi': [], 'None':[],  'TRB': [], 'TRG': [], 'IGH': []}
        pairingsChainGroup[element['barcode']][element['chain']].append(element)
        
if intermediateLogging:
    print("Cell Barcodes Observed: ",len(pairingsChainGroup.keys()))


# In[10]:


# Cells must have 1 heavy and one light Essential to later assumptions.
# This stage is currently required
pairingsChainGroup = filterToDict(getFilterOneToOne, pairingsChainGroup.items())

if intermediateLogging:
    print("One To One Cells", len(pairingsChainGroup.keys()))
        
# Only keep chains that are Productive, High Confidence, Full Length and Is Cell
oneToOneDict = pairingsChainGroup
pairingsChainGroup = filterToDict(getCellAllGood, pairingsChainGroup.items())
prodAndHigh = pairingsChainGroup

if intermediateLogging:
    print("Productive and high confidence and Full Length and Is Cell", len(pairingsChainGroup.keys()))


# Only keeps cells that contain transcripts if contamination is high this may need to be turned off.
if filters["removeAllWithExtraChains"]:
    pairingsChainGroup = filterToDict(getCellContainsNoDropableChain, pairingsChainGroup.items())

    if intermediateLogging:
        print("Drop Any cell which still contains a garbage chain", len(pairingsChainGroup.keys()))
    
# priorToUMIRemoval = pairingsChainGroup
if filters["minUmiFilter"]:
    pairingsChainGroup = filterToDict(getCellLimitedUmiChain, pairingsChainGroup.items())

    if intermediateLogging:
        print("Both transcripts have at least", minUmisToBeValid, "UMIs", len(pairingsChainGroup.keys()))
    
# Remove things that are not in your keept chains
if filters["isotypeFiltering"]:
    pairingsChainGroup = filterToDict(lambda x : set(map(lambda y : y['c_gene'] in isotypesToKeep, x[1]['IGH'])) & set([True]) == set([True]) , pairingsChainGroup.items())

    if intermediateLogging:
        print("Limit the search to cells containing a given heavy chain:", len(pairingsChainGroup))
    


# In[11]:


# Count the number of unique CDR3s that are still in cirulation at this point
cdr3s = []
for element in pairingsChainGroup.values():
    cdr3s.append(element['IGH'][0]['cdr3_nt'])

if intermediateLogging:
    print("total cdr3s:", len(cdr3s))
    print("Total Unique Cdr3s:", len(list(set(cdr3s))))


# In[12]:


# This change is going to break everything in the shorterm but it is the quickest way to hammer out the v3j change 
clonotypePairDict = {}
v3jPairDict = {}

for element in pairingsChainGroup.values():
    
    try:
        adding = (element['IGH'][0]['cdr3_nt'], element['IGK'][0]['cdr3_nt'])
        adding2 = (element['IGH'][0]['v_gene'], element['IGH'][0]['cdr3'],element['IGH'][0]['j_gene'], element['IGK'][0]['v_gene'], element['IGK'][0]['cdr3'],element['IGK'][0]['j_gene'])
        
    except:
        adding = (element['IGH'][0]['cdr3_nt'], element['IGL'][0]['cdr3_nt'])
        adding2 = (element['IGH'][0]['v_gene'], element['IGH'][0]['cdr3'],element['IGH'][0]['j_gene'], element['IGL'][0]['v_gene'], element['IGL'][0]['cdr3'],element['IGL'][0]['j_gene'])
        
    try:
        clonotypePairDict[adding].append(adding2) # offer the ability to perform lookups of every clonotype with a cdr3
        v3jPairDict[adding2].append((element['IGH'][0]['raw_clonotype_id']))

    except:
        clonotypePairDict[adding] = [adding2]
        v3jPairDict[adding2] = [element['IGH'][0]['raw_clonotype_id']]
        
if intermediateLogging:
    print("CDR3 Pairs (nucleotide): ", len(clonotypePairDict.keys()))
    print("V3J Pairs: ", len(v3jPairDict.keys()))
    print("Total Cells in play at this stage of the pipeline: ", functools.reduce(lambda x,y: x+y, list(map(lambda x: len(x[1]), v3jPairDict.items()))))


# In[13]:


# Order the list by Cell count. If you are ordering all of them then this may not be needed
orderedByCellCount = sorted(v3jPairDict.items(), key=lambda x: len(x[1]), reverse=True)

# Take the top numberOfTopSequences from the list.
selectedCdr3s = list(map(lambda x: (x[0],len(x[1]), len(set(x[1]))), orderedByCellCount))[:numberOfTopSequences]


# In[14]:


pulledClonotypes = list(functools.reduce(lambda x, y: set(x) | set(y), list(map(lambda x: x[1], v3jPairDict.items()))))


# In[15]:


folderClonotypes = {}

for element in pulledClonotypes:
    s = "_"
    folder = s.join(element.split("_")[:-1])
    clonotype = element.split("_")[-1]
    try:
        folderClonotypes[folder].append(clonotype)
    except:
        folderClonotypes[folder] = [clonotype]


# In[16]:


# Concatenate all of this stuff together
chains = {} 
foundElements =[]

for folder in folderClonotypes.keys():
    
    for record in SeqIO.parse(os.path.join(folder,"outs", "consensus.fastq"),"fastq"):
        
        # In the short term we are just going to run all of the clonotypes
        if record.name.split("_")[0] in folderClonotypes[folder]:
            
            try:
                record.id = folder + "_" + record.name.split("_")[0]
                chains[folder + "_" + record.name.split("_")[0]].append(record)

            except:
                record.id = folder + "_" + record.name.split("_")[0]
                chains[folder + "_" + record.name.split("_")[0]] = [record]

with open(os.path.join(exportsDir,"orderedClonotypes.fasta"), 'w') as handle:
    for element in chains.values():
        SeqIO.write(element,handle, "fasta")


# In[17]:


#  auxData,
with open(os.path.join(exportsDir, "allKeptChains.tsv"),"w") as f:
    p = subprocess.Popen(["igblastn", "-germline_db_V", humanVGeneReference, "-germline_db_D", humanDGeneReference, "-germline_db_J", humanJGeneReference,  "-num_threads", str(numThreads), "-outfmt",  "19", "-query", os.path.join(exportsDir, "orderedClonotypes.fasta")],stdout=f)
    stderr = p.communicate()
    p.wait()
    print(stderr[1])


# In[18]:


# Read in the results of the igblast processing

with open(os.path.join(exportsDir, "allKeptChains.tsv"), 'r') as f:
    dictreader = csv.DictReader(f, dialect='excel-tab')
    allItems = list(dictreader)


# In[19]:


sequencePairs = {}

for element in allItems:

    try:
        sequencePairs[element['sequence_id']][element['locus']] = element
        
    except:
        sequencePairs[element['sequence_id']] = {}
        sequencePairs[element['sequence_id']][element['locus']] = element


# In[20]:


# Produce a datastructure to lookup elements based on clonotype id 

keyedOnClonotypeId = {}

for key,  value in pairingsChainGroup.items():
   
    if not keyedOnClonotypeId.get(value['IGH'][0]['raw_clonotype_id'], False):
        keyedOnClonotypeId[value['IGH'][0]['raw_clonotype_id']] = {'IGL': [], 'IGH': []}
    
    if value['IGH'][0]['chain'] == "IGH":
        keyedOnClonotypeId[value['IGH'][0]['raw_clonotype_id']] = value


# In[21]:


# 
umisPerConsensus = {}

for folder in folders:
    fname = os.path.join(folder, "outs", "consensus_annotations.csv")
    with open(fname, 'r') as f:
        dictReader = csv.DictReader(f)
        for element in dictReader:
            
            element['clonotype_id'] = folder + "_" + element['clonotype_id']
            try:
                umisPerConsensus[element['clonotype_id']][element['chain']] = element['umis']
            except:
                umisPerConsensus[element['clonotype_id']] = {}
                umisPerConsensus[element['clonotype_id']][element['chain']] = element['umis']


# In[22]:


batchedOutput = [] 
clonotypeId =0
v3jIdenticalCount = {}


# v3jPairDict
for element in v3jPairDict.items():
    
    for element2 in element[1]:
        
        cell = sequencePairs[element2]
        heavy = cell['IGH']
        light = cell['IGK'] if cell.get('IGK', False) else cell.get("IGL", False)
        
        if not v3jIdenticalCount.get(element[0],False):
            v3jIdenticalCount[element[0]] = {}
            
        if not v3jIdenticalCount[element[0]].get((heavy['sequence_alignment_aa'], light['sequence_alignment_aa']), False):
            v3jIdenticalCount[element[0]][(heavy['sequence_alignment_aa'], light['sequence_alignment_aa'])] = 0
            
        v3jIdenticalCount[element[0]][(heavy['sequence_alignment_aa'], light['sequence_alignment_aa'])] += 1
        
# How should the clonotypes be sorted?
for key, value in v3jPairDict.items():
    v3jPairDict[key] = sorted(value, key=lambda x : getPercentIdentityBothChainsCell(sequencePairs[x]))
        
        
printedSequences = []

for element in v3jPairDict.items():

    for element2 in element[1]:

        cell = sequencePairs[element2]
        
        heavy = cell['IGH']
        light = cell['IGK'] if cell.get('IGK', False) else cell.get("IGL", False)
        
        if (heavy['sequence_alignment_aa'], light['sequence_alignment_aa']) not in printedSequences:
            
            printedSequences.append((heavy['sequence_alignment_aa'], light['sequence_alignment_aa']))

            result = createRow(sequencePairs[element2])
            result['number_of_exact_aa_matches'] = v3jIdenticalCount[element[0]][(heavy['sequence_alignment_aa'], light['sequence_alignment_aa'])]
            result['clonotype_group_characteristics'] = element[0]
            result['v3j id'] = clonotypeId
            batchedOutput.append(result)
        
    clonotypeId += 1
    
headers = list(batchedOutput[0].keys())

headers.remove('clonotype_group_characteristics')
headers.remove('v3j id')
headers.remove('number_of_exact_aa_matches')
headers = ['clonotype_group_characteristics', 'v3j id', "number_of_exact_aa_matches"] + headers

with open("pipelineOutputs/precursorOutput.csv","w") as f:
    dictwriter = csv.DictWriter(f, headers)
    dictwriter.writeheader()
    dictwriter.writerows(batchedOutput)


# In[23]:


# prep controls
with open(os.path.join(exportsDir, "controls.fasta"),'w') as controlFile:
    for key, value in controlsDict.items():
        for element in value:
            controlFile.write("> " + key + "\n")
            controlFile.write(element + "\n")    


# In[24]:


# subprocess.call(["igblastn", "-germline_db_V", humanVGeneReference, "-germline_db_D", humanDGeneReference, "-germline_db_J", auxData, "-num_threads", numThreads, "-outfmt 19", "-query ", os.path.join(exportsDir, "controls.fasta")," > " + os.path.join(exportsDir, "controls.tsv")])
#  auxData,
with open(os.path.join(exportsDir, "controls.tsv"),"w") as f:
    p = subprocess.Popen(["igblastn", "-germline_db_V", humanVGeneReference, "-germline_db_D", humanDGeneReference, "-germline_db_J", humanJGeneReference,  "-num_threads", str(numThreads), "-outfmt",  "19", "-query", os.path.join(exportsDir, "controls.fasta")],stdout=f)
    stderr = p.communicate()
    p.wait()
    print(stderr[1])


# In[25]:


with open(os.path.join(exportsDir, "controls.tsv"), 'r') as f:
    dictreader = csv.DictReader(f, dialect='excel-tab')
    controlList = list(dictreader)
    
controls = {}
for element in controlList:
    if not controls.get(element['sequence_id'], False):
        controls[element['sequence_id']] = {}
        
    if element['locus'] == 'IGH':
        controls[element['sequence_id']]['IGH'] = element
        
    if element['locus'] == 'IGK':
        controls[element['sequence_id']]['IGK'] = element
        
    if element['locus'] == 'IGL':
        controls[element['sequence_id']]['IGL'] = element
        
controlLines = []
for element in controls.values():
    
    temp = createRow(element)
    temp["Sequence name (Enter the names of your sequences in this column)  - optional"] = temp['10x_unique_identifier']
    controlLines.append(temp)


# In[26]:


# Generate and output the finaldownselection

toWrite = []
orderedAminoAcidSequences = [] 

# Iterate over the row of remaining items and format a dictionary that can be written out 
# If an amino acid sequence has already been ordered then skip it.
for key, value in v3jPairDict.items():
    element = createRow(sequencePairs[value[0]])
    if (element['heavy_aa'], element['light_aa']) not in orderedAminoAcidSequences:
        toWrite.append(element)
        orderedAminoAcidSequences.append((element['heavy_aa'], element['light_aa']))
    else:
        print("Sequence skipped because it was already ordered")

counter = outputIdStart

for element in toWrite:
    uId = str(counter)
    counter += 1
    element['Sequence name (Enter the names of your sequences in this column)  - optional'] = outputPrefix + uId
    
headers = ["Sequence name (Enter the names of your sequences in this column)  - optional"] + list(toWrite[0].keys())[0:-1]

with open("pipelineOutputs/finalOutputs.csv","w") as f:
    dictwriter = csv.DictWriter(f, headers)
    dictwriter.writeheader()
    dictwriter.writerows(controlLines)
    dictwriter.writerows(toWrite)


# In[27]:


with zipfile.ZipFile(exportsDir+".zip","w") as myzip:
    myzip.write(exportsDir)
    myzip.write(os.path.join(exportsDir,"testing"))
# We had a command to upload the outputs directory to a shared Box folder. Something else can be added in its place.


# In[28]:


with zipfile.ZipFile(exportsDir+".zip","w") as myzip:
    
    for roots, dirs, files in os.walk(exportsDir):

        myzip.write(roots)

        for filename in files:
            myzip.write(os.path.join(roots,filename))
            

