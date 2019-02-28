# (1) ALPHA-HELIX PARSER

# Number of helices in the protein.
helix_no=0
# Total length of all helices used to calculate the average length of a helix.
total_hel_length=0

# Here a PDB file of choice is opened and parsed accordingly.
PDB=open("trim21.pdb", "r")
file=PDB.readlines()
print("helixID" + "\t\t" + "length")
for line in file:
    if "HELIX" in line[0:6]: 
        helix_len=line[71:76]
        helix_id=line[11:14]
        total_hel_length+=int(helix_len)
        helix_no+=1
        print(helix_id + "\t\t" + helix_len)        
PDB.close()

print("\nThere are", helix_no, "alpha-helices in this protein")
print("\nMean length of alpha-helix is", "{:.2f}".format(total_hel_length/float(helix_no)))

#######################

# (2) BETA-SHEET PARSER

# Array of non-unique strand sequence numbers for all sheets where x=initSeqNumber, y=endSeqNumber, z=chain ID
strands=[]

# Dictionary with sheet IDs and corresponding initSeqNumber, endSeqNumber and chain ID for further calculations
dic={}

PDB=open("trim21.pdb", "r")
file=PDB.readlines()
for line in file:
    if "SHEET" in line[0:6]:
        initSeqNum=line[22:26]
        endSeqNum=line[33:37]
        chainID=line[21:22]
        sheetID=line[11:14]
        strands.append((initSeqNum, endSeqNum, chainID))
        if sheetID not in dic.keys():
            dic[sheetID]=[(initSeqNum, endSeqNum, chainID)]
        else:
            if (initSeqNum, endSeqNum, chainID) not in dic[sheetID]:
                dic[sheetID].append((initSeqNum, endSeqNum, chainID))
PDB.close()

print("The number of sheets in this protein is", len(dic))

"""
Note: some sheets contain identical strands but from different chains. 
For such sheets, the no. of strands remains as listed. 
For the sheets listing multiple identical strands from the same chain, 
only one copy of the strand is retained by the program.
"""

# Number of strands per sheet ID.
print("\nThe number of strands for each sheet is:\n")
for x in dic:
    print(x, ":", len(dic[x]))

print("\nTotal number of residues for each sheet is:\n")

# Total number of residues per each sheet ID. Adds 1 to get the count, not the difference.
for i in dic:
    res=sum([int(y) - int(x) + 1 for (x, y, z) in dic[i]])
    print(i, ":", res)

# Set of unique strands across all sheets. Contains unique initSeqNumber, endSeqNumber and chain ID.
unique_str=set(strands)    

# Average no. of strands per sheet: no. of unique strands/no. of sheets
mean_str_no=len(unique_str)/float(len(dic))

# Total no. of residues for unique strands across all sheets. Adds 1 to get the count, not the difference.
total_sheet_res=0
if unique_str==[]:
    total_sheet_res=total_sheet_res
else:
    total_sheet_res=sum([int(y) - int(x) + 1 for (x, y, z) in unique_str])

# Average number of residues in a sheet: total no. of residues from unique strands/no. of sheets 
mean_res_no=total_sheet_res/float(len(dic))

print("\nThe average number of unique strands in a sheet is", "{:.2f}".format(mean_str_no))
print("\nThe average number of residues in a sheet is", "{:.2f}".format(mean_res_no))

#######################

# (3) ATOM ENTRY PARSER

""""Please note that some residues are often missing in the ATOM entries.
    As specified by REMARK 465 in this file, 7 residues are missing in the ATOM entries for chain A.
    Therefore, the total residue count may not be fully reflected by these entries."""

# A list to store residue coordinates.
residues=[]

# A dictionary to store unique polypeptide chains and corresponding length.
chain_info={}

PDB=open("trim21.pdb", "r")
file=PDB.readlines()
for line in file: 
    if "ATOM" in line[0:6]:
        chainID=line[21:22]
        resName=line[17:20]
        resSeqNum=line[22:26]
        if chainID not in chain_info.keys():
            chain_info[chainID]=0     
        residues.append((resName, chainID , resSeqNum))  
PDB.close()

residue_set=set(residues)
#List of unique residue coordinates across all polypeptide chains.
residue_set=list(residue_set)

print("The number of polypeptide chains is", len(chain_info))

# Total length of each chain is computed here and stored as a dictionary value for its key.
print("\nThe residue length of each polypeptide chain is:\n")
for x in chain_info:
    length=0
    for y in residue_set:
        if x in y:
            length+=1
            chain_info[x]=length
    print(x, ":", chain_info[x])

# Percentage of residues in a helix across all polypeptide chains. total_hel_length obtained from Q1
helix_perc=100*(float(total_hel_length)/len(residue_set))

# Percentage of residues in a sheet across all polypeptide chains. total_sheet_res obtained in Q2
sheet_perc=100*(float(total_sheet_res)/len(residue_set))

# Percentage of residues in neither/irregular secondary structure
helix_sheet=total_hel_length + total_sheet_res
irreg_perc=100*((len(residue_set) - float(helix_sheet))/len(residue_set))

print("\nThe overall percentage of residues in different secondary structures is as follows:")
print("\nHelix: " + str("{:.2f}".format(helix_perc)) + "%" )
print("Sheet: " + str("{:.2f}".format(sheet_perc)) + "%" )
print("Neither/irregular: " + str("{:.2f}".format(irreg_perc)) + "%" )

#####################################

# (4) RESIDUE DISTRIBUTION CALCULATOR (accepts a single PDB file).

# This dictionary will store unique amino acids as keys and their frequencies as values.
amino_freq={}

# A complete list of amino acids only, obtained from ATOM entries for all popypeptide chains in the protein.
total_aa=[x for (x,y,z) in residue_set]
total_aa.sort()

# Add the individual amino acids present in the sequence as dictionary keys.
for aa in total_aa:
    if aa not in amino_freq:
        amino_freq[aa]=0
    else:
        continue 

# Write the header to a csv file. This will be followed by aa frequencies.
csv=open('trim21_aa_freq.csv', "a")
header="# Aminoacid" + "," + "p(aa)" + "\n"
csv.write(header)

#Compute the frequency of each aa. Add as a value for each key and write each row to the csv file.
for aa in amino_freq:
    amino_freq[aa]=total_aa.count(aa)/float(len(total_aa))
    amino_freq[aa]="{:.3f}".format(amino_freq[aa])
    row=aa + "," + amino_freq[aa] + "\n"
    csv.write(row)
csv.close()

###############################

# (5) GLOBAL RESIDUE CALCULATOR (accepts multiple PDB files; here 25 were chosen).

global_freq={}
residue_info=[]

fileList=["1b9e", "1bj4", "1cmi", "1f9q", "2piv", "1grt", "1hjm", "1ijk", "1ioc", "5og0", "2qpy", "5ons", 
     "1q2u", "4w93", "1umt", "1w24", "1za4", "2ccs", "2hgs", "3caa", "3cjb", "3jzk", "3tj2", "4fbx", "2v13"]

for fileName in fileList:
    PDB=open(fileName + ".pdb", "r")
    file=PDB.readlines()
    for line in file: 
        if "ATOM" in line[0:6]:
            chainID=line[21:22]
            resName=line[17:20]
            resSeqNum=int(line[22:26])
            # fileName appended so that, if by chance, there are identical residue coordinates between different files, they will be retained after running set()
            residue_info.append((fileName, resName, chainID , resSeqNum))  
PDB.close()

globalRes_set=set(residue_info)
#List of unique residue coordinates across all polypeptide chains for all listed files
globalRes=list(globalRes_set)

# Extract only the residues from the residue coordinates list and store in a new alphabetically sorted list
global_aa=[y for (x,y,z,w) in globalRes]
global_aa.sort()

# Add unique residues in the aggregated sequence as dictionary keys. Compute the global frequency of each aa.
for aa in global_aa:
        global_freq[aa]=global_aa.count(aa)/float(len(global_aa))
        
# Write the header to a csv file.
csv=open('global_aa_freq.csv', "a")
header="# Aminoacid" + "," + "p(aa)" + "\n"
csv.write(header)

# Write each amino acid and its corresponding frequency across the 25 files.
for aa in global_freq:
        row=aa + "," + "{:.3f}".format(global_freq[aa]) + "\n"
        csv.write(row)
csv.close()

###############################################

# (6) SECONDARY STRUCTURE PROPENSITY CALCULATOR
# * Computes residue frequency for each sequence segment separately.
# * Measures the propensity of each residue to be in a specific segment using pre-defined algorithms.

from math import log

# Calculates residue frequency for a specific segment (helix, beta, irregular).
def freq_count(seg):
    freq=seg.count(aa)/float(len(seg))
    return freq

# Calculates propensity of a residue to be in a specific segment.
def aa_prop(seg_freq): 
    prop=log(float(seg_freq[aa])/float(global_freq[aa]))
    return prop

# Alternative laplace smoothing method for calculating residue frequency in a specific segment.
def laplace_freq(seg):
    freq=(1+seg.count(aa))/float(len(seg))
    return freq

# Alternative laplace smoothing algorithm for calculating the propensity of each residue to be in a segment.
def laplace_prop(seg_freq, all_freq):
    prop=log(seg_freq/all_freq)

# Residue coordinates from helix segments for each file will be extracted here.
helix_seg=[]

# Residue coordinates from sheet segments for each file will be extracted here.
sheet_seg=[]

for fileName in fileList:
    PDB=open(fileName + ".pdb", "r")
    file=PDB.readlines()
    for line in file: 
        if "HELIX" in line[0:6]:
            initSeqNum=int(line[21:25])
            endSeqNum=int(line[33:37])
            chainID=line[19:20]
            seqRange=range(initSeqNum, endSeqNum+1)
            for res in globalRes:
                for n in seqRange:
                    if chainID in res and n in res and fileName in res:
                        if res not in helix_seg:
                            helix_seg.append(res)
                        else:
                            continue   
        if "SHEET" in line[0:6]:
            initSeqNum=int(line[22:26])
            endSeqNum=int(line[33:37])
            chainID=line[21:22]
            seqRange=range(initSeqNum, endSeqNum+1)
            for res in globalRes:
                for n in seqRange:
                    if chainID in res and n in res and fileName in res:
                        if res not in sheet_seg:
                            sheet_seg.append(res)
                        else:
                            continue
PDB.close()

# This list stores extracted helix residues from the helix residue coordinates list across the 25 files.
helixRes=[y for (x,y,z,w) in helix_seg]

# This list stores extracted sheet residues from the sheet residue coordinates list across the 25 files.
sheetRes=[y for (x,y,z,w) in sheet_seg]

# Coordinates from irregular structures are found by removing helix+sheet residues from the list of all residues.
irreg_seg=globalRes[:]
comb_seg=helix_seg[:]
comb_seg.extend(sheet_seg)
for res in comb_seg:
    irreg_seg.remove(res)

# List to store extracted residues found in irregular secondary structures.
irregRes=[y for (x,y,z,w) in irreg_seg]

# Dictionaries to store residue frequency for each secondary structure.
helix_freq={}
sheet_freq={}
irreg_freq={}
    
# Unique global residues used for references below.
ref_aa=list(set(global_aa))
ref_aa.sort()

# Residue frequency for each segment is calculated using the freq_count function defined earlier.
for aa in ref_aa:
    helix_freq[aa]=freq_count(helixRes)  
    sheet_freq[aa]=freq_count(sheetRes)
    irreg_freq[aa]=freq_count(irregRes)

# Dictionaries to store the relevant propensity of each aa.
propH={}
propS={}
propI={}

# Propensities for each residue are calculated using the aa_prop or laplace_prop (laplace smoothing) function.
# Laplace smoothing is only used if original algorithm results in errors caused by log operations on values of 0.
for aa in ref_aa:
    # re-calculated global aa frequencies where 1 is added to aa counts for laplace smoothing.
    global_laplace=(1+global_aa.count(aa))/float(len(global_aa))
    try:
        propH[aa]=aa_prop(helix_freq)
    except (ZeroDivisionError, ValueError):
        helix_laplace=laplace_freq(helixRes)
        propH[aa]=laplace_prop(helix_laplace, global_laplace)
    try:
        propS[aa]=aa_prop(sheet_freq)
    except (ZeroDivisionError, ValueError):
        sheet_laplace=laplace_freq(sheetRes)
        propS[aa]=laplace_prop(sheet_laplace, global_laplace)
    try:
        propI[aa]=aa_prop(irreg_freq)
    except (ZeroDivisionError, ValueError):
        irreg_laplace=laplace_freq(irregRes)
        propI[aa]=laplace_prop(irreg_laplace, global_laplace)
    
# Write all results to a new csv file.
csv=open('aa_stats.csv', "a")
header="# Aminoacid" + "," + "p(aa)" + "," + "p(aa|H)" + "," + "p(aa|S)" + "," + "p(aa|I)" + "," + "prop(H)" + "," + "prop(S)" + "," + "prop(I)" + "\n"
csv.write(header)

for aa in ref_aa:
    row=aa + "," + "{:.3f}".format(global_freq[aa]) + "," + "{:.3f}".format(helix_freq[aa]) + "," + "{:.3f}".format(sheet_freq[aa]) + "," + "{:.3f}".format(irreg_freq[aa]) + "," + "{:.3f}".format(propH[aa]) + "," + "{:.3f}".format(propS[aa]) + "," + "{:.3f}".format(propI[aa]) +  "\n"
    csv.write(row)

csv.close()

#############################

# (7) BARCHART GENERATOR TOOL: for residue distribution (single file) and propensity measures.

import matplotlib.pyplot as plt
import numpy as np

# A function to extract a list of values from each dictionary to be used by numpy.
def aa_list(dic_name):
    aa_list=[]
    for aa in dic_name:
        aa_list.append(aa)
    return aa_list

# A function to extract a list of aa keys from each dictionary to be used by numpy.
def value_list(dic_name):
    value_list=[]
    for aa in dic_name:
        value_list.append(float(dic_name[aa]))
    return value_list

# A function used to produce all barcharts plots.
def plot_bar(dic_name, ylabel, title):
    x_label=aa_list(dic_name)
    y_label=value_list(dic_name)
    index = np.arange(len(x_label))
    plt.bar(index, y_label)
    plt.xlabel("Amino Acid")
    plt.ylabel(ylabel)
    plt.xticks(index, x_label, rotation=90)
    plt.title(title)
    return plt.show()

# Compute a barchart plot of the residue distribution in a protein, here TRIM21 analysed earlier.
label="Frequency in Sequence"
bar_title="Amino Acid Frequencies in Protein"
plot_bar(amino_freq, label, bar_title)

# Compute barchart plots of the three residue propensities.
Hlabel="prop(H) = log(p(aa|H)/p(aa))"
Htitle="Propensity Measure of Amino Acids for Alpha-Helix"
plot_bar(propH, Hlabel, Htitle)

Slabel="prop(S) = log(p(aa|S)/p(aa))"
Stitle="Propensity Measure of Amino Acids for Beta-Sheet"
plot_bar(propS, Slabel, Stitle)

Ilabel="prop(I) = log(p(aa|I)/p(aa))"
Ititle="Propensity Measure of Amino Acids for Irregular Segments"
plot_bar(propI, Ilabel, Ititle)
