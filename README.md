# PDBParser

## About
PDBParser is a computational tool used for parsing and analysing protein structural information from [Protein Data Bank](https://www.rcsb.org "RCSB PDB Homepage") (PDB) files. The main purpose of the tool is to analyse the frequency and composition of the two main elements of protein secondary structure, namely alpha-helices and beta-sheets. Currently, PDBParser accepts PDB files containing single and multiple polypeptide chains and a single structural model, i.e. files that do not contain MODEL/ENDMDL tags. 

All PDB files used for testing purposes of the script are available as a .zip file in this repository.

## Main Features

### (1) Helix Tool
* Determines the number of alpha-helices in a protein.
* Computes helixID, length of each helix, and the mean length of a helix.

### (2) Sheet Tool
* Determines the number of beta-sheets in a protein.
* For each sheet, determines the sheetID, number of strands and total number of residues.
* Computes the mean number of strands and residues in a sheet.

### (3) ATOM Entry Parser
* Calculates the number of polypeptide chains in a protein and the length of each chain.
* Computes the overall percentage of residues in a helix, in a sheet and in neither structure.

### (4) Residue Distribution Calculator
* Accepts a single PDB file.
* Calculates the frequency of each amino acid in a protein sequence.

### (5) Global Residue Distribution Calculator
* Accepts multiple PDB files.
* Calculates the frequency of each amino acid across all provided files.

### (6) Secondary Structure Propensity Calculator
* Calculates the residue frequency separately for helical, sheet and irregular segments.
* Computes the propensity of each residue to be in each type of secondary structure.

### (7) Barchart Generator Tool
* Produces barchart plots of the residue distribution in a single protein (4) and the residue propensities (6).

## Pre-requisites
Python 2 or 3
