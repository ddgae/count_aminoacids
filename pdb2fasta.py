# -*- coding: utf-8 -*-
from __future__ import division
import sys
import re
import itertools
import numpy as np
import math
# author:david h. gae
# copyright © 2021 david gae Some right reserved

if len(sys.argv) <= 1:
    print ('usage: python pdb2fasta.py file.pdb (RNA concentration: integer) (concentration (<0.49 M): integer) (ion species: + integer - integer) > file.fasta')
    exit()
# Open PDB file
# Python: readline() read a single line as a string, read() read is n bytes (8 bits), readlines() list of strings. sys.argv(): first command is read as a string.
input_file1 = open(sys.argv[1])
input_file2 = open(sys.argv[1])
input_file3 = sys.argv[1]
f = open(input_file3,'r')
lines = f.readlines()
f.close()

# Open write file1
name=(sys.argv[1].split('.',1)[0])
f1= open(name+'.fasta',"w+")
f1.write(">P1;seq" + '\n')
letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','TER':'/'}

# Open write file2
name=(sys.argv[1].split('.',1)[0])
f2= open(name+'.fasta1',"w+")
f2.write(">P1,seq"+'\n')

# Empty Array
line3 = []
line4 = []
line5 = []
line2 = []
line6 = []
# Dictionary of amino acids and list based on dictionary.
aa = list(letters.keys())
aa2 = list(letters.values())
r = list()
z = list()

for line in input_file1:
    #print(line)
    if line[0:3] == 'TER' and line[13:13] == '' or line[13:14] == 'N' and line[13:16] != 'ND1' and line[13:16] != 'NE2' and line[13:15] != 'NE' and line[13:16] != 'NH1' and line[13:16] != 'NH2' and line[13:15] !='NZ':
       line3.append(line)
       # print(line)
       for i in line3:
                # print(i)
                # convert amino acid three letter to single letter
                for m in range(21):
                    # print(m)
                    if i[17:20] == aa[m]:
                       r = aa2[m]
                       z.append(r)
       line4 = ("".join(r))
       f1.write(line4)
f1.close()


# Biology: Sequences
# Python: Can not use the same file, unless you define a readlines() string and store it.
for line6 in lines:
    # print(line6)
    if line6[17:20] == 'ASP' and line6[13:14] == 'N' or line6[17:20] == 'GLU' and line6[13:14] == 'N' or line6[17:20] == 'HIS' and line6[13:14] == 'N' and line6[13:16] != 'ND1' and line6[13:16] != 'NE2' or line6[17:20] == 'ARG' and line6[13:14] == 'N' and  line6[13:15] != 'NE' and line6[13:16] != 'NH1' and line6[13:16] != 'NH2' or line6[17:20] == 'LYS' and line6[13:14] == 'N' and line6[13:15] != 'NZ':
        line3.append(line6)
        for i in line3:
            # print(i)
            # convert amino acid three letter to single letter
            for m in range(21):
                    # print(m)
                if i[17:20] == aa[m]:
                    r = aa2[m]
                    z.append(r)
        line5 = ("".join(r))
        f2.write(line5)
f2.close()
negative = 0
positive = 0

# Chemistry: Determination of Ionic molality
# 	     some formulation is defined on papers, and research findings.
y = list()
ya = list()
for line7 in input_file2:
    # 5/20 are charged amino acids
    if (line7[17:20] == 'ASP' and line7[13:14] == 'N' or line7[17:20] == 'GLU' and line7[13:14] == 'N'):
        negative= negative + 1
        #print(negative)
    #if (line7[17:20] == 'GLU' and line7[13:14] == 'N'):
    #    negative = negative + 1
    if (line7[17:20] == 'HIS' and line7[13:14] == 'N' and line7[13:16] != 'ND1' and line7[13:16] != 'NE2'  or
        line7[17:20] == 'LYS' and line7[13:14] == 'N' and  line7[13:15] != 'NZ' or
        line7[17:20] == 'ARG' and line7[13:14] == 'N' and  line7[13:15] !='NE' and line7[13:16] !='NH1' and line7[13:16] != 'NH2'):
        positive = positive +1
        #print(positive)
    #if (line7[17:20] == 'ARG' and line7[13:14] == 'N' and line7[13:16] != 'NE' and line7[13:16] != 'NH1' and line7[13:16] != 'NH2'):
    #    positive = positive +1
    #if (line7[17:20] == 'LYS' and line7[13:14] == 'N' and line7[13:16] != 'NZ'):
    #    positive = positive +1
print('- amino acid:', negative, '+ amino acid:', positive)
# Chemistry: Formal RNA properties can be found in Review Papers in bulk properties. 
# 	     Standard Room Temperature conditions are in 298K and dielectric of 78.54
# 	     The assumptions are in bulk solvent with NaCL (58 molarity in 1 L can be converted to 0.085 molality)
# 	     The equation are found in Chapter 5.8 of Cheng Raymond, Physical chemistry for the bioscience,sausalito CA, university science books,2005
#	     Mean activity coefficient = -0.509|[z+][z-]|Sqrt(I)  at standard room temperature 298K in bulk water. I = 1/2  Σ sqrt(M*(z+))^2
# 	     Result are for NaCl only.
#	     y = (-0.509 * (1) * 0.5) * float(math.sqrt((mol/kg*(1)**2))
y = (-0.509 * (-1)*(1)) * 0.5 * float(math.sqrt((0.085)*(1)**2))
# Math Property: log of negative is NaN
y1 = 10**y
print("mean activity coefficient with 0.5M y+/-:", y1)
# Hypothesis, Mean ionic activity with RNA may be a constant value.
# PLEASE NOTE addition of constant value is a HYPOTHESIS PLEASE DO NOT TAKE IT AS IS:
ya = (-0.509 * (1)* (-1)) * 0.5 * float(math.sqrt((0.085)*(1)**2)) + int(sys.argv[2])/8500
y2 = 10**ya
print("y+/-+ presence of RNA (hypothesis of RNA concentration in bulk water):",y2)
# Possible ionic strength of buffer concentration for PCR-based assay.
# Reverberi et al. Factor affecting the antigen-antibody reaction. Blood Transfusion 2007 Oct 5(4) 227-240
# I = 1/2 *(n stoichiometry) * Σ (Molarity (i)) * v (i) ^ 2)
# v = Na+ (1+) and Cl- (1-) ( n = 1), Note this is only for NaCl
# SIMPLE CALCULATION EXAMPLE:
#     Typically 250 ml breaker may be used for making a buffer solution:  58 g/mol * 0.02 mol/L * L = 1.16 gram, 5 g/250 mL,  = 0.02 g/ml
#     Remember that the standard volume for buffer is about 500 ml or PCR reservoir, 100/50 ml
#     Calculation: 58 g/mol molecular weight of NaCl * 1 mol/L = 58 g/L
#    		   34 plates (40 rxn extra = 0.01 L or 10 ml) * 96 well plate per 200 ul rxn) = 0.6628 L
#     		   58 g/mol * 0.02(mol/L) * 0.6628 L = X g (Volume required use for 96 well plate, if master mixature is purchased this may be the volume,
# 							    I might be wrong).
I = (0.5 * (1))  * float(sys.argv[3])* int(sys.argv[4])** 2 + float(sys.argv[3])* int(sys.argv[5]) ** 2
#I = (0.5 * (1)) * float((sys.argv[2])) (sys.argv[3])*(sys.argv[4])
print("Ionic strength:", I)
# Physical Chemistry: Debye-Hukel equation
#	   Constants: 1.824e10^6/ (eT)^3/2 |z+z-|sqrt(I)
# Ionic strength, univeristy of michigan chem241
# Formula:	     -0.509 *|z+z-| (sqrt(u))/ 1 + ( alpha (sqrt(u)/305))
# Cheng Raymond, Physical chemistry for the bioscience,sausalito CA, university science books,2005, chapter 5.11
# 		     -0.509 * z^2*sqrt(I)
# 		      alpha = hydration radius of ion since Na+ Cl- dissipate in Bulk water: ~0.138 nm or 138 pm
#ions +
A = sys.argv[4]
#ions -
B = sys.argv[5]
A2 = abs(int(A) ** 2)
B2 = abs(int(B) ** 2)
Ey = -0.509 * A2 * B2 * math.sqrt(I)
#Ey = -0.509 * float(sys.argv[3])*float(sys.argv[2])+float(sys.argv[4])*float(sys.argv[2])
Ey2 = 10 ** Ey
print("activity cofficient Ey+/-", Ey2)
#Chemistry:	      pH value related to ionic concentration
# 		      This theory hold true in low ionic concentration.
# 		      0.02 M is standard buffer concentration.
C2 = sys.argv[3]
#print(C2)
pH =  (I * float(C2))
#print(pH)
pH2 = -(math.log10(pH))
# Math Property: log of negative is NaN
if (float(pH) > 0 and float(pH2) > 0.4):
    print("pH is dependent on concentration please check script comments:", round(pH2))
