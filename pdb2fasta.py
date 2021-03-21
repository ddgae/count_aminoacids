from __future__ import division 
import sys
import re
import itertools
import numpy as np
import math

#author:david gae
#copyright © 2021 david gae Some right reserved.

if len(sys.argv) <= 1:
    print ('usage: python pdb2fasta.py file.pdb (RNA concentration: integer)  > file.fasta')
    exit()
#open of pdb file
#if use file then iterate only once
input_file = sys.argv[1]
f = open(input_file,'r')
lines = f.readlines()
f.close()

#open of write file
name=(sys.argv[1].split('.',1)[0])
f1= open(name+'.fasta',"w+")
f1.write(">seq" + '\n')
#f1.write("sequence:seq:  :: ::::-1.00:-1.00" + '\n')
letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','TER':'/'}

name=(sys.argv[1].split('.',1)[0])
f2= open(name+'.fasta1',"w+")
f2.write(">seq"+'\n')
#f2.write("sequence:seq:  :: ::::-1.00:-1.00" + '\n')

#empty char array of zeros.
line3 = []
line4 = []
line5 = []
line2 = []
line6 = []
#dicitionary table of AA
# define list class baesd on dictionary.
aa = list(letters.keys())
aa2 = list(letters.values())
r = list()
z = list()

for line in lines:
    if line[0:3] == 'TER' or line[13:14] == 'N' and line[13:15] != 'NE' and line[13:15] != 'NH1' and line[13:15] != 'NE2' and line[13:15] != 'NH2' and line[13:15] != 'NZ' and line[13:15] !='NH' and line[13:15] !='ND':
       line3.append(line)
       #print(line)
       for i in line3:
                #print(i)
                #convert amino acid three letter to single letter
                for m in range(21):
                    #print(m)
                    if i[17:20] == aa[m]:
                       r = aa2[m]
                       z.append(r)
       line4 = ("".join(r))
       f1.write(line4)
f1.close()


#charge residues determination.
#Assign new loop for a new file, [can not use loop again using the same file].
for line6 in lines:
    #print(line6)
    if line6[17:20] == 'ASP' and line6[13:14] == 'N' \
            or line6[17:20] == 'GLU' and line6[13:14] == 'N' \
            or line6[17:20] == 'HIS' and line6[13:14] == 'N' and line6[13:16] != 'ND1' and line6[13:16] != 'NE2' \
            or line6[17:20] == 'ARG' and line6[13:14] == 'N' and  line6[13:15] != 'NE' and line6[13:16] != 'NH1' and line6[13:16] != 'NH2' \
            or line6[17:20] == 'LYS' and line6[13:14] == 'N' and line6[13:15] != 'NZ':
        line3.append(line6)
        for i in line3:
            #print(i)
            #convert amino acid three letter to single letter
            for m in range(21):
                    #print(m)
                if i[17:20] == aa[m]:
                    r = aa2[m]
                    z.append(r)
        line5 = ("".join(r))
        f2.write(line5)
f2.close()
negative = 0
positive = 0

#Ionic molality
y = list()
ya = list()
for line7 in lines:
    # 5/20 is charged amino acid
    if (line7[17:20] == 'ASP' and line7[13:14] == 'N'):
        negative= negative + 1
    if (line7[17:20] == 'GLU' and line7[13:14] == 'N'):
        negative = negative + 1
    if (line7[17:20] == 'HIS' and line7[13:14] == 'N' and line7[13:16] != 'ND1' and line7[13:16] != 'NE2'):
        positive = positive +1
    if (line7[17:20] == 'ARG' and line7[13:14] == 'N' and line7[13:16] != 'NE' and line7[13:16] != 'NH1' and line7[13:16] != 'NH2'):
        positive = positive +1
    if (line7[17:20] == 'LYS' and line7[13:14] == 'N' and line7[13:16] != 'NZ'):
        positive = positive +1
print('- amino acid', negative, '+ amino acid', positive)
	#RNA formal charge
	#Mean ionic activity at 298K and dielectric of 78.54
	#The assumption is that it is bulk solvent in NaCL (58 molarity in 1 L therefore it is 0.058 molality)
	#The equation was found in Chapter 5.8 
	#Cheng Raymond, Physical chemistry for the bioscience,sausalito CA, university science books,2005
	#Mean ionic activity  = -0.509|[z+][z-]|Sqrt(I)  at standard 298K in bulk water. I = 1/2 sqrt(M*(z+))^2
	#Determined for NaCl only. 
y = (-0.509 * (1) * 0.5) * float(math.sqrt((58/1000)*(1)**2))
	#log of negative is NaN
y1 = 10**y
print("mean activity coefficient y+/-", y1)
	#Possible purposed hypothesis if an RNA was a constant.
	# PLEASE NOTE THIS EQUATION IS A HYPOTHESIS, DO NOT TAKE IT AS IS:
	# Constant maybe the number of charges of RNA.
ya = (-0.509 * (1) * 0.5) * float(math.sqrt((58/1000)*(1)**2)) + int(sys.argv[2])
y2 = 10**ya
print("y+/-+ presence of RNA (guess)",y2)
