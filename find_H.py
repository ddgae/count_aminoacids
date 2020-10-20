import sys
import re
import itertools
import numpy as np
import subprocess
from subprocess import Popen, PIPE
#author:gae_d

if len(sys.argv) <= 1:
    print ('usage: python3 pdb2fasta.py file.pdb > file.fasta')
    exit()
#open of pdb file
input_file = open(sys.argv[1])
#new = sys.argv[2]
#open of write file
name=(sys.argv[1].split('.',1)[0])
f= open(name+'.fasta',"w+")
f.write(">P1;seq" + '\n')
f.write("sequence:seq:  :: ::::-1.00:-1.00" + '\n')
letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
           'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
           'TYR':'Y','VAL':'V','TER':'/'}

#empty char array of zeros.
line3 = []
line4 = []
line2 = []
#dicitionary table of AA
# define list class baesd on dictionary.
aa = list(letters.keys())
aa2 = list(letters.values())
r = list()
z = list()

for line in input_file:
    if line[0:3] == 'TER' or line[13:14] == 'N' and line[13:15] != 'NE' and line[13:15] != 'NH1' and line[13:15] != 'NE2' and line[13:15] != 'NH2' and line[13:15] != 'NZ' and line[13:15] !='NH' and line[13:15] !='ND':
       line3.append(line)
       #sequence starts with 1 not, zero
       i = 1
       for i,x in enumerate(line3):
           for m in range(21):
              #for each character value add '/' and append
              if line[17:20] == aa[m]:
                  r = aa2[m]
                  r = '/'
                  z.append(r)
              #for each character substitute amino acid  with 1-letter code
              #elif line[17:18] == sys.argv[2]:
              if line[17:20] == sys.argv[2]:
                  r = 'H'
                  #r = aa2[m]
                  #print(r)
       line4 = (" ".join(r))
       f.write(line4)
f.close()


