# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 10:52:48 2021

@author: Sam
"""

from dtw import *
from timeit import default_timer as timer
import numpy as np

with open('validation_signals') as f:
    lines = f.readlines()
    
X = []
for i in range(0,4):
    X.append([float(x) for x in lines[i].split()])
 
#commented patterns not recognized despite being listed in stepPattern.py:    
 
patterns = [
    'symmetric1'                ,
    'symmetric2'                ,
    'asymmetric'                ,
 #   'symmetricVelichkoZagoruyko',
 #   'asymmetricItakura'         ,
    'symmetricP0'               ,
    'asymmetricP0'              ,
 #   'asymmetricP0b'             ,
    'symmetricP05'              ,
    'asymmetricP05'             ,
    'symmetricP1'               ,
    'asymmetricP1'              ,
    'symmetricP2'               ,
    'asymmetricP2'              ,
    'typeIa'                    ,
    'typeIb'                    ,
    'typeIc'                    ,
    'typeId'                    ,
    'typeIas'                   ,
    'typeIbs'                   ,
    'typeIcs'                   ,
    'typeIds'                   ,
    'typeIIa'                   ,
    'typeIIb'                   ,
    'typeIIc'                   ,
    'typeIId'                   ,
    'typeIIIc'                  ,
    'typeIVc'                   ,
    'mori2006'                  ]
    
start = timer()

warp = []
for i in range(0,len(patterns)-1):
    warp.append(dtw(X[0],X[1],step_pattern=patterns[i]))

end = timer()
print(end-start)

textfile = open("validation_dtwpython_output","w")
for i in range(0,len(patterns)-1):
    index1s_str = " ".join(str(j) for j in warp[i].index1s)
    index2s_str = " ".join(str(j) for j in warp[i].index2s)
    textfile.write(patterns[i] + " " + index1s_str + "\n")
    textfile.write(patterns[i] + " " + index2s_str + "\n")
    
#validation of open end/begin
warpo=[]

patternso = ['symmetric2','asymmetric']

for i in (0,1):
    for j in (0,1):
        warpo.append(dtw(X[2],X[3],
                         open_begin=i,open_end=j,step_pattern=patternso[i]))

    
    


