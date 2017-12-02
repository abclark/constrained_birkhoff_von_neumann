R = {0: [0,1,2], 1: [2,1,0], 3: [2,0,1]} 
 
import numpy as np 
 
m = np.ones(len(R[0])) 
P = {0: [0,0,0], 1: [0,0,0], 3: [0,0,0]} 
 
while max(m) > 0: 
  for key,value in R.items(): 
    R[key] = [i for i in value if m[i] != 0] 
  y = np.zeros(len(m)) 
  for key,value in R.items(): 
    y[value[0]] += 1 
  #time taken to deplete remaining mass 
  z = [max(i,0.001)/max(j,0.0001) for i,j in zip(m,y)] 
  #reduce masses  
  for key,value in R.items(): 
    m[value[0]] -= min(z) 
    P[key][value[0]] += min(z)           
else: 
  print(P) 
 
