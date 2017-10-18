#this python code tests generalized_birkhoff_von_neumann_decomposition.py

import matplotlib
import ggplot
import numpy as np
import generalized_birkhoff_von_neumann
from generalized_birkhoff_von_neumann import bihierarchy_test
from generalized_birkhoff_von_neumann import generalized_birkhoff_von_neumann_decomposition as decomp

#we first test starting matrices and trivial constraint structures (allowing each entry of basis
#matrices to be zero or one)

for n in range (2,3):
  for m in range (2,3):
    X = np.random.uniform(0, 1, (n,m))
    if decomp(X,{})[-1] == X:
      print("success")
    else:
      print("failure")
    








