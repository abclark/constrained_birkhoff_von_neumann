#generalized_birkhoff_von_neumann.py decomposes a matrix into a weighted average of basis matrices with integer coefficients
#satisfying imposed constraints. When the starting matrix is doubly stochastic and the basis matrices are restricted to
#be permutation matrices, this is the classical birkhoff_von_neumann decomposition. Formally, we implement the algorithm 
#identified in Budish, Che, Kojima, and Milgrom (2013). Thus, the constraint structure must form what they call a bihierarchy.
#
# Copyright 2017 Aubrey Clark.
#
# generalized_birkhoff_von_neumann is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

#: The current version of this package.
__version__ = '0.0.1-dev'

import networkx as nx
import numpy as np
import copy
import itertools
import math
from pprint import pprint

#numbers in some arrays will be rounded to the nearest integer if within the following tolerance
tolerance = np.finfo(np.float).eps * 10

#Below are two example startitng matrices X and two corresponding constraint structures. X is the matrix we wish to decompose 
#into a weighted average of basis matrices. constraint_structure is a dictionary whose keys are subsets of coordinates of the #basis matrices (the dimensions of which are the same as X) (e.g. frozenset({(0, 0), (0, 1), (0,2)})) refers to the (0,0),
#(0,1),(0,2) coordinates), and whose keys refer to the minimum and maximum sum of these entries in each of the basis matrices 
#(e.g. the value (1,1) means that the coordinates that this value's key represents sum to exactly one in each 
#of the basis matrices.)

#X = np.array([[.5, .2,.3], [.5,.5, 0], [.8, 0, .2], [.2, .3, .5]])
#constraint_structure = {frozenset({(0, 0), (0, 1), (0,2)}): (1,1), frozenset({(1, 0), (1, 1), (1,2)}):(1,1), frozenset({(2, 0), (2, 1), (2,2)}):(1,1), frozenset({(3, 0), (3, 1), (3,2)}):(1,1), frozenset({(0, 0), (1, 0), (2,0), (3,0)}):(1,2),  frozenset({(0, 1), (1, 1), (2,1), (3,1)}):(1,1), frozenset({(0, 2), (1, 2), (2,2), (3,2)}):(1,1), frozenset({(0, 0), (1, 0)}):(1,1)}

#X = np.array([[.3, .7], [.7,.3]])
#constraint_structure = {frozenset({(0, 1),(1,0)}): (1,2), frozenset({(1, 0),(1,1)}): (1,1)}

#bihierarchy_test decomposes the constraint structure into a bihierarchy, if it is one. If the constraint structure is not
#a bihierarchy or the starting matrix X does not abide by the constraint structure to begin with, bihierarchy_test will tell 
#you. bihierarcy_test will be invoked by generalized_berkhoff_con_neumann_decomposition, and so the latter function
#(which performs the decomposition) will warn about these issues.

def bihierarchy_test(X, constraint_structure):
  S = {index for index, x in np.ndenumerate(X)}
  if any(X[i]<0 or X[i]>1 for i in S):
    print("matrix entries must be between zero and one")
  for key, value in constraint_structure.items():
    if sum([X[i] for i in key]) < value[0] or sum([X[i] for i in key]) > value[1]:
      print("impossible constraint structure capacities")
  C = []
  for key, value in constraint_structure.items():
    C.append(set(key))
  permutations =  itertools.permutations(C)
  for c in permutations:
    listofA, listofB = [], []
    for idx, x in enumerate(c):
      if all( x < y or y < x or x.isdisjoint(y) for y in [c[i] for i in listofA]):
        target = listofA
      elif all(x < y or y < x or x.isdisjoint(y) for y in [c[i] for i in listofB]):
        target = listofB
      else:
        break
      target.append(idx)
    if len(listofA) + len(listofB) == len(c):
      return [[c[i] for i in listofA], [c[i] for i in listofB]]
    else:
      print("constraint structure cannot be decomposed into a bihierarchy")

#generalized_birkhoff_von_neumann_iterator is the core step in the decomposition. After the starting matrix X and the
#constraint structure have been represented as a weighted, directed graph G, this function takes as input a list H = [(G,p)]
#(where p is a probability, initially one) and decomposes the graph into two graphs, 
#each with an associated probability, and each of which are closer to representing a basis matrix.
#Seqential iteration, done in the main function generalized_birkhoff_von_neumann_decomposition, leads to the
#decomposition.
      
def generalized_birkhoff_von_neumann_iterator(X,H):
  tolerance = np.finfo(np.float).eps * 10
  S = {index for index, x in np.ndenumerate(X)}
  (G, p) = H.pop(0)
  #
  #remove edges with integer weights
  #extracts all edges satisfy the weight threshold:
  #
  eligible_edges = [(from_node,to_node,edge_attributes) for from_node,to_node,edge_attributes in G.edges(data=True) if from_node == frozenset(S) or to_node == frozenset(S) or all(i+tolerance < edge_attributes['weight'] or edge_attributes['weight'] < i - tolerance for i in range(0,math.floor(sum(sum(X) ) ) ) ) ]
  #
  K = nx.DiGraph()
  K.add_edges_from(eligible_edges)
  #
  #find a cycle and compute the push_forward and push_reverse probabilities and graphs
  cycle = nx.find_cycle(K, orientation='ignore')
  forward_weights = [(d['weight'],d['min_capacity'],d['max_capacity']) for (u,v,d) in K.edges(data=True) if (u,v,'forward') in cycle]
  reverse_weights = [(d['weight'],d['min_capacity'],d['max_capacity']) for (u,v,d) in K.edges(data=True) if (u,v,'reverse') in cycle]
  push_forward = min((x[2] - x[0] for x in forward_weights))
  push_reverse = min((x[2] - x[0] for x in reverse_weights))
  pull_forward = min((x[0] - x[1] for x in forward_weights))
  pull_reverse = min((x[0] - x[1] for x in reverse_weights))
  push_forward_pull_reverse = min(push_forward,pull_reverse)
  push_reverse_pull_forward = min(pull_forward,push_reverse)
  #
  #Construct the push_forward_pull_reverse graph
  #
  G1 = copy.deepcopy(G)
  for u,v,d in G1.edges(data=True):
    if (u,v,'forward') in cycle:
      d['weight']+=push_forward_pull_reverse
    if (u,v,'reverse') in cycle:
      d['weight']+=-push_forward_pull_reverse
  #
  #Construct the push_reverse_pull_forward graph
  #
  G2 = copy.deepcopy(G)
  for u,v,d in G2.edges(data=True):
    if (u,v,'reverse') in cycle:
      d['weight']+=push_reverse_pull_forward
    if (u,v,'forward') in cycle:
      d['weight']+=-push_reverse_pull_forward
  #
  gamma = push_reverse_pull_forward/(push_forward_pull_reverse + push_reverse_pull_forward)
  return([(G1,p*gamma), (G2,p*(1-gamma))])

#generalized_birkhoff_von_neumann_decomposition takes the primitives, a starting matrix X (a numpy array) and
#a constraint structure constraint_structure (a dictionary, whose structure is described in the examples above)
# and outputs a list of coefficients and basis matrices (each satisfying the constraints) such that the starting matrix is a
#weighted sum of the basis matrices. It first applies bihierarchy_test to these primitives (see above for what this does).
#Then, if all is okay at the first step, it represents these primitives as a weighted, directed graph and iteratively
# applies generalized_birkhoff_von_neumann_iterator. Finally, it cleans the solution, transforming the
#final iteration directed, weighted graphs to basis matrices, merging duplicate basis matrices and their probabilities,
#checking that the probabilities form a distribution, and checking that the average of the basis matricies under this
#distribution is indeed the starting matrix X

def generalized_birkhoff_von_neumann_decomposition(X,constraint_structure):
  #
  tolerance = np.finfo(np.float).eps * 10
  S = {index for index, x in np.ndenumerate(X)}
  #
  A,B = bihierarchy_test(X, constraint_structure)
  A.append(S), B.append(S)
  #
  for x in S:
    A.append({x}), B.append({x})
  #
  for x in S:
    constraint_structure.update({frozenset({x}):(0,1)})
  #
  R1 = nx.DiGraph()
  for x in A:
    for y in A:
      if x < y and not any(x < z < y for z in A):
        R1.add_edge(frozenset(y),frozenset(x),weight=sum([X[i] for i in x]), min_capacity = constraint_structure[frozenset(x)][0], max_capacity = constraint_structure[frozenset(x)][1])
  #
  R2 = nx.DiGraph()
  for x in B:
    for y in B:
      if y < x and not any(y < z < x for z in B):
        R2.add_edge((frozenset(y),'p'),(frozenset(x),'p'),weight = sum( [X[i] for i in y]), min_capacity = constraint_structure[frozenset(y)][0], max_capacity = constraint_structure[frozenset(y)][1])
  #
  G = nx.compose(R1,R2)
  #
  for index, x in np.ndenumerate(X):
    G.add_edge(frozenset({index}), (frozenset({index}),'p'), weight=x, min_capacity = 0, max_capacity = 1)
  #
  H=[(G,1)]
  solution=[]
  #
  while len(H) > 0:
    if any(tolerance < x < 1 - tolerance for x in [d['weight'] for (u,v,d) in H[0][0].edges(data=True) if u in [frozenset({x}) for x in S]]):
      H.extend(generalized_birkhoff_von_neumann_iterator(X,[H.pop(0)]))
    else:
      solution.append(H.pop(0))
  #
  solution_columns_and_probs = []
  #
  for y in solution:
    solution_columns_and_probs.append([[(u,d['weight']) for (u,v,d) in y[0].edges(data=True) if u in [frozenset({x}) for x in S]],y[1]])
  #
  solution_zeroed = []
  #
  for z in solution_columns_and_probs:
    list = []
    for y in z[0]:
      if y[1] < tolerance:
        list.append((y[0],0))
      elif y[1] > 1-tolerance:
        list.append((y[0],1))
    solution_zeroed.append([list,z[1]])
  #
  assgs = []
  coeffs = []
  #
  for a in solution_zeroed:
    Y = np.zeros(X.shape)
    for x in a[0]:
      for y in x[0]:
        Y[y]=x[1]
    assgs.append(Y)
    coeffs.append(a[1])
  #
  list = []
  #
  for idx, x in enumerate(solution_zeroed):
    if all(x[0]!= z[0] for z in [solution_zeroed[i] for i in list]):
      list.append(idx)
  #
  solution_simplified = []
  #
  for i in list:
    solution_simplified.append([solution_zeroed[i][0],sum([x[1] for x in solution_zeroed if x[0]==solution_zeroed[i][0]])])
  #
  assignments = []
  coefficients = []
  #
  for a in solution_simplified:
    Y = np.zeros(X.shape)
    for x in a[0]:
      for y in x[0]:
        Y[y]=x[1]
    assignments.append(Y)
    coefficients.append(a[1])
  #
  return([coefficients, assignments, sum(coefficients), sum(i[1]*i[0] for i in zip(coefficients, assignments))])
  pprint([coefficients, assignments, sum(coefficients), sum(i[1]*i[0] for i in zip(coefficients, assignments))])

