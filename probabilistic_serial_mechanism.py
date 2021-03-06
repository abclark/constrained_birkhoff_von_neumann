"""
=======================================
probabilistic_serial_mechanism.py
=======================================
Runs the probabilistic serial mechanism, taking as input a dictionary R. 
A key of this dictionary corresponds to an agent,
a value the agent's rank order list of the objects. 
Output is a dictionary P. 

A value of the output dictionary is a list whose entries are probabilities that the agent wins the 
corresponding objects.
E.g. for the dictionary R = {0: [0,1,2], 1: [2,1,0]}, there are two agents and three objects. 
Agent 1 ranks object 2 above object 1 above object 0. 
Output for this dictionary is P = {0: array([ 1. ,  0.5,  0. ]), 1: array([ 0. ,  0.5,  1. ])}. 
Agent 0 wins object 0 with probability 1 and object 1 with probability 0.5.

Copyright 2017 Aubrey Clark.

probabilistic_serial_mechanism.py is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the Free Software Foundation, 
either version 3 of the License, or (at your option) any later version.

"""

import numpy as np
import copy

TOLERANCE = 1.0e-5

def probabilistic_serial_mechanism(R, m):
  # define an empty dictionary to store the solution, and copies of R and m
  P={}
  Q = copy.deepcopy(R)
  t = copy.deepcopy(m)
  # give the empty dictionary P the structure of the solution
  for key, value in Q.items():
    P[key]= [0]*len(value)
  # eat probability mass while it remains
  while any(Q[key] != [] for key,values in Q.items()):
    # if an object has no remaining probability mass, remove it from the rank order lists
    for key,value in Q.items():
      Q[key] = [i for i in value if t[i] > TOLERANCE and P[key][i] < 1-TOLERANCE]
    # define a zero vector whose dimension equals the number of objects
    y = np.zeros(len(t))
    # define a vector of ones, one for each agent
    x = np.ones(len(Q.items()))
    # count how many agents rank each object first (under updated rank order lists)
    # for each agent's preferred object (under updated rank order lists), 
    # record in x how much more of that object's probability mass they can consume
    for key,value in Q.items():
      if value != []:
        y[value[0]] += 1
        x[key] = 1 - P[key][value[0]]
    # define a vector recording time taken until each object's probability mass is depleted without intervention
    z = [max(i,0.000001)/max(j,0.0000001) for i,j in zip(t,y)]
    # update probability masses m and record probability consumed in P
    for key,value in Q.items():
      if value != []:
        t[value[0]] -= min(min(z), min(x))
        P[key][value[0]] += min(min(z), min(x))
  # if all probaility masses are nil, the process is done—return the solution
  else:
    print([P,np.array([value for key, value in P.items()])])
