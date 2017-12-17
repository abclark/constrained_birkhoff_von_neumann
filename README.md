# generalized_birkhoff_von_neumann


<p>Decomposes a matrix into a weighted sum of basis matrices with binary entries satisfying user-imposed constraints. When the starting matrix is doubly stochastic and the basis elements are required to be permutation matrices, this is the classical Birkhoff von-Neumann decomposition.
Here we implement the algorithm identified in <a href="http://faculty.chicagobooth.edu/eric.budish/research/Budish-Che-Kojima-Milgrom-2013-AER.pdf">Budish, Che, Kojima, and Milgrom 2013</a>. 
The constraints must form what they call a bihierarchy.</p>

<h2>Installation</h2>
        
<pre><code>pip install generalized_birkhoff_von_neumann</code></pre>
        
<h2>Basic usage</h2>
        
<pre><code>import numpy as np
from generalized_birkhoff_von_neumann import generalized_birkhoff_von_neumann_decomposition

# Create a matrix whose entries are between 0 and 1, and a constraint structure. 
#
#For example
#
# X = np.array([[.5, .2,.3], [.5,.5, 0], [.8, 0, .2], [.2, .3, .5]])
# 
# constraint_structure = {frozenset({(0, 0), (0, 1), (0,2)}): (1,3)}
#
# The first, and only, entry of this constraint structure requires that the 
# (0, 0), (0, 1), (0,2) entries of each basis matrix sum to 1, 2, or 3.

generalized_birkhoff_von_neumann_decomposition(X, constraint_structure)
</code></pre>
        
<h2>Mathematical background</h2>
       
<p>See <a href="http://faculty.chicagobooth.edu/eric.budish/research/Budish-Che-Kojima-Milgrom-2013-AER.pdf">Budish, Che, Kojima, and Milgrom 2013</a>.</p>

<h2>API</h2>
        
<pre><code>>>> import numpy as np
>>> from generalized_birkhoff_von_neumann import generalized_birkhoff_von_neumann_decomposition
>>> X = np.array([[.5, .2,.3], [.5,.5, 0], [.8, 0, .2], [.2, .3, .5]])
>>> constraint_structure = {frozenset({(0, 0), (0, 1), (0,2)}): (1,1), frozenset({(1, 0), (1, 1), (1,2)}):(1,1), frozenset({(2, 0), (2, 1), (2,2)}):(1,1), frozenset({(3, 0), (3, 1), (3,2)}):(1,1), frozenset({(0, 0), (1, 0), (2,0), (3,0)}):(1,2), frozenset({(0, 1), (1, 1), (2,1), (3,1)}):(1,1), frozenset({(0, 2), (1, 2), (2,2), (3,2)}):(1,1), frozenset({(0, 0), (1, 0)}):(1,1)}
>>> generalized_birkhoff_von_neumann_decomposition(X,constraint_structure)
[[0.14285714285714288, 
0.3571428571428571, 
0.29999999999999999, 
0.057142857142857141, 
0.14285714285714282], 
[array([[ 0.,  1.,  0.],
       [ 1.,  0.,  0.],
       [ 1.,  0.,  0.],
       [ 0.,  0.,  1.]]), 
array([[ 1.,  0.,  0.],
       [ 0.,  1.,  0.],
       [ 1.,  0.,  0.],
       [ 0.,  0.,  1.]]), 
array([[ 0.,  0.,  1.],
       [ 1.,  0.,  0.],
       [ 1.,  0.,  0.],
       [ 0.,  1.,  0.]]), 
array([[ 0.,  1.,  0.],
       [ 1.,  0.,  0.],
       [ 0.,  0.,  1.],
       [ 1.,  0.,  0.]]), 
array([[ 1.,  0.,  0.],
       [ 0.,  1.,  0.],
       [ 0.,  0.,  1.],
       [ 1.,  0.,  0.]])], 
1.0, 
array([[ 0.5,  0.2,  0.3],
       [ 0.5,  0.5,  0. ],
       [ 0.8,  0. ,  0.2],
       [ 0.2,  0.3,  0.5]])]
>>> </code></pre>

<h2>Application: probabilistic serial mechanism with constraints</h2>        
        
<p>There are <code>n</code> objects to be allocated among <code>m</code> agents. Each agent has a ranking of the objects. The probability of agent <code>i</code> winning object <code>j</code> is derived from the following process: agents simultaneously and at the same speed 'eat' probability mass of their top ranked object for which probability mass still remains. The resulting probabilities can be collected in an <code>m</code> times <code>n</code> matrix <code>X</code> whose entry in its <code>i</code>'th row and <code>j</code>'th column is the probability that agent <code>i</code> wins object <code>j</code>.</p>
        
<p>To implement the probabilities of <code>X</code>, one needs to find an appropriate probability distribution over allocations (full descriptions of who gets what). An allocation may be represented as a matrix with binary entries, a <code>1</code> in the <code>i</code>'th row and <code>j</code>'th column meaning that agent <code>i</code> wins object <code>j</code> and a <code>0</code> meaning that agent <code>i</code> does not win object <code>j</code>. If <code>n</code> equals <code>m</code>, such a matrix is a permutation matrix. Further, <code>X</code> will necessarily be doubly stochastic and so the problem may be solved using the classical Birkhoff-von Neumann decomposition.</p>
  
<p><code>generalized_birkhoff_von_neumann_decomposition</code> may be used when <code>n</code> not equal to <code>m</code>, or when one wishes to impose additional constraints (e.g. each agent is allocated at least one object, no agent is allocated more than <code>k</code> objects, etc).</p>

<p>To run the probabilistic serial mechanism, follow the instructions below.</p> 
        
<h3>Installation</h3>
        
<pre><code>pip install probabilistic_serial_mechanism</code></pre>
        
<h3>Basic usage</h3>
        
<pre><code># Create a dictionary R with agents as keys and rank order lists as values.
#
#For example:
#
# R = {0: [0,1,2], 1: [2,1,0]}

probabilistic_serial_mechanism(R)</code></pre>
        
<h3>Mathematical background</h3>
       
<p>See <a href="http://faculty.chicagobooth.edu/eric.budish/research/Budish-Che-Kojima-Milgrom-2013-AER.pdf">Budish, Che, Kojima, and Milgrom 2013</a>.</p>

<h3>API</h3>
         
 <pre><code>>>> from probabilistic_serial_mechanism import probabilistic_serial_mechanism
>>> R = {0: [0,1,2], 1: [2,1,0]}
>>> probabilistic_serial_mechanism(R)
{0: array([ 1. ,  0.5,  0. ]), 1: array([ 0. ,  0.5,  1. ])}
>>> </code></pre>
