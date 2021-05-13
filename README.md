# matroid_bandits

This is a implementation of algorithms in [Chen et al (2014)](https://www.microsoft.com/en-us/research/uploads/prod/2016/06/nips2014_cpe.pdf) and [Chen, Gupta and Li (2016)](http://proceedings.mlr.press/v49/chen16a.pdf). 

These algorithms return an optimal solution with high probability from bases of matroid in multi-armed-bandit setting.

CLUCB was proposed by Chen et al (2014), and Exact-ExpGap was proposed by Chen, Gupta and Li (2016).

I experiment with graphical matroid, so these algorithms return the maximum spanning tree of Graph $G$. I tried the following two cases.

- Give large weight to only one spanning tree
    -  CLUCB returns the optimal solution in one loop.
    - The number of samples in Exact-ExpGap increases due to $\Delta_e$ if the number of edges is increased.
- Give randamized weight to edges
    - Exact-ExpGap retuns the optimal solution quickly than CLUCB.
    - CLUCB can returns the optimal solution in one loop if the weight are distributed in small overlap, otherwise it doesn’t converge easily.
    