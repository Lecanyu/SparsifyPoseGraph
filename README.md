# SparsifyPoseGraph
Sparsify pose graph and maintain the uncertainty. ***(Only support in Linux)***

This repository is based on several research papers 

1. [Nonlinear factor recovery for long-term SLAM](http://www2.informatik.uni-freiburg.de/~mazuran/papers/mazuran15ijrr.pdf). This paper proposed the NFR and re-implemented GLC method to sparsify pose graph. Both two method implementations have been included.

2. [Generic Node Removal for Factor-Graph SLAM](https://ri.cmu.edu/pub_files/2014/9/CarlevarisBianco14tro.pdf). This paper proposed the GLC method to sparsify pose graph. It addresses the issue that the direct composition between two nodes is not equivalent to marginalization (i.e. nodes delete).

3. [Estimating uncertain spatial relationships in robotics](https://www.researchgate.net/profile/Randall_Smith4/publication/221405213_Estimating_Uncertain_Spatial_Relationships_in_Robotics/links/0fcfd5141dc9b2bf2d000000.pdf). This paper gives an intuitive explanation for preliminary math.


You can check the test_xx.cpp in 'src/' about the algorithm details.


# Dependencies
1. [g2o](https://github.com/RainerKuemmerle/g2o)
2. [isam](http://people.csail.mit.edu/kaess/isam/). It can only be built on Linux.
3. several solver libraries

Refer to the CMakelists.txt for details
