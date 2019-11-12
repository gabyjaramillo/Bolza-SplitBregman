# Modified-Split-Bregman
We modify the Split Bregman algorithm to find minimizers of nonsmooth energy functionals in 1-d bounded domains.

# Instructions for use
There are three matlab files which are used to reproduced the examples presented in the paper by Gabriela Jaramillo and Shankar Venkataramani. 

## Minimizers
This is the main file that allows one to find minimizers for different energy functionals. One can choose between three non-convex potentials in the gradient variable representing a double well, a half-double well, and a tripple well potential. One can also choose between various examples of lower order potentials which are non-convex. 

## Obstacle
Calculates the convex envelop of non-convex potentials using the Split Bregman algorithm. This code is based on the paper by Tran and Osher: An L1 penalty method for general obstacle problems \url{https://epubs.siam.org/doi/abs/10.1137/140963303}

## Split_Bregman_combined
The algorithm in this file uses the Split Bregman to decompose the minimization into two subproblems linking the variable u and its gradient u_x via a constraint. The first subproblem is solved using Gauss-Seidel, while the second problem involving a nonsmooth functional is solved using a shrink operator (proximal gradient method).

