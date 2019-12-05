# Modified-Split-Bregman
We modify the Split Bregman algorithm to find minimizers of nonsmooth energy functionals in 1-d bounded domains.

# Instructions for use
There are two main folders. One where Dirichlet boundary conditions are used and that is labeled as Dirichlet_V3, and a second one where Natural boundary conditions are used and that is labeled Natural_V3. In each folder there are four matlab files which are used to reproduced the examples presented in the paper by Gabriela Jaramillo and Shankar Venkataramani. The other folders and files contained previous versions of the same code, and we suggest using version 3 as the most up to date version.

## Minimizers
This is the main file that allows one to find minimizers for different energy functionals. One can choose between three non-convex potentials in the gradient variable representing a double well, a half-double well, and a tripple well potential. One can also choose between various examples of lower order potentials which are non-convex. 

## Obstacle
Calculates the convex envelop of non-convex potentials using the Split Bregman algorithm. This code is based on the paper by Tran and Osher: An L1 penalty method for general obstacle problems \url{https://epubs.siam.org/doi/abs/10.1137/140963303}

## Split_Bregman_combined
The algorithm in this file uses the Split Bregman to decompose the minimization into two subproblems linking the variable u and its gradient u_x via a constraint. The first subproblem is solved using Gauss-Seidel, while the second problem involving a nonsmooth functional is solved using a shrink operator (proximal gradient method).

# Split_iteration
It adapts the file Minimizer.m so that one can iterate the calculations to compare the speed of the code for different values of the parameter $Delta x$
