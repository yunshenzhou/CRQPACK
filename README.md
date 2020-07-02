CRQ Matlab code
---------------

Purpose:

Solve the following constrained optimization problem
- min x'Ax
- st  x'x=1
- C'x=b

by various methods:

- Lanczos method [1]
- Projected power method [2]
- Explicit method [3]


Folders:
----------
- src:			Source code for CRQ solvers
- synthetic:		Driver routines to run the optimization problem for synthetic data
- imagecut:		Driver routines to run the optimization problem for image cut


Files in folder src:
----------
files for Lanczos method
- CRQ_Lanczos.m:	Solve the CRQ optimization problem by Lanczos method
- QEPmin.m		Solve QEPmin problem
- LGopt.m		Solve LGopt problem
- rLGopt.m		Solve rLGopt by solving secular equation

files for projected power method
- CRQ_ppm.m		Solve the CRQ optimization problem by projected power method [1]

files for explicit method
- CRQ_explicit.m	Solve the CRQ optimization problem by explicit method [2]
- f.m			Secular function
- fp.m			Derivative of secular function

Files in folder synthetic:
----------
- correct.m		Test the correctness and the history of the error for CRQopt
- CRQsharp.m		Examples where the bound is sharp
- CRQnotsharp.m		Examples where the bound is not sharp
- QEPres.m		Test the residual of the QEP

files in folder auxiliary:
- ErrorHistory.m	get the error history of the objective function and solution 				vector from the structure info
- UpperBound.m		get the upper bound for the error of objective function and 				solution vector in each iterations

Folders in folder imagecut:
----------
- data			Files for images and labels
- examples		Driver routine for image cut examples
- auxiliary		Code for transforming the image and label files to matrices in 				optimization problem

files in folder examples:
- demo_exist.m		Run image cut problems with labels loaded from files
- demo_new_pts.m	Run image cut problems with user defined labels
To run demo_new_pts.m, first click the labels of the background and press Enter, then click the labels of the object and press Enter, then the image cut result is shown.

files in folder auxiliary:
- a_times_b_cmplx.cpp	provide Matrix vector multiplications and solve triangular systems
- affinityic.cpp	determine the affinity of pixel pairs
- cimgnbmap.cpp		compute the neighborhood index matrix of an image with each 				neighborhood sampled
- computeEdges.m	compute the edge in imageX
- computeW.m		compute the weight matrix W
- createB.m		create the matrices for linear constraints
- doog1.m and doog2.m	make difference of offset gaussians kernel
- fft_filt_2.m		fft-based filtering
- gaussian.m		Evaluate the multi-variate density with mean vector m and 				covariance matrix C for the input vector x
- ICgraph.m		get similarity matrix based on intervening contours
- imread_ncut.m		read the image and label files
- make_filterbank_even2.m and make_filterbank_odd2.m 
			make filter banks
- mex_w_times_x_symmetric.cpp 
			multiply a sparse symmetric matrix with a vector
- ncut.m		compute normalized cut
- NcutImage.m		get the discretized normalized cut
- quadedgep.m		get locations and gradients of an ordered list of edges
- showmask.m		plot the result of image cut
- sparsifyc.cpp		sparsify the input matrix
- spmtimesd.cpp		compute a sparse matrix times a diagonal matrix

Note: 
1. files in folder auxiliary in image cut is adopted from code attached in [3].
2. To run image cut examples, you need to select imagecut/auxiliary as the current folder in MATLAB run the following code to compile cpp code:

mex affinityic.cpp
mex cimgnbmap.cpp
mex mex_w_times_x_symmetric.cpp
mex sparsifyc.cpp
mex spmtimesd.cpp


References:
----------
- [1] Zhou, Yunshen, Zhaojun Bai, and Ren-Cang Li. "Linear Constrained Rayleigh Quotient Optimization: Theory and Algorithms." arXiv preprint arXiv:1911.02770 (2019).
- [2] Xu, Linli, Wenye Li, and Dale Schuurmans. "Fast normalized cut with linear constraints." Computer Vision and Pattern Recognition, 2009. CVPR 2009. IEEE Conference on. IEEE, 2009.
- [3] Gander, Walter, Gene H. Golub, and Urs von Matt. "A constrained eigenvalue problem." Linear Algebra and its applications 114 (1989): 815-839.
- [4] Shi, Jianbo, and Jitendra Malik. "Normalized cuts and image segmentation." Departmental Papers (CIS) (2000): 107.
 
