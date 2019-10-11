This software computes the solution of Poisson problems via Hybrid
Variably Scaled Kernels (HVSKs). 

Authors:  S. De Marchi*, A. Martínez*, E. Perracchione*, M. Rossini^
                 *Universita' di Padova, Dipartimento di Matematica 
                 "Tullio Levi-Civita".
                 ^Universita' di Milano Bicocca, Dipartimento di Matematica
                 e Applicazioni

Contacts:  demarchi@math.unipd.it; acalomar@math.unipd.it;
                 emma.perracchione@math.unipd.it; milvia.rossini@unimib.it;

Last modified: 21/11/17.

The script Demo1.m provides an example for the solution of Poisson 
problems via RBF-PU collocation methods. The script Demo2.m is an example 
for the adaptive method. They make use of HVSKs to enanche the stability.

In what follows we give a brief description of the implemented Matlab 
functions. For more details concerning input and output parameters and the 
usage of single functions, see the comments within the Matlab routines.

The folder LaplaceSolver contains:
1. LaplaceMain.m: it is the main function that computes the RBF-PU 
    collocation via HVSKs
2. LaplaceAdaptMain.m: it is the main function that computes the adaptive 
    RBF-PU collocation via HVSKs

The folder DiffMatrixVSK contains:
3. RBF_PUM_diffmat_VSK_gauss.m: it constructs VSKs differentiation matrices
    for the gaussian 
4. gauss_PUM_VSK.m: it evaluates VSKs differentiation matrices for the 
    gaussian
5. RBF_PUM_diffmat_VSK_mat.m: it constructs VSKs differentiation matrices
    for the Matern C^3
6. mat_PUM_VSK.m: it evaluates VSKs differentiation matrices for the 
    Matern C^3
7. Scale_function.m: it defines the scale funtion used for VSKs

The folder DiffMatrix contains:
8. DistanceMatrix.m: it computes the distance matrix of two sets of points
9. RBF_PUM_diffmat_gauss.m: it constructs differentiation matrices via 
    standard bases for the gaussian
10. gauss_PUM.m: it evaluates the differentiation matrices via standard 
     bases for the gaussian
11. RBF_PUM_diffmat_mat.m: it constructs differentiation matrices via 
     standard bases for the matern C^3
12. mat_PUM.m: it evaluates the differentiation matrices via standard 
     bases for the matern C^3
13. W2weight.m: it construct the PU weight functions 
14. PU_weight.m: it evaluates the PU weight functions 

The folder Data_Structure contains:
15. BlockBased2D_ContainingQuery.m: it finds the corresponding block for 
     each patch
16. BlockBased2D_Structure.m: it constructs the structure to organize  
     points among patches
17. BlockBased2D_RangeSearch.m: it finds points lying on each patch
18. countingsort.m: it is an optimized searching procedure for integers
19. kNearestNeighbors.m: it finds the k nearest neighbors set

The folder Greedy:
20. greedy.m: it generates greedy points
21. ConstrCostFun.m: it is a constrained function for greedy points
22. CostFun.m: it defines the cost function for greedy points

The folder Data:
20. Data1.mat: it is the dataset used for Demo1
21. Data2.mat: it is the dataset used for Demo2
22. haltonseq.m: it generates Halton data

Remarks:
a. the function DistanceMatrix.m is from the book Meshfree Approximations  
    Methods with Matlab, Fasshauer, World Scientific, 2007.
b. routines 15--17 are from  [R. Cavoretto, A. De Rossi, E. Perracchione,
   Efficient computation of partition of unity interpolants through a 
   block-based searching technique,  Comput. Math. Appl. 71 (2016), 
   2568--2584] and can be downloaded at http://hdl.handle.net/2318/158790.
c. countingsort.m, kNearestNeighbors.m and haltonseq.m can be downloaded at 
   https://it.mathworks.com/matlabcentral/fileexchange/.
e. Routines 6--9 have been modified and are originally from 
    http://www.it.uu.se/research/scientific_computing/project/rbf...
    /software/rbfpum_convdiff
f. For the QR-RBF codes refer to 
    http://www.it.uu.se/research/scientific_computing/software/rbf_qr
