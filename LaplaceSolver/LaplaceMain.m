%-------------------------------------------------------------------------%
%
% File: LaplaceMain(N,x,neval,ep,fun,Lfun,phi,sys_standard,...
%         iter,cpu_t)
%
% Goal: script that solves the Poisson equation by means of HVSKs
%
% Inputs: N:       number of central points
%         x:       the central points
%         neval:   number of evaluation points in one direction
%         ep :     shape parameter
%         fun:     the test function
%         Lfun:    Laplacian of the test function
%         phi:     the kernel function
%
% Notes: See Demo.m for the remaining inputs that identify the metod used
%        to solve the equation
%
% Outputs: the RMSE calculated at the evaluation points
%
% Calls on: matlab functions cointained in the folders:
%           1. DataStructure. Such routines organize points among the
%           patches
%           2. DiffMatrix. Such routines compute differentiation matrices
%           via standard bases
%           3. DiffMatrixVSK. Such routines compute differentiation
%           matrices via VSKs
%           3. Data. The data used for the examples
%
% Authors:  S. De Marchi*, A. Martínez*, E. Perracchione*, M. Rossini^,
%           *Universita' di Padova, Dipartimento di Matematica
%           "Tullio Levi-Civita".
%           ^Universita' di Milano Bicocca, Dipartimento di Matematica
%            e Applicazioni
%
% References: 1. [M. Bozzini, L. Lenarduzzi, M. Rossini, R. Schaback,
%                 Interpolation with variably scaled kernels,
%                 IMA J. Numer. Anal. 35 (2015), 199--219]
%             2. [A. Heryudono, E. Larsson, A. Ramage, L. von Sydow,
%                 Preconditioning for radial basis function partition of
%                 unity methods, J. Sci. Comput. 67 (2016), 1089--1109]
%             3. [V. Shcherbakov, E. Larsson, Radial basis function
%                 partition of unity methods for pricing vanilla basket
%                 options, Comput. Math. Appl 71 (2016), 185--200]
%
% Last modified: 21/11/17.
%
%-------------------------------------------------------------------------%
function [RMSE] = LaplaceMain(N,x,neval,ep,fun,Lfun,phi,sys_standard,...
    iter,cpu_t)
% Define the evaluation points
grid = linspace(0,1,neval); [xee,yee] = meshgrid(grid);
epoints = [xee(:) yee(:)]; Ne = size(epoints,1);
% Define boundary nodes
sn = sqrt(N); bdylin = linspace(0,1,sn)';
bdy0 = zeros(sn-1,1); bdy1 = ones(sn-1,1);
xb = [bdylin(1:end-1) bdy0; bdy1 bdylin(1:end-1);...
    flipud(bdylin(2:end)) bdy1; bdy0 flipud(bdylin(2:end))];
% Define the global set of nodes (central and boundary points)
dsites = [xb;x];
% The number of patches in one dimension
npu =  ceil(sqrt(N)/2);
% Define the PU centres
xp = linspace(0,1,npu); cellradius = sqrt(2)/npu; [xp,yp] = meshgrid(xp);
puctrs = [xp(:) yp(:)];
npu = size(puctrs,1); Nb = size(xb,1); N = size(dsites,1); % Initialize
% Parameter for q-by-q cell-based partitioning structure
q = ceil(1./cellradius);
% Define the data structures for the nodes and evaluation points
if q>2
    idx_ds = BlockBased2D_Structure(dsites,q,cellradius,0);
    idx_ep = BlockBased2D_Structure(epoints,q,cellradius,0);
    idx_bd = BlockBased2D_Structure(xb,q,cellradius,0);
end
for i=1:npu
    % Find the block cointaining the j-th subdomain centre
    index = BlockBased2D_ContainingQuery(puctrs(i,:),q,cellradius,...
        0);
    % Find data sites located in the j-th subdomain
    idx_ns = BlockBased2D_RangeSearch(puctrs(i,:),cellradius,...
        dsites(idx_ds{index},:),idx_ds{index});
    locpts(i).ind = sort(idx_ns);
    % Find evaluation points located in the j-th subdomain
    eidx_ns = BlockBased2D_RangeSearch(puctrs(i,:),cellradius,...
        epoints(idx_ep{index},:),idx_ep{index});
    elocpts(i).ind = sort(eidx_ns);
    Ni(i) = length(locpts(i).ind); % Initialize
end
% Define the pu weight
pu = W2weight(dsites,locpts,puctrs,cellradius);
epu = W2weight(epoints,elocpts,puctrs,cellradius);
B = spalloc(N,N,sum(Ni.^2)); % Initialize
tic
for i = 1:npu
    idx = locpts(i).ind; % Compute points on the j-th subdomain
    flag_vsk = 0; % Initialize
    % Compute the interpolation matrix
    A = DistanceMatrix(dsites(idx,:),dsites(idx,:));
    A1 = phi(ep,A);
    % Define a treshold and check if A is close to be singular
    [~,s] = svd(A1);
    if ep < 1
        if min(diag(s))<10^(-16)/ep^4
            flag_vsk = 1;
        end
    else
        if min(diag(s))<10^(-17)/ep^4 || ...
                min(diag(s))>10^(-11)
            flag_vsk = 1;
        end
    end
    if  flag_vsk == 1
        h = 7e-006; h1 = 5; % Initialize
        % Compute the differentiation matrices
        [A,Ax,Axx,Ay,Ayy] = RBF_PUM_diffmat_VSK_gauss(...
            dsites(idx,:),dsites(idx,:),h,h1,1);
        A_L = Axx+Ayy;
        % Construct the local local matrix
        D_L = diag(pu(i).w)*(A_L) + ...
            2*diag(pu(i).wx)*(Ax) + 2*diag(pu(i).wy)*(Ay) + ...
            diag(pu(i).wxx + pu(i).wyy)*(A);
        % Accumulate into the global collocation matrix
        B(idx,idx) = B(idx,idx) + D_L/A;
    else
        % Compute the interpolation matrix
        [A,Ax,Axx,Ay,Ayy] = RBF_PUM_diffmat_gauss(dsites(idx,:),...
            ep);
        A_L = Axx+Ayy;
        % Construct the local local matrix
        D_L = diag(pu(i).w)*(A_L) + ...
            2*diag(pu(i).wx)*(Ax) + 2*diag(pu(i).wy)*(Ay) + ...
            diag(pu(i).wxx + pu(i).wyy)*(A);
        % Accumulate into the global collocation matrix
        B(idx,idx) = B(idx,idx) + D_L/A;
    end
    % Compute the evaluation points on the j-th subdomain
    eidx = elocpts(i).ind;
    if (~isempty(eidx))
        if flag_vsk == 1
            E(i).Ae = RBF_PUM_diffmat_VSK_gauss(epoints(eidx,:),...
                dsites(idx,:),h,h1,0);
            E(i).Ae = diag(epu(i).w)*E(i).Ae/A;
        else
            E(i).Ae = phi(ep,DistanceMatrix(epoints(eidx,:),...
                dsites(idx,:)));
            E(i).Ae = diag(epu(i).w)*E(i).Ae/A;
        end
    else
        E(i).Ae = zeros(0,length(idx));
    end
end
t_matrix = toc;
if cpu_t == 1
    fprintf('Time for constructing  matrix    %e secs \n',t_matrix);
end
% Organize boundary conditions
B(1:Nb,1:N) = eye(Nb,N); f = fun(xb(:,1),xb(:,2));
g = Lfun(x(:,1),x(:,2)); rhs = [f;g];
% Solve the system
if sys_standard == 1 % Standard method
    tic
    u = B\rhs;
    t_sys = toc;
    if cpu_t == 1
        fprintf('Time for solving system          %e secs \n',t_sys);
    end
    % Compute exact solution and RMSE
    ue  = zeros(Ne,1);
    for i=1:npu
        idx=locpts(i).ind;
        eidx = elocpts(i).ind;
        ue(eidx) = ue(eidx) + E(i).Ae*u(idx);
    end
    exact = fun(epoints(:,1),epoints(:,2));
    RMSE = norm(exact-ue,2)/neval;
end
if iter == 1 % Use preconditioning
    P = amd(B); B = B(P,P); rhs=rhs(P); % Permute the matrix
    % Compute preconditioner - ILU
    setup.type = 'ilutp';
    setup.droptol = 1e-9;
    setup.udiag = 1;
    % Compute incomplete LU factorization
    tic
    [L,U] = ilu(B,setup);
    t_prec = toc;
    % Use GMRES to compute the solution
    tolinner = 1e-8; inneritmax = 20;
    nrest = 20;% Iterations before restarting
    Aprec = @(x) B*(U\(L\x));
   [y,~] = gmres(Aprec,rhs,nrest,tolinner,inneritmax);
    u = U\(L\y);
    t_sys = toc;
    it= t_prec + t_sys;
    if cpu_t == 1
        fprintf('Time to compute precond          %e secs \n',t_prec);
        fprintf('Time for solving system          %e secs \n',t_sys);
        fprintf('Total time                       %e secs \n',it);
    end
    u(P) = u; % Permute the solution
    % Compute exact solution and RMSE
    ue  = zeros(Ne,1);
    for i = 1:npu
        idx = locpts(i).ind;
        eidx = elocpts(i).ind;
        ue(eidx) = ue(eidx) + E(i).Ae*u(idx);
    end
    exact = fun(epoints(:,1),epoints(:,2));
    RMSE = norm(exact-ue,2)/neval;
end
