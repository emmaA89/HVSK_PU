%-------------------------------------------------------------------------%
%
% File: LaplaceAdaptMain(N,x,neval,ep,fun,Lfun,phi,thetar,thetac)
%
% Goal: script that solves the Poisson equation by means of adaptive
%       HVSK-PU
%
% Inputs: N:       number of central points
%         x:       the central points
%         neval:   number of evaluation points in one direction
%         ep :     shape parameter
%         fun:     the test function
%         Lfun:    Laplacian of the test function
%         phi:     the kernel function
%         thetar:  lower bound for the residual tolerance
%         thetac:  upper bound for the residual tolerance
%
% Outputs: Max_points:  maximum number of points per patch
%          Iteration:   the number of steps needed to reach the accuracy
%          RES:         the final residual
%          Num_points:  the final number of nodes
%
% Calls on: matlab functions cointained in the folders:
%           1. DataStructure. Such routines organize points among the
%           patches
%           2. DiffMatrix. Such routines compute differantiation matrices
%           via standard bases
%           3. DiffMatrixVSK. Such routines compute differantiation
%           matrices via VSKs
%           3. Data. The data used for the examples
%
% Authors:  S. De Marchi*, A. Martínez*, E. Perracchione*, M. Rossini^,
%           *Universita' di Padova,
%            Dipartimento di Matematica "Tullio Levi-Civita".
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
function [Max_points,iteration,RES,Num_points]=LaplaceAdaptMain(N,...
    x,neval,ep,fun,Lfun,phi,thetar,thetac)
N_1 = N; % Initialize
% Define the evaluation points
grid = linspace(0,1,neval); [xee,yee] = meshgrid(grid);
xe = [xee(:) yee(:)]; Ne = size(xe,1);
% The number of patches in one dimension
Np =  2*ceil(sqrt(N)/2); % The number of patches in one dimension
xp = linspace(0,1,Np);
cellradius = sqrt(2)/Np; % Define the patch radius and PU centres
[xp,yp] = meshgrid(xp); Cp = [xp(:) yp(:)]; Np = size(Cp,1);
% Initialize
k = 1; ref = true; y = []; max_n =[];
while any(ref)
    % Compute boundary conditions
    N_2 = size(x,1); xb = [];
    sn = floor(sqrt(N_2)); bdylin = linspace(0,1,sn)';
    bdy0 = zeros(sn-1,1); bdy1 = ones(sn-1,1);
    xb = [bdylin(1:end-1) bdy0; bdy1 bdylin(1:end-1);...
        flipud(bdylin(2:end)) bdy1; bdy0 flipud(bdylin(2:end))];
    xc = [xb;x]; N = size(xc,1);
    if k==1 % Compute the test nodes
        ya = haltonseq(2*N_2,2); y = ya(N_2+1:end,:);
    else
        if size(y1,1) > N_1/10 % Compute the test nodes
            ya = haltonseq(2*N_2,2); y = ya(N_2+1:end,:);
        else
            % If few nodes must be added then use greedy nodes
            for j1 = 1:size(y1,1)
                ya = []; [neighbors,~,~] = kNearestNeighbors([xc; y],...
                    y1(j1,:),2*k+1);
                neigh = [y1(j1,:); neighbors(1:end,:)];
                seed = neigh; % Define the seed for greedy points
                [ya] = greedy(neigh,[xc;y],1);
                y = [y;ya(end,:)];
            end
        end
    end
    % Eventually remove nodes if they are on the boundary (this is due
    % to greedy algorithm)
    i4 = 1; id1 = []; % Initialize
    for i2 = 1:size(y,1)
        for i3 = 1:size(xc,1)
            if abs(y(i2,1)-xc(i3,1))<10^(-12) && ...
                    abs(y(i2,2)-xc(i3,2))<10^(-12)
                id1(i4) = i2;
                i4 = i4+1;
            end
        end
    end
    y(id1,:) = [];
    % Initialize
    Nb = size(xb,1); q = ceil(1./cellradius);
    % Define the data structures for the nodes and evaluation points
    if q>2
        idx_ds = BlockBased2D_Structure(xc,q,cellradius,0);
        idx_ep = BlockBased2D_Structure(xe,q,cellradius,0);
        idx_bd = BlockBased2D_Structure(xb,q,cellradius,0);
        idx_y = BlockBased2D_Structure(y,q,cellradius,0);
    end
    for i=1:Np
        index = BlockBased2D_ContainingQuery(Cp(i,:),q,cellradius,...
            0);
        % Find data sites located on the j-th subdomain
        idx = BlockBased2D_RangeSearch(Cp(i,:),cellradius,...
            xc(idx_ds{index},:),idx_ds{index});
        locpts(i).ind=sort(idx);
        eidx = BlockBased2D_RangeSearch(Cp(i,:),cellradius,...
            xe(idx_ep{index},:),idx_ep{index});
        elocpts(i).ind=sort(eidx);
        nidx = BlockBased2D_RangeSearch(Cp(i,:),cellradius,...
            xb(idx_bd{index},:),idx_bd{index});
        nlocpts(i).ind=sort(nidx);
        yidx = BlockBased2D_RangeSearch(Cp(i,:),cellradius,...
            y(idx_y{index},:),idx_y{index});
        ylocpts(i).ind=sort(yidx);
        Ni(i) = length(locpts(i).ind);
    end
    max_n(k) = max(Ni(:));
    % Define the pu weight
    pu = W2weight(xc,locpts,Cp,cellradius);
    epu = W2weight(xe,elocpts,Cp,cellradius);
    epuy = W2weight(y,ylocpts,Cp,cellradius);
    % Initialize
    B = [];
    B = spalloc(N,N,sum(Ni.^2));
    % Compute the Laplacian
    for i = 1:Np
        ind = locpts(i).ind;
        ind_y = ylocpts(i).ind;
        ind_e = elocpts(i).ind;
        flag_vsk = 0; % Initialize
        % Compute the interpolation matrix
        A1 = DistanceMatrix(xc(ind,:),xc(ind,:));
        A1 = phi(ep,A1);
        % Define a treshold and check if A is singular
        if rank(A1,10^(-16)/ep^4) < length(ind)
            flag_vsk = 1;
        end
        if  flag_vsk == 1
            h1 = 5; h = 7*1e-006;
            [A,Ax,Axx,Ay,Ayy] = RBF_PUM_diffmat_VSK_mat(xc(ind,:),...
                xc(ind,:),h,h1,1);
            A_L = Axx+Ayy;
            % Construct the local local matrix
            D_L = diag(pu(i).w)*(A_L) + ...
                2*diag(pu(i).wx)*(Ax) + 2*diag(pu(i).wy)*(Ay) + ...
                diag(pu(i).wxx + pu(i).wyy)*(A);
            % Accumulate into the global collocation matrix
            B(ind,ind) = B(ind,ind) + D_L/(A);
        else
            [A,Ax,Axx,Ay,Ayy] = RBF_PUM_diffmat_mat(xc(ind,:),ep);
            A_L = Axx+Ayy;
            % Construct the local local matrix
            D_L = diag(pu(i).w)*(A_L) + ...
                2*diag(pu(i).wx)*(Ax) + 2*diag(pu(i).wy)*(Ay) + ...
                diag(pu(i).wxx + pu(i).wyy)*(A);
            % Accumulate into the global collocation matrix
            B(ind,ind) = B(ind,ind) + D_L/(A);
        end
        % Build the evaluation matrices
        if (~isempty(ind_e))
            if  flag_vsk == 1
                E(i).Ae = RBF_PUM_diffmat_VSK_mat(xe(ind_e,:),xc(ind,:),h,...
                    h1,0);
                E(i).Ae = diag(epu(i).w)*E(i).Ae/A;
            else
                E(i).Ae  = DistanceMatrix(xe(ind_e,:),xc(ind,:));
                E(i).Ae = phi(ep,E(i).Ae);
                E(i).Ae = diag(epu(i).w)*E(i).Ae/A;
            end
        else
            E(i).Ae = zeros(0,length(ind));
        end
        if (~isempty(ind_y))
            if  flag_vsk == 1
                E1(i).Ae  = RBF_PUM_diffmat_VSK_mat(y(ind_y,:),xc(ind,:),...
                    h,h1,0);
                E1(i).Ae = diag(epuy(i).w)*E1(i).Ae/A;
            else
                E1(i).Ae  = DistanceMatrix(y(ind_y,:),xc(ind,:));
                E1(i).Ae = phi(ep,E1(i).Ae);
                E1(i).Ae = diag(epuy(i).w)*E1(i).Ae/A;
            end
        else
            E1(i).Ae = zeros(0,length(ind));
        end
    end
    % Organize boundary conditions
    B(1:Nb,1:N) = eye(Nb,N);
    f = fun(xb(:,1),xb(:,2));
    g = Lfun(x(:,1),x(:,2));
    rhs = [f;g];
    u = B\rhs; % Solve the system
    ue  = zeros(Ne,1);
    for i=1:Np
        ind = locpts(i).ind;
        ind_e = elocpts(i).ind;
        ue(ind_e) = ue(ind_e) + E(i).Ae*u(ind);
    end
    uey  = zeros(size(y,1),1);
    for i=1:Np
        ind = locpts(i).ind;
        ind_y = ylocpts(i).ind;
        if size(ind_y,1)>0
            uey(ind_y) = uey(ind_y) + E1(i).Ae*u(ind);
        end
    end
    resid=[]; % Compute residuals and remove unecessary nodes
    resid = abs(uey -fun(y(:,1),y(:,2)));
    coarsen = []; coarsen = find(resid<thetac);
    y_coar = y(coarsen,:);
    x_coar=[]; i2=1;
    for i1=1:size(y_coar,1)
        [~,neighborIds,~] = kNearestNeighbors(x,y_coar(i1,:),1);
        x_coar(i2)=neighborIds;
        i2=i2+1;
    end
    y_tcoar = [x(x_coar,:)];
    x(x_coar,:) = [];
    % Find nodes that need to be added
    ref = find(resid > thetar);
    y1 = []; y1 = [y(ref,:);y_tcoar];
    % The new set of nodes
    x = [x;y1];
    exact = fun(xe(:,1),xe(:,2));
    residual(k)= max(resid(:));
    points(k)=size([x; xb],1);
    k=k+1;
end
RES = residual(end);
Max_points = max_n(end);
iteration = k-1;
Num_points = points(end);
