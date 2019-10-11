
% The function is from
% http://hdl.handle.net/2318/158790
%
% Remarks: Refer also to
%          [R. Cavoretto, A. De Rossi, E. Perracchione, Efficient
%          computation of partition of unity interpolants through a
%          block-based searching technique, Comput. Math. Appl. 71 (2016),
%          2568--2584]

%-------------------------------------------------------------------------%
%
% File: BlockBased2D_Structure(dsites,q,puradius,min_dsites)
%
% Goal: find the data sites located in each of the q^2 cells and in the
%       neighbouring cells
%
% Inputs: dsites:     NX2 matrix representing a set of N data sites
%         q:          number of cells in one direction
%         puradius:   radius of PU subdomains
%         min_dsites: minimum of data sites among two directions
%
% Calls on: countingsort: by B. Moore, from MATLAB Central File Exchange
%
% Outputs: idx_dsites: multiarray containing the indices of the data points
%                      located in k-th block and in the neighbouring blocks
%
%-------------------------------------------------------------------------%
function [idx_dsites] = BlockBased2D_Structure(dsites,q,puradius,...
    min_dsites)
N = size(dsites,1); k = 1; i1 = 1; % Initialize
% Sort with respect to the first coordinate the data sites
[dsites_sortx,IX] = sortrows(dsites,1);
for k1 = 1:q
    % Build the k1-th strip parallel to the y-axis
    x(k1) = min_dsites + puradius*k1; i2 = 1; idx_x = []; % Initialize
    while (i1 <= N) && (dsites_sortx(i1,1) <= x(k1))
        % Find the points located in the k1-th strip parallel to the y-axis
        idx_x(i2) = IX(i1); i1 = i1 + 1; i2 = i2 + 1;
    end
    idx_sortx = countingsort(idx_x,N); % Sort the indices
    % Sort with respect to the second coordinate data sites
    % located in the k1-th strip parallel to the y-axis
    [dsites_k1x,IX_kx] = sortrows(dsites(idx_sortx(:),:),2);
    j1 = 1; %Initialize
    for k2 = 1:q
        % Build the k2-th strip parallel to the x-axis
        y(k2) = min_dsites + puradius*k2; j2 = 1; idx = []; % Initialize
        while (j1 < i2) && (dsites_k1x(j1,2) <= y(k2))
            % Find the points located in the k2-th strip parallel to
            % the x-axis
            idx_y = IX_kx(j1); idx(j2) = idx_sortx(idx_y);
            j1 = j1 + 1; j2 = j2 + 1;
        end
        % Define the multiarray containing the indices of the data sites
        % located in the k-th cell
        idx_dsites_k{k} = idx'; k = k + 1;
    end
end
% Find data sites located in the k-th cell and in the neighbouring cells
for k = 1:q^2
    idx_dsites{k} = idx_dsites_k{k}; % Initialize
    % Indices of neighbouring cells
    neigh = [k-1,k+1,k+q,k-q,k+q+1,k-q+1,k-q-1,k+q-1];
    % Reduce the number of neighbouring cells for border cells
    border_cell = 0;
    if k == 1 % Corner cell
        border_cell = 1; % Reduce the number of neighbouring cells
        idx_dsites{k} = [idx_dsites_k{k};idx_dsites_k{k+1};...
            idx_dsites_k{k+q+1};idx_dsites_k{k+q}];
    end
    if k == q % Corner cell
        border_cell = 1; % Reduce the number of neighbouring cells
        idx_dsites{k} = [idx_dsites_k{k};idx_dsites_k{k-1};...
            idx_dsites_k{k+q};idx_dsites_k{k+q-1}];
    end
    if k == q^2 % Corner cell
        border_cell = 1; % Reduce the number of neighbouring cells
        idx_dsites{k} = [idx_dsites_k{k};idx_dsites_k{k-1};...
            idx_dsites_k{k-q};idx_dsites_k{k-q-1}];
    end
    if k == q^2-q+1 % Corner cell
        border_cell = 1; % Reduce the number of neighbouring cells
        idx_dsites{k} = [idx_dsites_k{k};idx_dsites_k{k+1};...
            idx_dsites_k{k-q};idx_dsites_k{k-q+1}];
    end
    for j = 2:q-1
        if (k == j) % Left border cells
            border_cell = 1; % Reduce the number of neighbouring cells
            for i = 1:size(neigh,2)
                if (neigh(i)>0)
                    idx_dsites{k} = [idx_dsites{k};idx_dsites_k{neigh(i)}];
                end
            end
        end
    end
    for j = 2:q-1
        if k == j*q % Upper border cells
            border_cell = 1; % Reduce the number of neighbouring cells
            for i = 1:size(neigh,2)
                if (neigh(i) ~= k+1) && (neigh(i) ~= k+q+1) ...
                        && (neigh(i) ~= k-q+1)
                    idx_dsites{k} = [idx_dsites{k};idx_dsites_k{neigh(i)}];
                end
            end
        end
    end
    for j = 1:q-2 % Right border cells
        if (k == q^2-j)
            for i = 1:size(neigh,2)
                border_cell = 1; % Reduce the number of neighbouring cells
                if  (neigh(i) <= q^2)
                    idx_dsites{k} = [idx_dsites{k};idx_dsites_k{neigh(i)}];
                end
            end
        end
    end
    for j = 1:q-2 % Lower border cells
        if (k == j*q+1)
            border_cell = 1; % Reduce the number of neighbouring cells
            for i = 1:size(neigh,2)
                if  (neigh(i) ~= k-1) && (neigh(i) ~= k+q-1) ...
                        && (neigh(i) ~= k-q-1)
                    idx_dsites{k} = [idx_dsites{k};idx_dsites_k{neigh(i)}];
                end
            end
        end
    end
    if (border_cell == 0)
        for i = 1:size(neigh,2)
            idx_dsites{k} = [idx_dsites{k};idx_dsites_k{neigh(i)}];
        end
    end
end