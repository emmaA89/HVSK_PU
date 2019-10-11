
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
% File: BlockBased2D_RangeSearch(puctr,puradius,dsites,index)
%
% Goal: find the data sites located in a given subdomain and the distances
%       between the subdomain centre and data sites
%
% Inputs: puctr:     subdomain centre
%         puradius:  radius of PU subdomains
%         dsites:    NX2 matrix representing a set of N data sites
%         index:     vector containing the indices of the data points
%                    located in the k-th block (the cell containing the
%                    subdomain centre) and in the neighbouring cells
%
% Outputs: idx:  vector containing the indices of the data points located
%                in a given PU subdomain
%          dist: vector containing the distances between the data sites
%                and the subdomain centre
%
%-------------------------------------------------------------------------%
function [idx, dist] = BlockBased2D_RangeSearch(puctr,puradius,dsites,...
    index)
N = size(dsites,1); dist1 = zeros(1,N); dist = []; idx = []; %Initialize
% Compute distances between the data sites and the centre
for i = 1:N
    dist1(i) = sqrt((puctr(1,1)-dsites(i,1))^2 + ...
        (puctr(1,2)-dsites(i,2))^2);
end
% Use a sort procedure to order distances
[sort_dist,IX] = sort(dist1);
N1 = size(sort_dist,2); j1 = 1; j2 = 1; %Initialize
% Find the data sites located in the given subdomain
while (j2 <= N1) && (sort_dist(j2) <= puradius)
    if nargin == 3
        idx(j1) = IX(j2); dist(j1) = dist1(IX(j2));
        j1 = j1 + 1; j2 = j2 + 1;
    else
        idx(j1) = index(IX(j2)); dist(j1) = dist1(IX(j2));
        j1 = j1 + 1; j2 = j2 + 1;
    end
end