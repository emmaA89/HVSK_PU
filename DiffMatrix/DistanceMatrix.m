
% The function is from
% [G.E. Fasshauer, Meshfree Approximation Methods with Matlab,
% World Scientific, Singapore, 2007]

%-------------------------------------------------------------------------%
%
% File: DistanceMatrix(dsites,ctrs)
%
% Goal: computes the distance matrix between two sets of points
%
% Inputs: dsites:    a set of points
%         ctrs:      a set of RBF centres
%
% Outputs: DM:  the distance matrix
%
%-------------------------------------------------------------------------%
function DM = DistanceMatrix(dsites,ctrs)
[M,~] = size(dsites); [N,s] = size(ctrs);
DM = zeros(M,N);
% Accumulate sum of squares of coordinate differences
% The ndgrid command produces two MxN matrices:
%   dr, consisting of N identical columns (each containing
%       the d-th coordinate of the M data sites)
%   cc, consisting of M identical rows (each containing
%       the d-th coordinate of the N centers)
for d = 1:s
    [dr,cc] = ndgrid(dsites(:,d),ctrs(:,d));
    DM = DM + (dr-cc).^2;
end
DM = sqrt(DM);