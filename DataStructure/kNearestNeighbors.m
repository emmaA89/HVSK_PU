
% The function is from Matlab Central Fileexchange (by Ani)
% https://it.mathworks.com/matlabcentral/fileexchange/15562-k-nearest...
% -neighbors

%-------------------------------------------------------------------------%
%
% File: kNearestNeighbors(dataMatrix,queryMatrix,k)
%
% Goal: script that computes the k Nearest Neighbors set
%
% Inputs: dataMatrix:  NxD-N vectors with dimensionality
%                      D (within which we search for the nearest
%                      neighbors)
%         queryMatrix: MxD-M query vectors with
%                      dimensionality D
%        k: Number of nearest neighbors desired
%
% Output: neighbors:         k nearest neighbors
%         neighborIds:       indices of k nearest neighbors
%         neighborDistances: distance between neighbors
%
%-------------------------------------------------------------------------%
function [neighbors, neighborIds, neighborDistances] = ...
    kNearestNeighbors(dataMatrix,queryMatrix,k)
neighborIds = zeros(size(queryMatrix,1),k);
neighborDistances = neighborIds;
numDataVectors = size(dataMatrix,1);
numQueryVectors = size(queryMatrix,1);
for i = 1:numQueryVectors
    dist = sum((repmat(queryMatrix(i,:),numDataVectors,1)- ...
        dataMatrix).^2,2);
    [sortval, sortpos] = sort(dist,'ascend');
    neighborIds(i,:) = sortpos(1:k);
    neighborDistances(i,:) = sqrt(sortval(1:k));
end
for i = 1:k
    neighbors(i,:) = dataMatrix(neighborIds(1,i),:);
end