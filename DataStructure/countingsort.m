
% The function is from Matlab Central Fileexchange (by Brian Moore)
% https://it.mathworks.com/matlabcentral/fileexchange/

%-------------------------------------------------------------------------%
%
% File: countingsort(x,r)
%
% Goal: script that sorts a vecor of integers
%
% Inputs: x:  a vector of length containing integers in the range
%             [1,...,r]
%         r:  an upper bound on max(x)
%
% Outputs: sx: the sorted (ascending) version of x
%
%-------------------------------------------------------------------------%
function [sx] = countingsort(x,r)
% Compute histogram
n = numel(x);
C = zeros(r,1);
for j = 1:n
    C(x(j)) = C(x(j)) + 1;
end
% Convert to cumulative values
for i = 2:r
    C(i) = C(i) + C(i - 1);
end
% Sort the array
sx = nan(n,1);
for j = n:-1:1
    sx(C(x(j))) = x(j);
    C(x(j)) = C(x(j)) - 1;
end