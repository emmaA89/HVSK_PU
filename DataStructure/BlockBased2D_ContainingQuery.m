
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
% File: BlockBased2D_ContainingQuery(puctr,q,puradius,min_dsites)
%
% Goal: script that given a subdomain centre returns the index of
%       the square cell containing the subdomain centre
%
% Inputs: puctr:       subdomain centre
%         q:           number of square cells in one direction
%         puradius:    radius of PU subdomains
%         min_dsites:  minumim of data sites among two directions
%
% Outputs: index: the index of the square cell containing the subdomain
%                 centre
%
%-------------------------------------------------------------------------%
function [index] = BlockBased2D_ContainingQuery(puctr,q,puradius,...
    min_dsites)
k1 = 1;
while k1 <= q
    % Build the k1-th strip parallel to the y-axis
    x(k1) = min_dsites + puradius*k1;
    if  (puctr(1,1) <= x(k1))
        % Find the index of the strip parallel to the y-axis containig
        % the first coordinate of the subdomain centre
        idx_x = k1;
        break
    end
    k1 = k1 + 1;
end
k2 = 1; % Initialize
while k2 <= q
    % Build the k2-th strip parallel to the x-axis
    y(k2) = min_dsites + puradius*k2;
    if (puctr(1,2) <= y(k2))
        % Find the index of the strip parallel to the x-axis containig
        % the second coordinate of the subdomain centre
        idx_y = k2;
        break
    end
    k2 = k2 + 1;
end
% Define the index of the cell containing the subdomain centre
index = (idx_x-1)*q + idx_y;