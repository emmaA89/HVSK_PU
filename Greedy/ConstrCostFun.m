%-------------------------------------------------------------------------%
%
% File: ConstrCostFun(x,x0,neigh)
%
% Goal: script that defines constraints for the construction of greedy
%       points
%
% Inputs: x:     the greedy points, parameters to be optimized
%         x0:    the original point (the seed)
%         neigh: the neighbors
%
% Outputs: val:  the value of the cost function
%
% Calls on: 1. CostFun.m: the cost function
%
% Authors:  S. De Marchi*, A. Martínez*, E. Perracchione*, M. Rossini^,
%           *Universita' di Padova, Dipartimento di Matematica
%           "Tullio Levi-Civita".
%           ^Universita' di Milano Bicocca, Dipartimento di Matematica
%            e Applicazioni
%
% Last modified: 21/11/17.
%
%-------------------------------------------------------------------------%
function [val] = ConstrCostFun(x,x0,neigh)
% Impose constraints
if  x(1,1)>max(neigh(:,1)) || x(1,2)>max(neigh(:,2)) || ...
        x(1,1)<min(neigh(:,1)) || x(1,2)<min(neigh(:,2)) ...
        || x(1)<0 || x(2)<0 || x(1)>1 || x(2)>1
    val = Inf;
else
    % Evaluate the cost functions
    val = CostFun(x,x0);
end