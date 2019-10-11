%-------------------------------------------------------------------------%
%
% File: CostFun(x,x0)
%
% Goal: script that computes greedy points by finding the maximum of
%       separatrion distances
%
% Inputs: x:     the greedy points, parameters to be optimized
%         x0:    the original point (the seed)
%
% Outputs: gw:  the value of the cost function
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
function [gw] = CostFun(x,x0)
a = DistanceMatrix(x0,x);
gw = -min(min(a)+eye(size(a)));