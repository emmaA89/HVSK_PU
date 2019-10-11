%-------------------------------------------------------------------------%
%
% File: greedy(neigh,x0,N)
%
% Goal: script that computes greedy points in a neighborhood
%
% Inputs: neigh:   the neighbors
%         x0:      the original point (the seed)
%         N:       number of greedy points
%
% Outputs: x:  the greedy points
%
% Calls on: 1. ConstrCostFun.m constrained cost function
%
% Authors:  S. De Marchi*, A. Martínez*, E. Perracchione*, M. Rossini^,
%           *Universita' di Padova,
%            Dipartimento di Matematica "Tullio Levi-Civita".
%           ^Universita' di Milano Bicocca, Dipartimento di Matematica
%            e Applicazioni
%
% Last modified: 21/11/17.
%
%-------------------------------------------------------------------------%
function [x] = greedy(neigh,x0,N)
% Set options for optimization
options=optimset('Display','off','MaxFunEval',300,'TolFun',1e-06,...
    'TolCon',1e-06,'TolX',1e-06);
for i=1:N
    xe = fminsearch(@(x) ConstrCostFun(x,x0,neigh),[sum(neigh(:,1)) ...
        sum(neigh(:,2))]./size(neigh,1),options);
    x0 = [x0;xe];
end
x=x0;