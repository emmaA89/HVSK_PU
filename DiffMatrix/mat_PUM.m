%-------------------------------------------------------------------------%
%
% File: mat_PUM(ep,xc,xe)
%
% Goal: script that evaluates the differentiation matrices for the Matern
%       C^3 kernel
%
% Inputs: ep:         the shape parameter
%         xc:         the data sites
%         xe:         the evaluation points
%
% Outputs: w,wx,wy,wxx,wyy: evaluation of the C^3 Matern function and of
%                           its the derivatives
%
% Notes: The original codes here modified are from
%        http://www.it.uu.se/research/scientific_computing/project/rbf
%        /software
%
% References: see e.g.
%          1. [A. Heryudono, E. Larsson, A. Ramage, L. von Sydow,
%          Preconditioning for radial basis function partition of unity
%          methods, J. Sci. Comput. 67 (2016), 1089--1109]
%          2. [V. Shcherbakov, E. Larsson, Radial basis function partition
%          of unity methods for pricing vanilla basket options, Comput.
%          Math. Appl 71 (2016), 185--200]
%
% Last modified: 21/11/17.
%
%-------------------------------------------------------------------------%
function [w,wx,wy,wxx,wyy] = mat_PUM(ep,xc,xe)
dx = xe(1,1) - xc(:,1); dy = xe(1,2) - xc(:,2); % Initialize
r = sqrt(dx.^2 + dy.^2); r2= dx.^2 + dy.^2; 
% Evaluate the RBF and its derivatives
we = exp(-ep*r).*(15+15*ep*r+6*ep^2*r2+ep^3*r.^3); w=we';
wex = -ep^2*dx.*exp(-ep*r).*(3+3*ep*r+ep^2*r2); wx=wex';
wey = -ep^2*dy.*exp(-ep*r).*(3+3*ep*r+ep^2*r2); wy=wey';
wexx =ep^2*exp(-ep*r).*(ep^3*dx.^2.*r-ep^2*dy.^2-3*ep*r-3); wxx=wexx';
weyy =ep^2*exp(-ep*r).*(ep^3*dy.^2.*r-ep^2*dx.^2-3*ep*r-3); wyy=weyy';