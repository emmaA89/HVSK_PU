%-------------------------------------------------------------------------%
%
% File: gauss_PUM(ep,xc,xe)
%
% Goal: script that evaluates the differentiation matrices for the Gaussian
%       kernel
%
% Inputs: ep:         the shape parameter
%         xc:         the data sites
%         xe:         the evaluation points
%
% Outputs: w,wx,wy,wxx,wyy: evaluation of the gaussian function and of
%                           its derivatives
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
function [w,wx,wy,wxx,wyy] = gauss_PUM(ep,xc,xe)
dx = xe(1,1) - xc(:,1); dy = xe(1,2) - xc(:,2); % Initialize
r2 = dx.^2 + dy.^2; 
% Evaluate the RBF and its derivatives
we = exp(-ep^2*r2); w=we'; 
wex = (-2*ep^2*dx).*we; wx=wex'; wey = (-2*ep^2*dy).*we; wy=wey';
wexx =(-2*ep^2).*we+(4*dx.^2*ep^4).*we; wxx=wexx';
weyy =(-2*ep^2).*we+(4*dy.^2*ep^4).*we; wyy=weyy';
