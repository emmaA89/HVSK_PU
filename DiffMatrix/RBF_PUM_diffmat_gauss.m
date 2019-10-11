%-------------------------------------------------------------------------%
%
% File: RBF_PUM_diffmat_gauss(data,ep)
%
% Goal: script that computes the differentiation matrices for the Gaussian
%       kernel
%
% Inputs: data:       the collocation points
%         ep:         the shape parameter
%
% Calls on: gauss_PUM
%
% Outputs: A,Ax,Axx,Ay,Ayy: the differentiation matrices
%
% Calls on: gauss_PUM
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
function [A,Ax,Axx,Ay,Ayy] = RBF_PUM_diffmat_gauss(data,ep)
N = size(data,1); [A,Ax,Ay,Axx,Ayy] = deal(zeros(N,N)); % Initialize
% Generate differentiation matrices
for i = 1:N
    [psi,psix,psiy,psixx,psiyy] = gauss_PUM(ep,data,data(i,:));
    % Compute the differentiation matrices
    A(i,:) = psi; Ax(i,:) = psix; Ay(i,:) = psiy; Axx(i,:) = psixx;
    Ayy(i,:) = psiyy;
end
