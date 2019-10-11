%-------------------------------------------------------------------------%
%
% Usage: run Demo2.m
%
% Goal: script that gives an example for the adaptive collocation via
%       RBF-PU that makes use of HVSK
%
% Calls on: LaplaceAdaptMain (script that performs adaptive partition of
%           unity collocation via HVSKs)
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
clear all; close all; clc; format short e; warning off
addpath('DataStructure'); addpath('DiffMatrix'); addpath('DiffMatrixVSK')
addpath('Data'); addpath('LaplaceSolver'); addpath('Greedy')
N = 100; % Number of central nodes
load('Data2'); % Define the central nodes
neval = 40; % Number of evaluation points
ep = 0.3; % The shape parameter
% The test function
fun = @(x,y) 1./20.*exp(4.*x).*cos(2.*x+y);
Lfun = @(x,y) 1./20.*(16.*exp(4.*x).*cos(2.*x+y)-...
    8.*exp(4.*x).*sin(2.*x+y)-8.*exp(4.*x).*sin(2.*x+y)-...
    4.*exp(4.*x).*cos(2.*x+y) -exp(4.*x).*cos(2.*x+y));
phi = @(ep,r)  exp(-ep*r).*(15+15*ep*r+6*ep.^2*r.^2+ep^3*r.^3); % The RBF
thetar = 1*10^(-5); thetac = 1*10^(-9); % The tolerances
[Max_points,iteration,RES,Num_points] = LaplaceAdaptMain(N,x,...
    neval,ep,fun,Lfun,phi,thetar,thetac);
fprintf('Final number of nodes                     %g \n',Num_points);
fprintf('Maximum number point per patch            %g \n',Max_points);
fprintf('Number of iteration                       %g \n',iteration);
fprintf('Residual at final iteration               %e \n',RES);
