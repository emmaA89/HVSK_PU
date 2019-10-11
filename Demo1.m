%-------------------------------------------------------------------------%
%
% Usage: run Demo1.m
%
% Goal: script that gives an example for the computation of the
%       HVSK-PU collocation method
%
% Calls on: LaplaceMain (script that performs partition of
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
clear all; close all; clc; warning off; format short e
addpath('DataStructure'); addpath('DiffMatrix'); addpath('DiffMatrixVSK');
addpath('Data'); addpath('LaplaceSolver')
N = 81; % Number of central nodes
load('Data1'); % Define the central nodes
neval = 40; % Number of evaluation points
p = 20; % Number of test shape parameters
ep = logspace(-3,2,p); % Define the shape parameters
sys_standard = 1; % If sys_standard = 1 the code uses standard matlab
% backslash to solve the system
iter = 0; % If iter = 1 the code uses iterative methods (GMRES) to solve
% the system
% The test function and its laplacian
fun = @(x1,x2) sin(x1.^2+2*x2.^2)-sin(2*x1.^2+(x2-0.5).^2);
Lfun = @(x1,x2) 6*cos(x1.^2+2*x2.^2)-6*cos((x2-1/2).^2+2*x1.^2) ...
    -4*x1.^2.*sin(x1.^2+2*x2.^2)-16*x2.^2.*sin(x1.^2+2*x2.^2) ...
    +sin((x2-1/2).^2+2*x1.^2).*(2*x2-1).^2 ...
    + 16*x1.^2.*sin((x2-1/2).^2+2*x1.^2);
phi = @(ep,r) exp(-ep.^2.*r.^2); % The RBF
RMSE = []; % Initialize
% Below examples for HVSK and iterative methods (for direct ones set
% sys_standard = 1; and iter = 0;
figure
hold on
disp('------------------- Collocation via HVSK --------------------------')
for i=1:length(ep)
    cpu_t = 0;
    if i == 1 
        cpu_t = 1; % Display CPU times
    end
    RMSE(i) = LaplaceMain(N,x,neval,ep(i),fun,Lfun,phi,...
        sys_standard,iter,cpu_t);
end
[r1, r2] = min(RMSE);
fprintf('Optimal RMSE                     %e \n',RMSE(r2));
fprintf('At epsilon                       %e \n',ep(r2));
loglog(ep,RMSE,'-k')
axis square
grid on
