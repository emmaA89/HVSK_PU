%-------------------------------------------------------------------------%
%
% File: RBF_PUM_diffmat_VSK_gauss(epoints,dsites,h,h1,p)
%
% Goal: script that constructs VSKs differentiation matrices for the
%       Gaussian kernel
%
% Inputs: epoints: the evaluation points
%         dsites:  the data sites
%         h:       parameter for the scale function
%         h1:      parameter for the scale function
%         p:       if p=0 the code only computes the kernel matrix,
%                  otherwise evaluates also the differentiation matrices
%
% Outputs: A:      kernel matrix
%          Ax:     first order partial derivatives of A along the first
%                  dimension
%          Axx:    second order partial derivatives of A along the first
%                  dimension
%          Ay:     first order partial derivatives of A along the second
%                  dimension
%          Ayy:    second order partial derivatives of A along the second
%                  dimension
%
% Calls on: gauss_PUM_VSK: evaluates the differentiation matrices with the
%                          scale function
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
function [A,Ax,Axx,Ay,Ayy] = RBF_PUM_diffmat_VSK_gauss(epoints,dsites,...
    h,h1,p)
N1 = length(epoints); N = length(dsites); % Initialize
[A,Ax,Ay,Axx,Ayy] = deal(zeros(N1,N));
% Generate differentiation matrices
if p == 1
    for i = 1:N1
        [psi,psix,psiy,psixx,psiyy] = gauss_PUM_VSK(dsites,epoints(i,:),...
            h,h1,p);
        A(i,:) = psi; Ax(i,:) = psix; Ay(i,:) = psiy; Axx(i,:) = psixx;
        Ayy(i,:) = psiyy;
    end
end
% If p = 0 only compute the distance matrix
if p == 0
    for i = 1:N1
        [psi] = gauss_PUM_VSK(dsites,epoints(i,:),h,h1,p);
        A(i,:) = psi;
    end
    Ax = 0; Ay = 0; Axx = 0; Ayy = 0;
end