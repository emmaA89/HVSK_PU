%-------------------------------------------------------------------------%
%
% File: PU_weight(xc,locpts,Cp,cellradius)
%
% Goal: script that evaluates Wendland's function and its derivatives
%
% Inputs: xc:         the collocation points
%         locpts:     the collocation points on the subdomains
%         Cp:         centres of PU subdomains
%         cellradius: the radius of patches
%
% Outputs: phi,phix,phiy,phixx,phiyy: evaluation of the PU weights and
%                                     the derivatives
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
function [phi,phix,phiy,phixx,phiyy] = PU_weight(xc,locpts,Cp,cellradius)
Np = size(Cp,1);
for i = 1:length(Np)
    N(i) = length(locpts(i).ind);
end
zs = spalloc(size(xc,1),Np,sum(N));
[phi,phix,phiy,phixx,phiyy] = deal(zs);
for i = 1:Np
    rho = 1/cellradius;
    dx = xc(locpts(i).ind,1)-Cp(i,1);
    dy = xc(locpts(i).ind,2)-Cp(i,2);
    r = sqrt(dx.^2 + dy.^2);
    r1 = r.*rho;
    phi(locpts(i).ind,i) = (4*rho*r+1).*((1-rho*r)).^4;
    phix(locpts(i).ind,i) = -20*rho^2*dx.*((1-rho*r)).^3;
    phiy(locpts(i).ind,i) = -20*rho^2*dy.*((1-rho*r)).^3;
    if r == 0
        phixx=20/cellradius.^2; phiyy=20/cellradius.^2;
    else
        phixx(locpts(i).ind,i)=20*((1-r1)).^2.*(-r1+4*...
            (dx.^2./cellradius.^2)+(dy.^2./cellradius.^2))./...
            cellradius.^2./r1;
        phiyy(locpts(i).ind,i)=20*((1-r1)).^2.*(-r1+(dx.^2....
            /cellradius.^2)+4*(dy.^2./cellradius.^2))./cellradius.^2./r1;
    end
end