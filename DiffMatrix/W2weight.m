%-------------------------------------------------------------------------%
%
% File: W2weight(xc,locpts,Cp,cellradius)
%
% Goal: script that computes the PU weights
%
% Inputs: xc:         the collocation points
%         locpts:     the collocation points on the subdomains
%         Cp:         centres of the PU subdomains
%         cellradius: the radius of patches
%
% Outputs: pu: the PU weights
%
% Calls on: PU_weight
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
function [pu] = W2weight(xc,locpts,Cp,cellradius)
Np = size(Cp,1); % Initialize
% Construct the weight functions
[phi,phix,phiy,phixx,phiyy] = PU_weight(xc,locpts,Cp,cellradius);
% Compute the sums of the generating functions and their derivatives
s = sum(phi,2); sx = sum(phix,2); sy = sum(phiy,2); sxx = sum(phixx,2);
syy = sum(phiyy,2);
% Compute the weights for PU collocation
for i = 1:Np
    loc = locpts(i).ind;
    if isempty(loc)
        loc = find([0;0;0]); s(loc)=find([0;0;0]);
    end
    pu(i).w = phi(loc,i)./s(loc);
    pu(i).wx = phix(loc,i)./s(loc) - phi(loc,i).*sx(loc)./s(loc).^2;
    pu(i).wy = phiy(loc,i)./s(loc) - phi(loc,i).*sy(loc)./s(loc).^2;
    pu(i).wxx = -2*phix(loc,i).*sx(loc)./s(loc).^2 + ...
        phixx(loc,i)./s(loc) + phi(loc,i).*(2*sx(loc).^2./s(loc).^3 - ...
        sxx(loc)./s(loc).^2);
    pu(i).wyy = -2*phiy(loc,i).*sy(loc)./s(loc).^2 + ...
        phiyy(loc,i)./s(loc) + phi(loc,i).*(2*sy(loc).^2./s(loc).^3 - ...
        syy(loc)./s(loc).^2);
end