%-------------------------------------------------------------------------%
%
% File: scale_function(z1,z2,xc,h,h1,p)
%
% Goal: script that evaluates the scale function
%
% Inputs: z1:      first coordinate of the centre
%         z2:      second coordinate of the centre
%         h:       parameter for the scale function
%         h1:      parameter for the scale function
%         p:       if p=0 the code only computes the kernel matrix,
%                  otherwise evaluates also the differentiation matrices
%
% Outputs: val:      the value of the scale function
%          valx:     first order partial derivatives of val along the
%                    first dimension
%          valxx:    second order partial derivatives of val along the
%                    first dimension
%          valy:     first order partial derivatives of val along the
%                    second dimension
%          valyy:    second order partial derivatives of val along the
%                    second dimension
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
function [val,valx,valy,valxx,valyy] = scale_function(z1,z2,xc,h,h1,p)
if p == 1
    val =  sum(abs(1/pi*(atan(h.*(z1-xc(:,1)))).*...
        exp(-h1.*(z2-xc(:,2)).^2)));
    valx = sum(abs(val)./(val).*(1/pi*(h./(1+h.^2.*(z1-xc(:,1))...
        .^2)).*exp(-h1.*(z2-xc(:,2)).^2)));
    valy = (sum(abs(val)./(val).*(-2*h1.*(z2-xc(:,2))).*(1/(pi)*...
        atan(h.*(z1-xc(:,1)))).*exp(-h1.*(z2-xc(:,2)).^2)));
    valxx =  sum(abs(val)./(val).*((-2*h.^2.*(z2-xc(:,2)))./...
        ((1+h.^2.*(z1-xc(:,1)).^2)).^2).*(1./pi).*...
        exp(-h1.*(z2-xc(:,2)).^2));
    valyy = sum((1/pi*atan(h.*(z1-xc(:,1)))).*(-2.*h1.*exp(-h1.*...
        (z2-xc(:,2)).^2)+(2*h1.*(z2-xc(:,2))).^2.*...
        exp(-h1.*(z2-xc(:,2)).^2)));
else
    val =  sum(abs(1/pi*(atan(h.*(z1-xc(:,1)))).*...
        exp(-h1.*(z2-xc(:,2)).^2)));
end