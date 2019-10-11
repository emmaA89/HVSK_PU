%-------------------------------------------------------------------------%
%
% File: gauss_PUM_VSK(dsites,epoints,h,h1,p)
%
% Goal: script that evaluates VSKs differentiation matrices for the
%       Gaussian kernel
%
% Inputs: epoints: the evaluation points
%         dsites:  the data sites
%         h:       parameter for the scale function
%         h1:      parameter for the scale function
%         p:       if p=0 the code only computes the kernel matrix,
%                  otherwise evaluates also the differentiation matrices
%
% Outputs: w:      kernel function
%          wx:     first order partial derivatives of w along the first
%                  dimension
%          wxx:    second order partial derivatives of w along the first
%                  dimension
%          wy:     first order partial derivatives of w along the second
%                  dimension
%          wyy:    second order partial derivatives of w along the second
%                  dimension
%
% Calls on: scale_function: evaluates the scale function and its
%                           derivatives
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
function [w,wx,wy,wxx,wyy] = gauss_PUM_VSK(dsites,epoints,h,h1,p)
% If h is not a scalar (faster if h is a scalar)
if length(h) > 1
    % If p=0 only compute the distance matrix
    if p == 0
        for i = 1:size(dsites,1)
            % Evaluate the scale function
            [val_ep] = scale_function(epoints(1,1),epoints(1,2),...
                dsites,h(i),h1,0);
            [val_ds] = scale_function(dsites(1,1),dsites(1,2),dsites,...
                h(i),h1,0);
            a1 = val_ep-val_ds; dx=epoints(1,1) - dsites(i,1);
            dy = epoints(1,2) - dsites(i,2);
            r22(i) = (dx.^2+dy.^2+(a1).^2);
        end
        % Compute the VSK matrix
        we = exp(-r22);  w = we';
        wx = 0; wy = 0; wxx = 0; wyy = 0;
    end
    if p == 1
        for i = 1:size(dsites,1)
            % Evaluate the scale function and its derivatives
            [val_ep,valx_ep,valy_ep,valxx_ep,valyy_ep] = scale_function...
                (epoints(1,1),epoints(1,2),dsites,h(i),h1,1);
            [val_ds] = scale_function(dsites(1,1),dsites(1,2),...
                dsites,h(i),h1,0);
            a1 = val_ep-val_ds;
            dx = epoints(1,1) - dsites(i,1);
            dy = epoints(1,2) - dsites(i,2);
            r(i) = sqrt(dx.^2+dy.^2+(a1).^2);
            rx(i) = dx+a1*valx_ep;
            ry(i) = dy+a1*valy_ep;
            rxx(i) = 1+valx_ep.^2+(a1).*valxx_ep;
            ryy(i) = 1+valy_ep.^2+(a1).*valyy_ep;
        end
        % Compute the VSK differentiation matrices
        we = exp(-r.^2); w=we';
        wex = -2*rx.*we; wx=wex';
        wey = -2*ry.*we; wy=wey';
        wexx = 4.*we.*rx.^2-2.*we.*rxx; wxx=wexx';
        weyy =4.*we.*ry.^2-2.*we.*ryy; wyy=weyy';
    end
else
    if p == 0
        % Evaluate the scale function
        [val_ep] = scale_function(epoints(1,1),epoints(1,2),dsites,h,h1,0);
        [val_ds] = scale_function(dsites(1,1),dsites(1,2),dsites,h,h1,0);
        a1 = val_ep-val_ds;
        dx = epoints(1,1) - dsites(:,1);
        dy = epoints(1,2) - dsites(:,2);
        r22 = ((dx).^2+(dy).^2+(a1).^2);
        % Compute the VSK matrix
        we = exp(-r22);  w = we';
        wx = 0; wy = 0; wxx = 0; wyy = 0;
    end
    if p == 1
        % Evaluate the scale function and its derivatives
        [val_ep,valx_ep,valy_ep,valxx_ep,valyy_ep] = scale_function...
            (epoints(1,1),epoints(1,2),dsites,h,h1,1);
        [val_ds] = scale_function(dsites(1,1),dsites(1,2),dsites,h,h1,0);
        a1 = val_ep-val_ds;
        dx = epoints(1,1) - dsites(:,1);
        dy = epoints(1,2) - dsites(:,2);
        r22 = (dx.^2+dy.^2+(a1).^2);
        rx = dx+a1*valx_ep;
        ry = dy+a1*valy_ep;
        rxx = 1+valx_ep.^2+(a1).*valxx_ep;
        ryy = 1+valy_ep.^2+(a1).*valyy_ep;
        % Compute the VSK differentiation matrices
        we = exp(-r22); w=we';
        wex = -2*rx.*we; wx=wex';
        wey = -2*ry.*we; wy=wey';
        wexx = 4.*we.*rx.^2-2.*we.*rxx; wxx=wexx';
        weyy =4.*we.*ry.^2-2.*we.*ryy; wyy=weyy';
    end
end