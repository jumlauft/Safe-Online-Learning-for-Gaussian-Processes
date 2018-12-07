function [K,dKdhyp] = covSEaddmultj(hyp,xi,xj,u,uj)
%COVSEADDJ linear combination of R Squared Exponential covariance functions
% with Automatic Relevance Detemination (ARD) distance measure.
%
%       k(xi,xj) = g_1(xi)*k_1(xi,xj)*g_1(xj) + g_2(xi)*k_2(xi,xj)*g_2(xj) + ...
% where
%     k_r(xi,xj) = sf_r^2 * exp(-0.5*sum_k((xi - xj)^2/lk_r^2))    r=1...R
%
% hyp = [ log(l1_1); log(l2_1);..; log(sf_1);log(l1_2); log(l2_2);..; log(sf_2);... ]
%
% In:
%    xi            E  x n           Data points
%    xj            E  x m           Data points
%    hyp           R*(E+1) x 1(1xE+1)   Hyperparameter
%    g             @{R x 1}          Cell array of function handles
% Out:
%    K             n  x m        Covariance Matrix
%    dKdhyp        n  x m x R*(E+1)  Derivative wrt hyp
%    dKdxi         n x m x E x n  Derivative wrt xi (not implemented yet)
%
% E    : input dimensions
% n,m  : number of points
% R    : number of kernels added up
% Copyright (c) by Jonas Umlauft under BSD License
% Last modified: Jonas Umlauft 10/2018


hyp = hyp(:);
[E, n] = size(xi);
R = length(hyp)/(E+1);

if nargin >3
    if length(u) ==1 && ~iscell(u)
        u = {u};
    end
    if R~=size(u,1) && R~=size(u,2)
        error('size mismatch');
    end
end
% Copy in case only one input is provided
if nargin < 3
    xj = xi;
end

% Verify Sizes
if  size(xj,1)~=E
    error('size mismatch');
end

m = size(xj,2);

Kall = zeros(n,m,R);
ximxj2il2 = zeros(n,m,E,R);
ximxj = xi-permute(xj,[1 3 2]);
ui = zeros(R,n);
for r = 1:R
    ell = exp(hyp((r-1)*(E+1)+1:(r-1)*(E+1)+E));  % characteristic length scale
    sf2 = exp(2*hyp((r-1)*(E+1)+E+1));
    
    if isinf(sf2) || any(isinf(ell)), error('hyp is inf'); end
    
    il2 = 1./(ell.^2);
    ximxj2il2(:,:,:,r) = permute((ximxj.^2).*il2,[2 3 1]);
    Kall(:,:,r) =  sf2 * exp(-sum(ximxj2il2(:,:,:,r),3)./2);
    if nargin >3
        if ~isnumeric(u)
            ui(r,:) = u{r}(xi);
            uj(r,:) = u{r}(xj);
        else
            ui=u;
            if isequal(xi,xj), uj=u;end
        end
        multi = ui(r,:)'.*uj(r,:);
        Kall(:,:,r) = multi.* Kall(:,:,r);
    end
end
K = sum(Kall,3);

% Derivatives wrt hyperparameter
if nargout > 1
    dKdhyp = zeros(n,m,(E+1)*R);
    for r = 1:R
        
        % Derivative wrt ell
        for i=1:E
            dKdhyp(:,:,(r-1)*(E+1)+i) = Kall(:,:,r).*ximxj2il2(:,:,i,r);
        end
        % Derivative wrt sf
        dKdhyp(:,:,(r-1)*(E+1)+E+1) = 2*Kall(:,:,r);
    end
end



