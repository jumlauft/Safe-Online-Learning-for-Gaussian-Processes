function [dfdx] = gradestj(fun,x0, eps)
% GRADESTJ Computes numerical gradient of scalar function
% In:
%    x0      E x N      Point where gradient computed
%    fun     fhandle    function  E x 1 -> scalar
%    eps     scalar     distance between two points for slope calculation
% Out:
%  dfdx      E x N
% Copyright (c) by Jonas Umlauft under BSD License 
% Last modified: Jonas Umlauft 10/2018
[E, N] = size(x0);
if ~exist('eps','var'), eps = 1e-3; end

xpme = repmat(x0,1,E,2) + kron(eye(E),ones(1,N)).*permute([eps -eps],[1 3 2]);
V = fun(reshape(xpme,E,E*N*2))';
dfdx = reshape((V(1:E*N) - V(N*E+1:end))./(2*eps),N,E)';
end
