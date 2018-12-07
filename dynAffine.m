function dxdt = dynAffine(t,x,ctrl,p)
%dynAffine Dynamics of arbitrary input affine System of the form
%           dx_1 dt = x(2)
%           dx_2 dt = x(3)
%                ....
%           dx_E dt = f(x) + g(x)u
% IN: 
%   t     1 x 1
%   x     E x 1 
%   ctrl  @fun 
%   p               Parameter struct
%    .f   @fun      function handle
%    .g   @fun      function handle
% OUT: 
%   dxdt  E x 1
% Copyright (c) by Jonas Umlauft under BSD License
% Last modified: Jonas Umlauft 10/2018
E = size(x,1);
dxdt = zeros(E,1);
dxdt(1:E-1) = x(2:E);
dxdt(E) = p.f(x) + p.g(x)*ctrl(t,x);