function u = ctrlFeLi(t,x,p,reffun)
%CTRLFELI Feedback Linearizating Controller
% IN: 
%   t       1 x 1   current time
%   x       E x N   current state
%   p               Parameter struct
%     .f    @fun    Estimate for f
%     .g    @fun    Estimate for g
%     .kc   1 x 1   controller gain
%     .lam  E-1 x 1 filter coefficient (must be Hurwitz)
%   ref     @fun    reference trajectory
% OUT: 
%   u  1 x 1
% E: state space dimension
% Copyright (c) by Jonas Umlauft under BSD License
% Last modified: Jonas Umlauft 10/2018

E = size(x,1);
xd = reffun(t); 
e = x-xd(1:E);

r = [p.lam' 1]*e;
rho = p.lam'*e(2:E) - xd(E+1);
nu = -p.kc*r - rho;
u = 1./p.g(x) .* ( -p.f(x) +nu);
end