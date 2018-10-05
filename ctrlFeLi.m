function u = ctrlFeLi(t,x,p,reffun)
%ctrlFeLi Feedback Linearization Controller with PD Controller
% IN: 
%   t     1 x 1
%   x     E x N   position 1..E, velcoities E+1...2*E
%   p               Parameter struct
%    .f  @fun
%    .g  @fun
%    .k  E x 1
%   ref  @fun 
% OUT: 
%   u  1 x 1
% Author: Jonas Umlauft, Last modified: 02/2017
E = size(x,1);
xd = reffun(t); 
e = x-xd(1:E);
% u = 1./p.g(x) .* ( -p.f(x) + xd(E+1) - p.k'*e);

r = [p.lam' 1]*e;
rho = p.lam'*e(2:E) - xd(E+1);
nu = -p.kc*r - rho;
u = 1./p.g(x) .* ( -p.f(x) +nu);


end