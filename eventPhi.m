function  [value,isterminal,direction]  = eventPhi(t,x,reffun,beta,p,r_min)
%EVENTPHI ode event to stop simulation if variance of GP is too large, 
%   more formally, the event is triggered if
%            Phi = sigma/r  >   kc/beta; 
%   The event is also triggered if |r|<r_min
% IN: 
%   t       1 x 1       time
%   x       E x 1       state
%   varfun  fhandle     variance function 
%   reffun  fhandle     reference function 
%   beta    1 x 1       scaling factor
%   p. 
%     kc    1 x 1       control gain
%     lam   E-1 x 1     filtering gain
% OUT: 
%  see ode event documentation

% Copyright (c) by Jonas Umlauft under BSD License 
% Last modified: Jonas Umlauft 10/2018

E = size(x,1);
xd = reffun(t);
r = [p.lam' 1]*(x-xd(1:E));

value(1) =  getPhi(t,x,p,reffun) - p.kc/beta;
value(2) = abs(r)-r_min; 


isterminal =[1 1];
direction = [0 0];
end

