function  [value,isterminal,direction]  = eventPhi(t,x,reffun,beta,p,r_min)
%EVENTPHI ode event to stop simulation if variance of GP is too large, 
%   more formally, the event is triggered if
%            Phi = sigma/r  >   kc/beta; 
%   The event is also triggered if |r|<r_min
% IN: 
%   t       1 x 1       time
%   x       E x 1       state
%   varfun  @fun     variance function 
%   reffun  @fun     reference function 
%   beta    1 x 1       scaling factor
%   p. 
%     kc    1 x 1       control gain
%     lam   E-1 x 1     filtering gain
% OUT: 
%  see ode event documentation

% Copyright (c) by Jonas Umlauft under BSD License 
% Last modified: Jonas Umlauft 10/2018


value(1) =  getPhi(t,x,p,reffun,r_min) - p.kc/beta;
isterminal = 1;
direction = 0;
end

