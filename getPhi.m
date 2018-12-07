function phi = getPhi(t,x,p,reffun,r_min)
%getPhi Computes Phi=sigma(x)/r(x)
%   where x sigma is GP standard deviation and r = [lambda 1](x-xd)
% In: 
%   t           1 x N   time
%   x           E x N   state
%   p  
%    .lam       E-1 x 1 filtering gain
%     .varfun   @fun    GP variance function
%   reffun      @fun    reference trajectory
%   r_min       1 x 1   small constant to avoid numerical difficulties
% Out: 
%      phi      N x 1   ratio sigma/r
%
% Copyright (c) by Jonas Umlauft under BSD License 
% Last modified: Jonas Umlauft 10/2018

xd = reffun(t);   
r = [p.lam' 1]*(x-xd(1:end-1,:));
r = max(r_min,r);
phi = sqrt(p.varfun(x))./abs(r)';
end

