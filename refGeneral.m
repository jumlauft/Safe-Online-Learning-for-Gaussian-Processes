function r = refGeneral(t,E,reffun)
%refGeneral General reference function
%   returns up to the n-th derviative of a refence fun to all t
%   e.g.for sinosoidal call refGeneral(t,E+1,@(x) sin(x))
% In:
%    t       N x 1 (1 x N)      time steps
%    E       1 x 1              number of derivatives
%    reffun  fhandle            Reference function 
% Out:
%    r       E x N              Reference and E derivatives for all t

% Copyright (c) by Jonas Umlauft under BSD License 
% Author: Jonas Umlauft, Last modified: 11/2017

if ~isscalar(E) || ~isvector(t)
    error('wrong input dimension');
end
t = t(:)';
r = zeros(E,length(t));
r(1,:) = reffun(t);
for e = 2:E
    reffun = @(x) gradestj(reffun,x);
    r(e,:) = reffun(t);
end


end
