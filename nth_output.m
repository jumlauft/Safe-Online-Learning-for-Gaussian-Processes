function out = nth_output(N, fcn, varargin)
%NTH_OUTPUT returns Nth output value of a function fcn
% In:
%    N      1 x 1   integer of desired output
%   fcn     fhandle
% Out:
%   out     N-th output of fcn
% Copyright (c) by Jonas Umlauft under BSD License
% Last modified: Jonas Umlauft 10/2018

[out{1:N}] = fcn(varargin{:});
out = out{N};
end

