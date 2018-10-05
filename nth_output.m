function [ value ] = nth_output(N, fcn, varargin )
%NTH_OUTPUT returns Nth output value of a function fcn
% In:
%    N      1 x 1   integer of desired output
%   fcn     fhandle 
% Out: 
%   value   any output
% Last modified: Jonas Umlauft 10/2018

  [value{1:N}] = fcn(varargin{:});
  value = value{N};
end

