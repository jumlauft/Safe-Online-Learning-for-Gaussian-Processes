function [grid, max_dist] = ndgridj(grid_min, grid_max,ns)
% generates a grid and returns all combnations of that grid in a list
% In:
%   grid_min   1  x D  lower bounds of grid for each dimension separately
%   grid_max   1  x D  upper bounds of grid for each dimension separately
%   ns         1  x D  number of points for each dimension separately
% Out:
%   grid       Prod(ns) x D
%   max_dist   1  x 1   maximum distance a point can have in the grid
%
% Copyright (c) by Jonas Umlauft under BSD License
% Last modified: Jonas Umlauft 10/2018


D = numel(ns);

if isscalar(grid_max), grid_max = repmat(grid_max,1,D); end
if isscalar(grid_min), grid_min = repmat(grid_min,1,D); end

if numel(grid_max) ~= D ||  numel(grid_min) ~= D
    error('grid_max, grid_min and ns must have same dimensions');
end

if ~isvector(grid_max)  ||  ~isvector(grid_min) || ~isvector(ns)
    error('grid_max, grid_min and ns must be vectors');
end

if any(grid_max-grid_min<0)
    error('grid_max ist not always larger than grid_min');
end

gg = cell(1,D);

for i=1:D
    gg{i} = linspace(grid_min(i),grid_max(i),ns(i));
end
[gg{:}] = ndgrid(gg{:});

grid = reshape(permute(cat(D+1, gg{:}),[D:-1:1 D+1]),[],D)';

if nargout >1
    dist = (grid_max-grid_min)./ns;
    max_dist = sqrt(sum(dist.^2))/2;
end


