function y = pp0val(p, x)
% function y = pp0val(p, x)
%
%   p is in ppform
%   x is a row vector of nodes
%
%   For column valued functions p only
%
%   Like ppval or fnval, but 
%   is zero outside [min(p.breaks), max(p.breaks)]

    if (isempty(p)); y = 0*x; return; end
    
    y = ppval(p, x);
    y(:, ~[(min(p.breaks) <= x) & (x <= max(p.breaks))]) = 0;
end
