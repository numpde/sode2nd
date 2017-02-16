function p1 = pp0add(p1, varargin)
% function p = pp0add(p1, p2, ...)
%
% 	pi are in ppform
%
%   At least one argument is expected
%
%   If p1 is of type cell then pp0add(p1{:}) is returned
%
%   Computes the sum of of p1, p2, p3, ..., with 
% 	a) p1 or p2 can be [] to mean the zero function
%   b) p1, p2 are zero outside their [min(break), max(break)]

%   R. Andreev, 20150820

    if (iscell(p1)); p1 = pp0add(p1{:}); return; end
    if (isempty(varargin)); return; end
    p1 = pp0add(pp0cmb(p1, '+', varargin{1}), varargin{2:end});
end
