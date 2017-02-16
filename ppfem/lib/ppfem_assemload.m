function F = ppfem_assemload(f, spaces)
% function F = ppfem_assemload(f, spaces)
%
%   spaces = {X1, X2, ...}
%   f is a function handle f(x1, x2, ...)
%
%   F is a Nx1 vector with
%       x1 the slowest index, 
%       x2 the 2nd-slowest, etc.

%   R. Andreev, 2015-08-21

    space = spaces{1}; 
    if (length(spaces) > 1)
        f = @(t) ppfem_assemload(@(varargin)f(t,varargin{:}), spaces(2:end));
    end
    
    F = [];
    for p = eachof(space.pp)
        sumover = @(index, g) sum(cell2mat(foreach(index, g)), 2);
        pieces = foreach(1:p.pieces, @(k) ppbrk(p, k));
        F = [F; 
            sumover(pieces, ...
                @(piece) sumover(space.QR(piece.breaks), ...
                	@(xw) f(xw.x) * ppval(piece, xw.x) * xw.w));
            ];
    end
end

function fran = foreach(ran, f)
    if (isnumeric(ran)); ran = num2cell(ran); end
    if (isstruct(ran)); apply = @arrayfun; else apply = @cellfun; end
    fran = apply(f, ran, 'UniformOutput', false);
end

function array = eachof(array)
    array = [array{:}];
end
