function E = ppfem_eval(space, X, optional_sign)
% function E = ppfem_eval(space, x, sign)
%
%   If u is a column vector of coefficients for ``space''
%   then E * u evaluates the function at nodes x = [x1, x2, ...]
%
%   If space contains discontinuous piecewise polynomials,
%   evaluation may fail at break points
%
%   sign is optional and can be 
%       '+' or +1   from above
%       '-' or -1   from below
%       'a'         average
%       'j'         jump (above - below)
%       0           ignored
%       otherwise   error

%   R. Andreev, 2015-08-19 -- 2016-04-22

    E = sparse(numel(X), space.dim);
    E(:) = cell2mat(foreach(space.pp, @(p) transpose(pp0val(p,X))));
    
    if (nargin >= 3)
        if (isequal(optional_sign, 0)); return; end
        spdiag = @(v) spdiags(v(:), 0, numel(v), numel(v));
        for x = X
            % Basis filter for evaluation from from Above / Below
            A = spdiag(double(cell2mat(foreach(space.pp, @(p)(max(p.breaks) > x)))));
            B = spdiag(double(cell2mat(foreach(space.pp, @(p)(min(p.breaks) < x)))));
            switch optional_sign
                case {'a'}
                    A = A/2; B = B/2;   % average
                case {-1, '-'}
                    A = 0;              % don't use values from above
                case {+1, '+'}
                    B = 0;              % don't use values from below
                case 'j'
                    B = -B;             % Compute jump = above - below
                otherwise
                    error('Invalid sign. Valid values are -, +, a, j, 0');
            end
            % Apply filter
            E(x==X,:) = (E(x==X,:) * A) + (E(x==X,:) * B); 
        end
    end
end

function fran = foreach(ran, f)
    fran = cellfun(f, ran, 'UniformOutput', false);
end
