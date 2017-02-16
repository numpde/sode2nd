function M = ppfem_gramian(U, V, w)
% function M = ppfem_gramian(U, V, w)
%
%   Construct the Gramian for the spaces U and V
%
%   The function handle w to
%   the integration weight w : IR --> IR
%   is optional and defaults to w = 1

%   R. Andreev, 2015-08-19 -- 2016-08-20

    if (nargin < 2); V = U; end
    
    % Find the sparsity pattern of the mass matrix
    pp = {}; I = []; J = [];
    for i = 1 : U.dim; pi = U.pp{i};
    for j = 1 : V.dim; pj = V.pp{j};
        if (pi.breaks(end) <= pj.breaks(1)); continue; end
        if (pj.breaks(end) <= pi.breaks(1)); continue; end
        pp{end+1} = pp0cmb(pi, '*', pj);
        I(end+1) = i;
        J(end+1) = j;
    end
    end
    
    if (isempty(pp))
        M = sparse(U.dim, V.dim); 
        return; 
    end
    
    % Default integration weight w = 1
    if (nargin <= 2)
        pp = cellfun(@(p)pp0integral(p), pp, 'UniformOutput', false);
        M = sparse(I, J, [pp{:}], U.dim, V.dim);
        return;
    end
    
    % General integration weight 
    Z.x = sort(unique(union(U.x, V.x)));
    Z.pp = pp;
    Z = ppfem_filldetails(Z);
    M = sparse(I, J, ppfem_assemload(w, {Z}), U.dim, V.dim);
end
