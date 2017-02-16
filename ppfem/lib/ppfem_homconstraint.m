function [space, R, V] = ppfem_homconstraint(space, E)
% function [newspace, R, V] = ppfem_homconstraint(oldspace, E)
%
%   Construct a subspace newspace defined by the condition
%       E * b = 0
%   for any (column) vector of coefficients b for oldspace
%
%   E is a matrix with size2 = #columns = oldspace.dim
%
%   The subspace newspace is obtained
%   by transforming oldspace with R:
%       newspace = ppfem_recombine(R, oldspace);
%
%   Hence, R' is the prolongation operator:
%       oldcoeff = R' * newcoeff
%
%   Each column v of V is a vector of coefficients
%   with respect to oldspace such that
%       (E * v)_i = 1 for at least one i 
%   and E(i,:) represents a nonredundant constraint
%
%   The entries e of E satisfying
%       |e| < 1e-10 * max|E(:)|
%   are rounded to zero

%   R. Andreev, 20150819

    if (nnz(E) == 0); R = 1; V = zeros(space.dim, 0); return; end
    assert(size(E,2) == space.dim, 'Number of columns of E should equal space.dim');
    
    E(abs(E) < 1e-10 * max(abs(E(:)))) = 0;
    E_old = E;
    
    % Sort in descending order
    E = sortrows(E, -(1:size(E,2)));
    
    % Reduced row echelon form represents the same constraints
    % jb contains indices of the bound variables
    [E, jb] = rref(E);
    if (issparse(E_old)); E = sparse(E); end
    
    % Rank
    r = length(jb);
    % Remove trivial constraints
    E = E(1:r, :);
    
    % Rescale columns of V (cf. function description)
    V = pinv(full(E));
    [~, I] = sort(abs(E_old * V), 1, 'descend');
    V = V * diag(1 ./ diag(E_old(I(1,:),:) * V));
    
    % jb is increasing (due to sort of E above)
    assert(all(diff(jb) > 0));
    
    % post: E(:, jb) == Id;
    
    R = speye(space.dim);
    E(:, jb) = 0;
    R(:, jb) = -E';
    R(jb, :) = [];
    
    % Impose the homogeneity constraint in the new space:
    % 	newspace.pp = R * oldspace.pp
    % Note that R' is the prolongation: 
    %   oldcoeff = R' * newcoeff

    space = ppfem_recombine(R, space);
end
