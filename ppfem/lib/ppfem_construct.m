function space = ppfem_construct(x, P, type)
% function space = ppfem_construct(x, p, type)
%
%   x is a row vector defining the 1d mesh
%   p is the piecewise polynomial degree 
%
%   Constructs a space with basis functions:
%       B-splines if type == 'B' (default)
%       P^p space if type == 'P'
%       Local L2-orthogonal Legendre polynomials if type == 'L' or 'D'
%       Forward / backward facing Radau [1] if type = 'R' / 'R*'
%       Jumplets if type == 'J'
%       
%
%   Variants for type
%       'B0': "Open" B-splines (knots don't accumulate at boundary)
%       'L1': Like 'L' with Legendre polynomials normalized to L(right) = 1
%
%   The ordering is
%    - In case of 'B': from left to right
%    - In case of 'P': hats first, then interval-wise
%    - In case of 'L': interval-wise
%
%   [1] See Andreev/Schweitzer, ETNA, 41 (2014)

%   R. Andreev, 2015-08-19 -- 2016-08-18

    if (nargin < 3); type = 'B'; end
    assert(P >= 0, 'Piecewise polynomial degree should be nonnegative');
    
    switch type
        case 'B'
            % Knots for all the bsplines
            K = [x(1)*ones(1,P), x, x(end)*ones(1,P)];
            % Construct bsplines
            space.pp = foreach(1:(length(K)-P-1), @(a) bspline(K(a+[0:P+1])));
        case {'B0', 'BO'} % B-zero or B-Oh
            K = [x(end+[-P:-1])-x(end)+x(1), x, x(1+[1:P])-x(1)+x(end)];
            space.pp = foreach(1:(length(K)-P-1), @(a) bspline(K(a+[0:P+1])));
            space = ppfem_transform(space, @(p) fnbrk(p, [max(min(p.breaks),min(x)), min(max(p.breaks),max(x))]));
        case 'P'
            space = ppfem_construct(x, min(P,1), 'B');
            intervals = foreach(struct('a', num2cell(x(1:end-1)), 'b', num2cell(x(2:end))), @(i) [i.a, i.b]);
            space.pp = flatten({space.pp, foreach(intervals, @(I) BabuskaShen(I, P))});
        case {'L', 'D'}
            intervals = foreach(struct('a', num2cell(x(1:end-1)), 'b', num2cell(x(2:end))), @(i) [i.a, i.b]);
            space.pp = flatten(foreach(intervals, @(I) Legendre(I, P))');
        case 'L1'
            space = ppfem_construct(x, P, 'L');
            space.pp = foreach(space.pp, @(p) pp0scale(p, 1 / fnval(p, p.breaks(end))));
        case {'R', 'R*'}
            % Radau space (adjoint)
            space = ppfem_construct(x, P, 'L');
            if (isequal(type, 'R')); f = 1; else f = -1; end
            N = (length(x) - 1); % Number of elements
            R = sparse(sort([1:(N*P), (1:N)*P]), 1:(N*(P+1)), repmat([ones(1,P), f * sqrt((2*P+1)/(2*P-1))], [1, N]));
            space = ppfem_recombine(R, space);
        case 'J'
            % Jumplets
            L = @(I, x) pp0add(foreach(Legendre(I, P), @(p) pp0scale(p, fnval(p, x))));
            space.pp = flatten(foreach([2:(length(x)-1)], @(n) pp0add(L(x([n-1,n]), x(n)), pp0scale(L(x([n,n+1]), x(n)), -1))));
            space.pp = foreach(space.pp, @(p) pp0scale(p, 1 / sqrt(pp0integral(pp0cmb(p, '*', p)))));
        otherwise
            error('Invalid basis type specified');
    end
    
    space = ppfem_filldetails(space, x);
end

function B = flatten(A)
% Flatten a cell array of possible cell arrays
    B = {};
    for i = 1:numel(A)  
        if(iscell(A{i}))
            F = flatten(A{i});
            B = [B, F{:}];
        else
            B = [B, A{i}];
        end
    end
end

function pp = BabuskaShen(I, P)
% Babuska--Shen basis on the interval I = [a, b]
% starting with polynomial degree 2 up to polynomial degree P >= 2
    if (P < 2); pp = {}; return; end
    pp = Legendre(I, P-1);
    pp = foreach(pp(2:end), @fnint);
end

function pp = Legendre(I, P)
% L2-orthogonal Legendre polynomial basis on the interval I = [a, b]
% up to and including polynomial degree P
    assert(all(size(I) == [1, 2]), 'The interval should be of the form I = [a, b]');
    if (P < 0); pp = {}; return; end
    
    h = diff(I);
    pp = {mkpp(I, 0), mkpp(I, 1/sqrt(h))};
    for p = 1:P
        pp{end+1} = pp0scale(pp0cmb(pp0scale(pp0cmb(pp{end}, '*', mkpp(I, [2/h, -1])), sqrt(2*p-1)/p*sqrt(2*p+1)), '+', pp0scale(pp{end-1}, -(p-1)/p/sqrt(2*p-3)*sqrt(2*p+1))), 1);
    end
    pp(1) = [];
    
%     % Test orthogonality:
%     space.pp = pp;
%     M = ppfem_gramian(ppfem_filldetails(space));
%     assert(max(max(abs(M - speye(size(M))))) <= 1e-6);
end
    
function fran = foreach(ran, f)
    if (isnumeric(ran)); ran = num2cell(ran); end
    if (isstruct(ran)); apply = @arrayfun; else apply = @cellfun; end
    fran = apply(f, ran, 'UniformOutput', false);
end
