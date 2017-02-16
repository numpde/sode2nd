function space = ppfem_transform(space, f)
% function space = ppfem_transform(space, f)
%
%   Apply f to each basis function of space
%   The function handle f should accept/return a ppform

%   R. Andreev, 20150819

    space.pp = foreach(space.pp, @(p) f(p));
    space = ppfem_filldetails(space);
end

function fran = foreach(ran, f)
    fran = cellfun(f, ran, 'UniformOutput', false);
end
