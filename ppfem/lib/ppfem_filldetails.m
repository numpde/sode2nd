function space = ppfem_filldetails(space, x)
% function space = ppfem_filldetails(space, x)
%
%   TODO

%   R. Andreev, 20150819

    space.dim = length(space.pp);
    space.QR = ppfem_gauleg(max(getfields(eachof(space.pp), 'order')));
    if (nargin < 2); x = sort(unique(getfields(eachof(space.pp), 'breaks'))); end
    space.x = x;
end

function array = eachof(array)
    array = [array{:}];
end

function f = getfields(s, name)
    f = [s(:).(name)];
end
