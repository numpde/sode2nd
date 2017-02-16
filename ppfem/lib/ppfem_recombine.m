function s = ppfem_recombine(M, s)
% function Ms = ppfem_recombine(M, s)
%
%   Computes a new ppfem space Ms
%   such that Ms.pp = M * s.pp

%   R. Andreev, 20150821

    assert(size(M,2) == s.dim, 'size(M,2) should equal s.dim');

    rows = num2cell(M,2)';
    
    s.pp = ...
        foreach(rows, ...
            @(row) pp0add( ...
                foreach(find(row), ...
                    @(j) pp0scale(s.pp{j}, full(row(j))))));
                
    s = ppfem_filldetails(s, s.x);
end

function fran = foreach(ran, f)
    if (isnumeric(ran)); ran = num2cell(ran); end
    fran = cellfun(f, ran, 'UniformOutput', false);
end
