function q = pp0rfn(p, br)
% function q = pp0rfn(p, br)
%
%   p is in ppform
%   br is a vector of knots
%
%   Same semantics as pprfn but
%   q is zero outside [min(p.break), max(p.break)]

%   R. Andreev, 20150819

    q = pprfn(p, br);
    % Make q = 0 outside the original support of p
    q.coefs(q.breaks < min(p.breaks), :) = 0;
    q.coefs(q.breaks(1:end-1) >= max(p.breaks), :) = 0;
end
