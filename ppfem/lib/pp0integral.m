function I = pp0integral(p)
% function I = pp0integral(p)
%
%   p is in ppform
%
%   Computes the integral I of p
%   from p.break(1) to p.break(end)

    I = ppval(fnint(p), p.breaks(end));
end
