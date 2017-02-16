function p = pp0scale(p, lambda)
% function q = pp0scale(p, lambda)
%
%   Scale ppform p by scalar lambda

%   R. Andreev, 20150819

    assert(isnumeric(lambda), 'lambda should be a number');
	p = fncmb(p, lambda);
end
