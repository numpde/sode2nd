function p = pp0cmb(p1, op, p2)
% function p = pp0cmb(p1, op, p2)
%
% 	p1 and p2 are in ppform
%
%   Same semantics as fncmb(p1, op, p1) but
% 	a) p1 or p2 can be [] to mean the zero function
%   b) p1, p2 are zero outside [min(break), max(break)]

%   R. Andreev, 20150819


    assert(nargin == 3, 'Expected three arguments: p1, op, p2');
    assert(~(isempty(p1) & isempty(p2)), 'p1 and p2 cannot both be empty');
    if (isempty(p1)); p1 = fncmb(p2, 0); end
    if (isempty(p2)); p2 = fncmb(p1, 0); end
    breaks = unique([p1.breaks, p2.breaks], 'sorted');
    p = fncmb(pp0rfn(p1, breaks), op, pp0rfn(p2, breaks));
end

function p = fncmb_(p1, op, p2)
% function fncmb_(p1, op, p2)
% 
%   Like fncmb but with less checks
%   and assuming that p1 and p2 are in pp form
%   with the same break sequence

    assert(nargin == 3);
    assert(all(size(p1.breaks) == size(p2.breaks)));
    assert(all(p1.breaks == p2.breaks));

    breaks = p1.breaks;
    
    c1 = p1.coefs;
    c2 = p2.coefs;
    assert(size(c1,1) == size(c2,1));
    c1 = [zeros(size(c1,1), size(c2,2)-size(c1,2)), c1];
    c2 = [zeros(size(c2,1), size(c1,2)-size(c2,2)), c2];

    switch op
        case '+'
            c = c1 + c2;
        case '-'
            c = c1 - c2;
%         case '*'
%               warning('May be wrong.');
%             c = [];
%             for i = 1 : size(c1,1)
%                 c = [c; conv(c1(i,:), c2(i,:))];
%             end
        otherwise
            p = fncmb(p1, op, p2);
            return;
    end
    
    p = mkpp(breaks, c);
end
