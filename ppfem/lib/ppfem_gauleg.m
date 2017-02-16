function QR = ppfem_gauleg(N, I)
% function QR = ppfem_gauleg(N, I)
%
%   N node Gauss-Legendre quadrature 
%   for each interval (I(k,1), I(k,2))
%   where the matrix I has size K x 2 
%
%   If I is not specified, returns 
%   a function handle QR such that setting
%       qr = QR(I); 
%   gives struct arrays
%       qr.x, qr.w
%   for which the row vectors
%       [qr(k,:).x], [qr(k,:).w]
%   are the nodes and weights on the k-th interval and
%       [qr.x], [qr.w]
%   are all the nodes and weights
%
%   If I is specified, returns qr

%   R. Andreev, 2015-08-19 -- 2016-08-20

    assert(((1 <= N) && (N <= 2^10)), 'Expect 1 <= N <= 2^10 for the number of nodes N');

    d = 1 ./ sqrt(4 - 1 ./ ([1:N-1].^2));
    [V, D] = eig(diag(d,-1) + diag(d,1)); 
    X = diag(D) / 2; W = abs(V(1,:))' .^2;
    QR = @(I) struct('x', num2cell((ones(size(X))*mean(I')+(X*diff(I')))'), 'w', num2cell((W*diff(I'))'));
    
    if (nargin > 1); QR = QR(I); end
end
