function ret = iszero(A, tol)
%ISZERO Summary of this function goes here
%   Detailed explanation goes here

if nargin <= 1
    tol = 1e-8;
end

ret = all(abs(A) < tol, 'all');
end
