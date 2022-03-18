%---------------------------------------------------------------------------------------------------
% For Paper
% "Stochastic Packet Loss in Multi-Agent Systems: An Empirical and Theoretical Analysis"
% by C. Hespe, A. Datar, D. Schneider, H. Saadabadi, H. Werner and H. Frey
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function ret = iszero(A, tol)
%ISZERO Summary of this function goes here
%   Detailed explanation goes here

if nargin <= 1
    tol = 1e-8;
end

ret = all(abs(A) < tol, 'all');
end
