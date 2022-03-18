%---------------------------------------------------------------------------------------------------
% For Paper
% "Stochastic Packet Loss in Multi-Agent Systems: An Empirical and Theoretical Analysis"
% by C. Hespe, A. Datar, D. Schneider, H. Saadabadi, H. Werner and H. Frey
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function [H2, time] = h2norm_massioni(sysD, sysC, Kd, Kc, L)
%H2NORM_MASSIONI Summary of this function goes here
%   Detailed explanation goes here

timer  = tic;

lambda = eig(L);
H2     = 0;
for i = 2:length(lambda)
    sys = addparts(sysD, sysC, lambda(i));
    K   = addparts(Kd, Kc, lambda(i));
    CL  = lft(sys, K);
    H2  = H2 + norm(CL)^2;
end

H2   = sqrt(H2);
time = toc(timer);
end
