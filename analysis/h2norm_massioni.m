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
