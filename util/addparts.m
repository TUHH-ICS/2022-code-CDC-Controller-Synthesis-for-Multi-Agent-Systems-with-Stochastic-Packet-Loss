function sys = addparts(sysD, sysC, lambda)
[Ad, Bd, Cd, Dd] = ssdata(sysD);
[Ac, Bc, Cc, Dc] = ssdata(sysC);
sys = ss(Ad+lambda*Ac, Bd+lambda*Bc, Cd+lambda*Cc, Dd+lambda*Dc, 1);     
end
