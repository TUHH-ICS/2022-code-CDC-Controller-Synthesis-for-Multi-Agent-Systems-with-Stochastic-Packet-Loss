%---------------------------------------------------------------------------------------------------
% For Paper
% "Stochastic Packet Loss in Multi-Agent Systems: An Empirical and Theoretical Analysis"
% by C. Hespe, A. Datar, D. Schneider, H. Saadabadi, H. Werner and H. Frey
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function sys = addparts(sysD, sysC, lambda)
[Ad, Bd, Cd, Dd] = ssdata(sysD);
[Ac, Bc, Cc, Dc] = ssdata(sysC);
sys = ss(Ad+lambda*Ac, Bd+lambda*Bc, Cd+lambda*Cc, Dd+lambda*Dc, 1);     
end
