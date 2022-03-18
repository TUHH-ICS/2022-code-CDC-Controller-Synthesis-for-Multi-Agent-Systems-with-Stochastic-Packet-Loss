function G = assemble_mas(sysD, sysC, sysP, ny, nu, L, L0)
%ASSEMBLE_MAS Summary of this function goes here
%   Detailed explanation goes here
hat = @(K, M, N) kron(eye(size(L)), K) + kron(L, M) + kron(L0, N);

[Ad, Bd, Cd, Dd] = ssdata(sysD);
[Ac, Bc, Cc, Dc] = ssdata(sysC);
[Ap, Bp, Cp, Dp] = ssdata(sysP); 

nw = size(Bc,2) - nu;
nz = size(Cc,1) - ny;

Bd_w  = Bd(:,1:nw);
Bc_w  = Bc(:,1:nw);
Bp_w  = Bp(:,1:nw);
Bd_u  = Bd(:,nw+1:end);
Bc_u  = Bc(:,nw+1:end);
Bp_u  = Bp(:,nw+1:end);
Cd_z  = Cd(1:nz,:);
Cc_z  = Cc(1:nz,:);
Cp_z  = Cp(1:nz,:);
Cd_y  = Cd(nz+1:end,:);
Cc_y  = Cc(nz+1:end,:);
Cp_y  = Cp(nz+1:end,:);
Dd_zw = Dd(1:nz,1:nw);
Dc_zw = Dc(1:nz,1:nw);
Dp_zw = Dp(1:nz,1:nw);
Dd_zu = Dd(1:nz,nw+1:end);
Dc_zu = Dc(1:nz,nw+1:end);
Dp_zu = Dp(1:nz,nw+1:end);
Dd_yw = Dd(nz+1:end,1:nw);
Dc_yw = Dc(nz+1:end,1:nw);
Dp_yw = Dp(nz+1:end,1:nw);
Dd_yu = Dd(nz+1:end,nw+1:end);
Dc_yu = Dc(nz+1:end,nw+1:end);
Dp_yu = Dp(nz+1:end,nw+1:end);

G = ss(hat(Ad, Ac, Ap), [hat(Bd_w, Bc_w, Bp_w), hat(Bd_u, Bc_u, Bp_u)],...
       [hat(Cd_z, Cc_z, Cp_z); hat(Cd_y, Cc_y, Cp_y)],...
       [hat(Dd_zw, Dc_zw, Dp_zw), hat(Dd_zu, Dc_zu, Dp_zu);
        hat(Dd_yw, Dc_yw, Dp_yw), hat(Dd_yu, Dc_yu, Dp_yu)], 1);
end
