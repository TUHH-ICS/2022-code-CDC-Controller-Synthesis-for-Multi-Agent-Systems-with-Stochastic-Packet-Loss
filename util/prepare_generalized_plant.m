function [sysD, sysC, sysP, nay, nau] = prepare_generalized_plant(G, R)
%PREPARE_GENERALIZED_PLANT Summary of this function goes here
%   Detailed explanation goes here

% The feedthrough term is assumed to be 0
[A, B, C, ~] = ssdata(ss(G));
npx = size(A,1);
npu = size(B,2);
npy = size(C,1);

%% Assemble the generalized plant
Ad    = A;
Ac    = zeros(npx);
Ap    = zeros(npx);
Bd_u  = B;
Bc_u  = zeros(npx,npu);
Bp_u  = zeros(npx,npu);
Bd_w  = B;
Bc_w  = zeros(npx,npu);
Bp_w  = zeros(npx,npu);
Cd_y  = zeros(npy,npx);
Cc_y  = -C;
Cp_y  = zeros(npy,npx);
Cd_z  = zeros(npy+npu,npx);
Cc_z  = zeros(npy+npu,npx);
Cp_z  = [ -C              ;
           zeros(npu,npx) ];
Dd_yu = zeros(npy,npu);
Dc_yu = zeros(npy,npu);
Dp_yu = zeros(npy,npu);
Dd_yw = zeros(npy,npu);
Dc_yw = zeros(npy,npu);
Dp_yw = zeros(npy,npu);
Dd_zu = [ zeros(npy,npu); sqrtm(R) ];
Dc_zu = zeros(npy+npu,npu);
Dp_zu = zeros(npy+npu,npu);
Dd_zw = zeros(npy+npu,npu);
Dc_zw = zeros(npy+npu,npu);
Dp_zw = zeros(npy+npu,npu);

sysD = ss(Ad, [Bd_w, Bd_u], [Cd_z; Cd_y], [Dd_zw, Dd_zu; Dd_yw, Dd_yu], 1);
sysC = ss(Ac, [Bc_w, Bc_u], [Cc_z; Cc_y], [Dc_zw, Dc_zu; Dc_yw, Dc_yu], 1);
sysP = ss(Ap, [Bp_w, Bp_u], [Cp_z; Cp_y], [Dp_zw, Dp_zu; Dp_yw, Dp_yu], 1);
nay  = size(Cd_y,1);
nau  = size(Bd_u,2);
end
