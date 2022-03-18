function [Fd, Fc, H2, Q, solver_stats] = h2syn_lossy(sysD, sysC, sysP, ny, nu, L0, p, backoff)
%H2SYN_LOSSY Summary of this function goes here
%   Detailed explanation goes here

offset = 1e-8;
opts   = sdpsettings('verbose', 0);

solver_stats        = struct;
solver_stats.prep   = NaN;
solver_stats.solver = NaN;
solver_stats.trans  = NaN;
solver_stats.total  = NaN;
solver_stats.yalmip = NaN;

Fd    = [];
Fc    = [];
H2    = Inf;
Q     = [];
timer = tic;

if nargin <= 7
    backoff = 0.05;
elseif backoff < 0
    error('Backoff must be non-negative!')
end

%%
[Ad, Bd, Cd, Dd] = ssdata(sysD);
[Ac, Bc, Cc, Dc] = ssdata(sysC);
[Ap, Bp, Cp, Dp] = ssdata(sysP);

% Read performance channel size
nx = size(Ac,1);
nw = size(Bc,2) - nu;
nz = size(Cc,1) - ny;
nL = size(L0,1);

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

if ~iszero(Dd_yu) || ~iszero(Dc_yu) || ~iszero(Dp_yu)
    error('There may not be feedthrough from input to output (D_yu)!')
end

if ~iszero(Ap) || ~iszero(Bp_u) || ~iszero(Cp_y)
    error('The plant may not have deterministic interconnections outside of the performance channels!')
end

P_in  = ~iszero(Bc_u) || ~iszero(Bp_u) || ~iszero(Dc_zu) || ~iszero(Dp_zu);
P_out = ~iszero(Cc_y) || ~iszero(Cp_y) || ~iszero(Dc_yw) || ~iszero(Dp_yw);
if P_in && P_out
    error('The plant may not be distributed on both control input and measured output at the same time')
end

lambda = eig(L0);

%%
alpha = sdpvar(1);
beta  = sdpvar(1);
X     = sdpvar(nx, nx, 'symmetric');
Y     = sdpvar(nx, nx, 'symmetric');
Z     = sdpvar(nw, nw, nL-1, 'symmetric');
Kd    = sdpvar(nx, nx, 'full');
Kc    = sdpvar(nx, nx, 'full');
Ld    = sdpvar(nx, ny, 'full');
Lc    = sdpvar(nx, ny, 'full');
Md    = sdpvar(nu, nx, 'full');
Mc    = sdpvar(nu, nx, 'full');
Nd    = sdpvar(nu, ny, 'full');
Nc    = sdpvar(nu, ny, 'full');
cost  = 0;

if ~P_in && ~P_out
    Constraints = [];
elseif ~P_in
    Constraints = [ Nc == 0, Lc == 0 ];
else
    Constraints = [ Nc == 0, Mc == 0 ];
end

for i = 2:nL
    li  = lambda(i);
    lit = sqrt(p*(1-p)*li);
    Zi  = Z(:,:,i-1);
    
    A   = Ad    + p*li*Ac    + li*Ap;
    Bu  = Bd_u  + p*li*Bc_u  + li*Bp_u;
    Bw  = Bd_w  + p*li*Bc_w  + li*Bp_w;
    Cy  = Cd_y  + p*li*Cc_y  + li*Cp_y;
    Cz  = Cd_z  + p*li*Cc_z  + li*Cp_z;
    Dyw = Dd_yw + p*li*Dc_yw + li*Dp_yw;
    Dzu = Dd_zu + p*li*Dc_zu + li*Dp_zu;
    Dzw = Dd_zw + p*li*Dc_zw + li*Dp_zw;
    K   = Kd    + p*li*Kc;
    L   = Ld    + p*li*Lc;
    M   = Md    + p*li*Mc;
    N   = Nd    + p*li*Nc;
    
    lmi_51 = lit*(Ac*Y + Bc_u*Md + Bd_u*Mc);
    lmi_52 = lit*(Ac + Bc_u*Nd*Cd_y + Bd_u*Nc*Cd_y + Bd_u*Nd*Cc_y);
    lmi_62 = lit*(X*Ac + Lc*Cd_y + Ld*Cc_y);
    lmi_81 = lit*(Cc_z*Y + Dd_zu*Mc + Dc_zu*Md);
    lmi_82 = lit*(Cc_z + Dc_zu*Nd*Cd_y + Dd_zu*Nc*Cd_y + Dd_zu*Nd*Cc_y);
    trc_41 = lit*(Bc_w + Bc_u*Nd*Dd_yw + Bd_u*Nc*Dd_yw + Bd_u*Nd*Dc_yw);
    trc_51 = lit*(X*Bc_w + Lc*Dd_yw + Ld*Dc_yw);
    trc_71 = lit*(Dc_zw + Dc_zu*Nd*Dd_yw + Dd_zu*Nc*Dd_yw + Dd_zu*Nd*Dc_yw);
    
    LMI = [ Y           eye(nx)      Y*A'+M'*Bu'    K'            lmi_51'       lit*Kc'       Y*Cz'+M'*Dzu'    lmi_81'      ; 
            eye(nx)     X            A'+Cy'*N'*Bu'  A'*X+Cy'*L'   lmi_52'       lmi_62'       Cz'+Cy'*N'*Dzu'  lmi_82'      ;
            A*Y+Bu*M    A+Bu*N*Cy    Y              eye(nx)       zeros(nx)     zeros(nx)     zeros(nx,nz)     zeros(nx,nz) ;
            K           X*A+L*Cy     eye(nx)        X             zeros(nx)     zeros(nx)     zeros(nx,nz)     zeros(nx,nz) ;
            lmi_51      lmi_52       zeros(nx)      zeros(nx)     Y/2           eye(nx)/2     zeros(nx,nz)     zeros(nx,nz) ;
            lit*Kc      lmi_62       zeros(nx)      zeros(nx)     eye(nx)/2     X/2           zeros(nx,nz)     zeros(nx,nz) ;
            Cz*Y+Dzu*M  Cz+Dzu*N*Cy  zeros(nz,nx)   zeros(nz,nx)  zeros(nz,nx)  zeros(nz,nx)  eye(nz)          zeros(nz)    ;
            lmi_81      lmi_82       zeros(nz,nx)   zeros(nz,nx)  zeros(nz,nx)  zeros(nz,nx)  zeros(nz)        eye(nz)/2    ];
    TRC = [ Zi             Bw'+Dyw'*N'*Bu'  Bw'*X+Dyw'*L'  trc_41'       trc_51'       Dzw'+Dyw'*N'*Dzu'  trc_71'      ;
            Bw+Bu*N*Dyw    Y                eye(nx)        zeros(nx)     zeros(nx)     zeros(nx,nz)       zeros(nx,nz) ;
            X*Bw+L*Dyw     eye(nx)          X              zeros(nx)     zeros(nx)     zeros(nx,nz)       zeros(nx,nz) ;
            trc_41         zeros(nx)        zeros(nx)      Y/2           eye(nx)/2     zeros(nx,nz)       zeros(nx,nz) ;
            trc_51         zeros(nx)        zeros(nx)      eye(nx)/2     X/2           zeros(nx,nz)       zeros(nx,nz) ;
            Dzw+Dzu*N*Dyw  zeros(nz,nx)     zeros(nz,nx)   zeros(nz,nx)  zeros(nz,nx)  eye(nz)            zeros(nz)    ;
            trc_71         zeros(nz,nx)     zeros(nz,nx)   zeros(nz,nx)  zeros(nz,nx)  zeros(nz)          eye(nz)/2    ];
    cost = cost + trace(Zi);
    
    % Add additional conditions to improve numerical conditioning
    NUM_1 = [ K, L; M, N ];
    NUM   = [ alpha * eye(nx+nu), NUM_1; NUM_1', alpha * eye(nx+ny) ]; 
    
    try
        Constraints = [ Constraints                    ,...
                        LMI >= offset * eye(size(LMI)) ,...
                        TRC >= offset * eye(size(TRC)) ,...
                        NUM >= offset * eye(size(NUM)) ];
    catch ME
        warning(ME.message)
        solver_stats.total = toc(timer);
        return
    end
end

% Add additional conditions to improve numerical conditioning
Constraints = [ Constraints                                               ,...
                beta >= 1                                                 ,...
                X    <= alpha * eye(nx)                                   ,...
                Y    <= alpha * eye(nx)                                   ,...
                [ Y, beta*eye(nx); beta*eye(nx), X] >= offset * eye(2*nx) ];

%% Implement three stage procedure to improve numerical conditioning
solver_stats.prep = toc(timer);
sol = optimize(Constraints, cost, opts);

if backoff > 0
    if sol.problem ~= 0
        warning('YALMIP return an error: %s', sol.info)
        solver_stats.total = toc(timer);
        return
    end
    solver_stats.yalmip = sol.yalmiptime;
    solver_stats.solver = sol.solvertime;
    
    % Backoff is squared since the cost is gamma^2
    c_opt = (1+backoff)^2 * value(cost);
    sol   = optimize([Constraints, cost <= c_opt], alpha, opts);
    if sol.problem ~= 0
        warning('YALMIP return an error: %s', sol.info)
        solver_stats.total = toc(timer);
        return
    end
    solver_stats.yalmip = solver_stats.yalmip + sol.yalmiptime;
    solver_stats.solver = solver_stats.solver + sol.solvertime;
    
    a_opt = (1+backoff) * value(alpha);
    sol   = optimize([Constraints, cost <= c_opt, alpha <= a_opt], -beta, opts);
else
    solver_stats.yalmip = 0;
    solver_stats.solver = 0;
end

%% Calculate the controller from the final solution variables
if sol.problem ~= 0
    warning('YALMIP return an error: %s', sol.info)
else
    solver_stats.yalmip = solver_stats.yalmip + sol.yalmiptime;
    solver_stats.solver = solver_stats.solver + sol.solvertime;
    trans = tic;
    
    [U, Sigma, V] = svd(eye(nx) - value(X*Y));
    U = U*sqrtm(Sigma);
    V = V*sqrtm(Sigma);
    
    H2 = sqrt(value(cost));
    Pi = value([ Y   eye(nx)   ;
                 V'  zeros(nx) ]);
    Q = Pi' \ value([ Y, eye(nx); eye(nx), X ]) / Pi;
    
    Dd_K = value(Nd);
    Cd_K = value(Md - Nd*Cd_y*Y)/V';
    Bd_K = U\value(Ld - X*Bd_u*Nd);
    Ad_K = U\value(Kd - X*(Ad+Bd_u*Nd*Cd_y)*Y)/V' - Bd_K*value(Cd_y*Y)/V' - U\value(X*Bd_u)*Cd_K;
    Fd = ss(Ad_K, Bd_K, Cd_K, Dd_K);
    
    Dc_K = value(Nc);
    Cc_K = value(Mc - (Nc*Cd_y+Nd*Cc_y)*Y)/V';
    Bc_K = U\value(Lc - X*(Bc_u*Nd+Bd_u*Nc));
    Ac_K = U\value(Kc - X*(Ac+Bc_u*Nd*Cd_y+Bd_u*Nc*Cd_y+Bd_u*Nd*Cc_y)*Y)/V' ...
           - value((Bc_K*Cd_y+Bd_K*Cc_y)*Y)/V' - U\value(X*(Bc_u*Cd_K+Bd_u*Cc_K));
    Fc = ss(Ac_K, Bc_K, Cc_K, Dc_K);
    
    solver_stats.trans = toc(trans);
end

solver_stats.total = toc(timer);
end
