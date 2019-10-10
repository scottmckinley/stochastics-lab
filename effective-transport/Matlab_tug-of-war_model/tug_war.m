%% Tug-of-war model with maximum 2 plus and 2 minus motors

clear
close all

global eps0m eps0p Fdm Fdp Fsp Fsm vbm vfm vbp vfp pi0p pi0m Fext

% max number of plus and minus motors
Np = 2;
Nm = 2;

ind_diff = 0; % 0 if paused stationary state, 1 if diffusive
D        = 0.0; % diffusion coefficient

% Indices for motors modeled
indp = 1; % plus motor
indm = 1; % minus motor

Fs1 = 6; %[0.5:0.5:10];

% Kinesin 1, Dynein (Ref 22)
Fs_vect   = [Fs1    1.1  ];   % pN, stall force
Fd_vect   = [3      0.75 ];   % pN, detachment force
eps0_vect = [1      0.27 ];   % 1/s, unbinding rate
pi0_vect  = [5      1.6  ];   % h 1/s, binding rate
vf_vect   = [1      0.65 ];   % micron/s, forward velocity
vb_vect   = [0.006  0.072];   % micron/s, backward velocity

% Extract plus motor parameters
Fsp   = Fs_vect(indp);   Fdp  = Fd_vect(indp);
eps0p = eps0_vect(indp); pi0p = pi0_vect(indp);
vfp   = vf_vect(indp);   vbp  = vb_vect(indp);

% Extract minus motor parameters
Fsm   = Fs_vect(indm);   Fdm  = Fd_vect(indm);
eps0m = eps0_vect(indm); pi0m = pi0_vect(indm);
vfm   = vf_vect(indm);   vbm  = vb_vect(indm);

%% Set up transition matrix P

P = zeros((Np+1)*(Nm+1),(Np+1)*(Nm+1));
[v1,v2,v3,v4] = vectors();
P = P + diag(v1,-3) + diag(v2,3) + diag(v3,-1) + diag(v4,1);

% Set up absorbing state matrix Ptilde/Q
Pabs = P;
Pabs(1,:) = [1,zeros(1,(Np+1)*(Nm+1)-1)];
Q = Pabs(2:end,2:end);
Pinit = P(1,2:end);

% External force
Fext = 0;

%% Set up necessary moment vectors

% Set up speeds and rates out of each transient state
v_states = [vc(1,0) vc(2,0) vc(0,1) vc(1,1) vc(2,1) vc(0,2) vc(1,2) vc(2,2)];
rout_states = [                 epsp(1,0)+pip(2,1)+pim(2,0)           epsp(2,0)+pim(2,0)           ...
    epsm(0,1)+pip(2,0)+pim(2,1) epsm(1,1)+pip(2,1)+pim(2,1)+epsp(1,1) epsm(2,1)+pim(2,1)+epsp(2,1) ...
    epsm(0,2)+pip(2,0)          epsm(1,2)+pip(2,1)+epsp(1,2)          epsm(2,2)+epsp(2,2)
    ];
% Speed and rate out of the absorbing state
rout_state0 = pim(2,0) + pip(2,0);
v_state0    = 0;

% First, cross and second moments of times and rewards in each transient
% state
rvect  = v_states./rout_states;
tvect  = 1./rout_states;
rtvect = 2.*v_states./(rout_states.^2);
r2vect = 2.*(v_states.^2)./(rout_states.^2);
t2vect = 2./(rout_states.^2);

% Diagonal matrices of moments
diagr = diag(rvect);
diagt = diag(tvect);


%% Calculate moments given start in each transient state

% Identity matrix
imat = eye((Np+1)*(Nm+1)-1);

% First moments
imatminq = imat-Q;
uvectr = imatminq\rvect';
uvectt = imatminq\tvect';
% Second moments
vvectr = imatminq\(r2vect'+2.*diagr*Q*uvectr);
vvectt = imatminq\(t2vect'+2.*diagt*Q*uvectt);
% Cross moment
uvectrt = imatminq\(rtvect'+diagr*Q*uvectt+diagt*Q*uvectr);

%% Final calculations using classical theory

% Effective velocity
ex = sum(0 + uvectr'.*Pinit);
et = 1/rout_state0 + sum(uvectt'.*Pinit);
Vinf = ex/et;

et0x0 = 0;
ex0   = 0;
et0   = (1/rout_state0);
if ind_diff == 0
    ex02 = 0;               % if fully paused/stationary state, 0 second moment
else
    ex02 = 2*D/rout_state0; % if diffusive state, non-zero second moment
end
ext = et0x0 + ex0*dot(uvectt,Pinit) + et0*dot(uvectr,Pinit) + dot(uvectrt,Pinit);
ex2 = ex02 + 2*ex0*dot(uvectr,Pinit) + dot(vvectr,Pinit);
et2 = 2/rout_state0^2 + 2*(1/rout_state0)*dot(uvectt,Pinit) + dot(vvectt,Pinit);

varx  = ex2 - ex^2;
vart  = et2 - et^2;
vartx = ext - ex*et;

% Effective diffusivity
Dinf = (Vinf^2*vart+varx-2*Vinf*vartx)/(2*et);

%% Display effective velocity and diffusivity

disp(['Effective velocity is ' num2str(Vinf)]);
disp(['Effective diffusivity is ' num2str(Dinf)]);

