% Eissa et al, The Cross-Scale Effects of Neural Interactions during Human 
% Neocortical Seizure Activity, PNAS, 2017

%% Discretization

dt = 0.001;         % time step size [s]
T = 10.5/dt;      

N = 90;             
h = 9/5*4/(N-1);    % spatial discretization step size [mm]

%% Parameters

% Synaptic decay
alpha_E = 125;      % [s^-1]
alpha_I = 25;       % [s^-1]

% Firing rate
theta_E = 14;
theta_I = 9;
zeta_E = 6;
zeta_I = 3;

% Level of Feedforward inhibition
gamma = 1/2;

% Connectivity
eta_EE = 125;
eta_EI = -140;
eta_IE = 120;
eta_II = -62.5;

mu_EE = 0.4;        % [mm]
mu_EI = 0.5;        % [mm]
mu_IE = 0.5;        % [mm]
mu_II = 0.4;        % [mm]

%% Connection Matrices

C = zeros(N^2);

for i = 1:N^2
    
    for j = 1:N^2
        
        hd = abs(floor((i-1)/N)-floor((j-1)/N));
        vd = abs(mod(i-1,N)-mod(j-1,N));
        
        C(i,j) = sqrt(hd^2+vd^2)*h;
        
    end
    
end

for j = 1:N
    
    % Correction for boundary elements
    C(:,j) = 1/2*C(:,j);
    C(:,(N-1)*N+j) = 1/2*C(:,(N-1)*N+j);
    C(:,(j-1)*N+1) = 1/2*C(:,(j-1)*N+1);
    C(:,j*N) = 1/2*C(:,j*N);
    
end

J_EE = eta_EE*exp(-C./mu_EE);
J_EI = eta_EI*exp(-C./mu_EI);
J_IE = eta_IE*exp(-C./mu_IE);
J_II = eta_II*exp(-C./mu_II);

%% Preallocation

U_E = zeros(N^2,1);
U_I = zeros(N^2,1);
U_Eo = zeros(N^2,1);
U_Io = zeros(N^2,1);

u_E = zeros(1,T);
u_I = zeros(1,T);

SP = zeros(50,50,T);
LFP = zeros(50,50,T);

u_E(1:2) = 0.1502;
u_I(1:2) = 0.0655;

figure()

%% Time stepping

ce = alpha_E*dt;
ci = alpha_I*dt;

for t = 2:T
    t
    if t*dt > 1 && t*dt < 1.01

        W = u_E(t)*ones(N);
        W(end-9:end,1:10) = 10*ones(10);
        
    else
        
        W = u_E(t)*ones(N);
        
    end
    
        
    R_E = h^2*(J_EE*U_E+J_EI*U_I)+W(:);
    R_I = h^2*(J_IE*U_E+J_II*U_I)+gamma*W(:);
    
    r_E = 10*u_E(t)-14*u_I(t)+20*mean(mean(U_E))+5; 
    r_I = 12*u_E(t)-5*u_I(t)+gamma*(20*mean(mean(U_E))+5);

    F_E = exp(-((R_E-theta_E)./zeta_E).^2);
    F_I = exp(-((R_I-theta_I)./zeta_I).^2);
    
    f_E = exp(-((r_E-theta_E)./zeta_E).^2);
    f_I = exp(-((r_I-theta_I)./zeta_I).^2);
    
    U_En = (2-ce^2)/(1+ce)*U_E-(1-ce)/(1+ce)*U_Eo+ce^2/(1+ce)*F_E;
    U_In = (2-ci^2)/(1+ci)*U_I-(1-ci)/(1+ci)*U_Io+ci^2/(1+ci)*F_I;
    
    u_E(t+1) = (2-ce^2)/(1+ce)*u_E(t)-(1-ce)/(1+ce)*u_E(t-1)+ce^2/(1+ce)*f_E;
    u_I(t+1) = (2-ci^2)/(1+ci)*u_I(t)-(1-ci)/(1+ci)*u_I(t-1)+ci^2/(1+ci)*f_I;
    
    U_Eo = U_E;
    U_E = U_En;
    U_Io = U_I;
    U_I = U_In;
    
    SPo = reshape(F_E,N,N);
    rLFPo = reshape(-(J_EE+J_IE)*U_E,N,N);
    
    SP(:,:,t) = SPo(21:70,21:70);
    rLFP(:,:,t) = LFPo(21:70,21:70);
    
end