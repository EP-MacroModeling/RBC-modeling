% =========================================================================
% Baseline NON-Linear RBC GHH - Macro II, PhD - Albonico Slides Notation
% =========================================================================
% The GHH Utility specification:
%
%           u(c_t,n_t) = 1/1-sigma(c_t - n_t^(1+phi)/(1+phi))^(1-sigma) 
%
% --> removes the wealth/income effects (i.e. affect intratemporal 
% substitutability between consumption and labor). To see it properly, recall
% that from FOCs w_t = -Un_t/Uc_t. Here the utility is such that w_t = n_t^phi


% -------------------------------------------------------------------------
% BRIEF CODE DESCRIPTION:
% This code presents a NONLINEAR version of a baseline RBC model.
% The structure of the Dynare file includes:
%
% 1) VARIABLES AND PARAMETERS
%    1.1) Declare endogenous and exogenous variables
%    1.2) Declare and calibrate parameters (Notice: now also SS variables
%         must be included in the parameters)
%
% 2) STEADY STATE RELATIONSHIPS
%    We need to express all steady state variables as a function of 
%    parameters
%
% 3) MODEL BLOCK
%    Use "model;" to input nonlinear equations, ending with "end;"
%
% 4) INITIAL VALUES
%    With the command "initval;" assign initial values (i.e. steady state 
%    values) to all the variables in you model.
%
% 5) CHECKS
%    Use "steady;" and "check;" to validate the model's solution
%
% 5) SHOCKS
%    Define the stochastic structure using "shocks;" and "end;"
%
% 6) SIMULATION
%    Use "stoch_simul" to simulate impulse response functions (IRFs)
% -------------------------------------------------------------------------

clc;
close all;

% =========================================================================
% 1) VARIABLES AND PARAMETERS
% =========================================================================

% Endogenous variables
var c n w rk k y inv A;

% Exogenous shock
varexo eta_a;

% Parameters
parameters sigma phi Nss Kss Yss Iss Css Ass rho beta delta alpha rho_a;

% Parameter calibration
sigma = 1;          % Typicall ranges from 1 to 5 (2 is common choice)
rho     = 0.005;    % Discount rate (beta = 1/(1+rho))
beta    = 1/(1+rho);
delta   = 0.025;    % Depreciation rate (10% annual)
alpha   = 0.33;     % Capital share
rho_a   = 0.9;      % AR(1) parameter for TFP

% =========================================================================
% 2) STEADY STATE RELATIONSHIPS
% =========================================================================


Nss = 0.33;
Kss = ((delta + rho)/alpha)^(1/(alpha-1))*Nss;
Yss = Kss^alpha*Nss^(1-alpha);
Iss = delta*Kss;
Css = Yss - Iss;
Ass = 1;
rkss = delta + rho;
wss = (1-alpha)*Yss/Nss;
phi = log(wss - Nss);

% =========================================================================
% 3) MODEL
% =========================================================================

model;

  // (1) Labor Supply
  w = n^phi;

  // (2) Euler Equation
  (c(+1) - n(+1)^(1+phi)/(1+phi)) = (c - n^(1+phi)/(1+phi))*beta*(rk(+1) + 1 - delta)^(1/sigma);

  // (3) Capital Accumulation
  k = (1 - delta) * k(-1) + inv;

  // (4) Production Function
  y = A*k(-1)^alpha*n^(1-alpha);

  // (5) Labor Demand
  w = (1-alpha)*y/n;

  // (6) Capital Demand
  rk = alpha*y/k(-1);

  // (7) Market Clearing
  y = c + inv;

  // (8) TFP Shock Process
  log(A/Ass) = rho_a*log(A(-1)/Ass) + eta_a;

end;

% =========================================================================
% 4) INITIAL VALUES
% =========================================================================

initval;

c = Css;
inv = Iss;
y = Yss;
k = Kss;
n = Nss;
w = wss;
A = Ass;
rk = rkss;

end;



% =========================================================================
% 5) CHECKS
% =========================================================================

steady;
check;

% =========================================================================
% 5) SHOCKS
% =========================================================================

shocks;
  var eta_a; stderr 0.009;
end;

% =========================================================================
% 6) SIMULATION
% =========================================================================

% order = 1 will approximate the model to the first order around SS
% hp_filter = 1600 simulate time series with lambda = 1600
% ar = 4 reports autocorrelation of variables up to lag 4

%stoch_simul(order=1,irf=50,hp_filter=1600,ar=4) c n w rk k y inv;

stoch_simul(irf=50) c n w rk k y inv;

% Notice that when you run just stoch_simul(irf=50) dynare automatically gives
% order = 1. Moreover, hp_filter and ar do not affect IRFs, they just affect
% other stuff like simulated time series


% =========================================================================
% 7) PLOT IRFs AS PERCENT DEVIATIONS FROM STEADY STATE
%    
% =========================================================================

T=linspace(0,49,50);

figure
subplot(2,2,1)
plot(T,oo_.irfs.c_eta_a/oo_.mean(1,1),'LineWidth',2);
title('c deviation')

% Notice: mean(1,1) is because steady states are stored in a one column vector
% and the order is the one specified in "stoch_simul(irf=50) c n w rk k y inv;"
% such that consumption c is the first one (i.e. mean(1,1))

subplot(2,2,2)
plot(T,oo_.irfs.n_eta_a/oo_.mean(2,1),'LineWidth',2);
title('n  deviation')

subplot(2,2,3)
plot(T,oo_.irfs.w_eta_a/oo_.mean(3,1),'LineWidth',2);
title('w deviation')

subplot(2,2,4)
plot(T,oo_.irfs.y_eta_a/oo_.mean(6,1),'LineWidth',2);
title('y deviation')

% "y" is the sixth variable in the order so that to get the SS --> mean(6,1)

% =========================================================================
% 8) STORE NORMALIZED IRFs - TFP Shock
%    
% =========================================================================

// Consumption
ghh_c_eta_a = oo_.irfs.c_eta_a/oo_.mean(1,1);

// Labor
ghh_n_eta_a = oo_.irfs.n_eta_a/oo_.mean(2,1);

// Wages
ghh_w_eta_a = oo_.irfs.w_eta_a/oo_.mean(3,1);