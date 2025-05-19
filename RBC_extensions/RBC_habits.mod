% =========================================================================
% Baseline NON-Linear RBC Habits - Esposito, Pollastri (2025)
% =========================================================================
% The GHH Utility specification:
%
%           u(c_t,n_t) = log(c_t - phi*c_t-1) + chi*log(1-n_t)
%
% Moreover, the budget constraint (not needed for habits):
%       c_t + inv_t + b_t+1 = w_t*n_t + r^k_t*k_t + (1+r_t-1)*b_t

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
var c n w rk r k y inv A B lambda; % (11 variables)

% Exogenous shock
varexo eta_a;

% Parameters
parameters phi chi lambdass Nss Kss Yss Iss Css Ass Bss rkss rss wss rho beta delta alpha rho_a;

% Parameter calibration
rho     = 0.005;    % Discount rate (beta = 1/(1+rho))
phi     = 0.8;
beta    = 1/(1+rho);
delta   = 0.025;    % Depreciation rate (10% annual)
alpha   = 0.33;     % Capital share
rho_a   = 0.9;      % AR(1) parameter for TFP

% =========================================================================
% 2) STEADY STATE RELATIONSHIPS
% =========================================================================


Nss = 0.33;
rkss = delta + rho;
rss = rho;
wss = (1-alpha)*((delta + rho)/alpha)^(alpha/(alpha -1));
lambdass = ((delta + rho)/alpha)^(alpha/1-alpha)*((delta + rho)/(rho + delta - alpha*delta))*1/(1-phi)*(1 + beta*phi);
chi = lambdass*wss*(1-Nss);
Kss = ((delta + rho)/alpha)^(1/(alpha-1))*Nss;
Yss = Kss^alpha*Nss^(1-alpha);
Iss = delta*Kss;
Css = Yss - Iss;
Ass = 1;
Bss = (1/rss)*(Css + Iss - wss*Nss - rkss*Kss);

% =========================================================================
% 3) MODEL
% =========================================================================

model;

  // (1) Labor Supply
  w = chi/(lambda*(1-n));

  // (2) FOC Consumption
  1/(c - phi*c(-1)) + beta*phi*(1/(c(+1)-phi*c)) = lambda;

  // (3) FOC Capitale
  lambda = beta*lambda(+1)*(rk(+1) + 1 - delta);

  // (4) Real interest rate on bonds
  r = rk(+1) - delta;

  // (5) Capital accumulation
  k = (1 - delta) * k(-1) + inv;

  // (6) Production Function
  y = A*k(-1)^alpha*n^(1-alpha);

  // (7) Labor Demand
  w = (1-alpha)*y/n;

  // (8) Capital Demand
  rk = alpha*y/k(-1);

  // (9) Market Clearing
  y = c + inv;

  // (10) Budget Constraint CHECK --> SHOULD I WRITE B AND B(-1)?
  c + inv + B(+1) = w*n + rk*k(-1) + (1+r(-1))*B;
% Chatgpt tells me that it is correct to write B and B(+1) but honestly I don't
% get why



  // (11) TFP Shock Process
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
r = rss;
B = Bss;
lambda = lambdass;

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

subplot(2,2,2)
plot(T,oo_.irfs.n_eta_a/oo_.mean(2,1),'LineWidth',2);
title('n  deviation')

subplot(2,2,3)
plot(T,oo_.irfs.w_eta_a/oo_.mean(3,1),'LineWidth',2);
title('w deviation')

subplot(2,2,4)
plot(T,oo_.irfs.y_eta_a/oo_.mean(7,1),'LineWidth',2);
title('y deviation')



% =========================================================================
% 8) STORE NORMALIZED IRFs - TFP Shock
%    
% =========================================================================

// Consumption
habit_c_eta_a = oo_.irfs.c_eta_a/oo_.mean(1,1);

// Labor
habit_n_eta_a = oo_.irfs.n_eta_a/oo_.mean(2,1); 

// Wages
habit_w_eta_a = oo_.irfs.w_eta_a/oo_.mean(3,1);
