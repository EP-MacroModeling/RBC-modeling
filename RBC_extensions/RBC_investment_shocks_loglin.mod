% =========================================================================
% Baseline Log-Linear RBC - Macro II, PhD - Albonico Slides Notation
% =========================================================================

% -------------------------------------------------------------------------
% BRIEF CODE DESCRIPTION:
% This code presents a log-linearized version of a baseline RBC model.
% The structure of a log-linear Dynare file includes:
%
% 1) VARIABLES AND PARAMETERS
%    1.1) Declare endogenous and exogenous variables
%    1.2) Declare and calibrate parameters
%
% 2) STEADY STATE RELATIONSHIPS
%    No steady_state_model block needed for log-linear model, but we must
%    assign values to steady-state ratios used in the log-linear equations.
%
% 3) MODEL BLOCK
%    Use "model(linear)" to input log-linear equations, ending with "end;"
%
% 4) CHECKS
%    Use "steady;" and "check;" to validate the model's solution
%
% 5) SHOCKS
%    Define the stochastic structure using "shocks;" and "end;"
%
% 6) SIMULATION
%    Use "stoch_simul" to simulate impulse response functions (IRFs)
% -------------------------------------------------------------------------

clc;
%close all;

% =========================================================================
% 1) VARIABLES AND PARAMETERS
% =========================================================================

% Endogenous variables
var c n w rk k y a z inv; % z is the investment shock!

% Exogenous shock
varexo eta_a eta_z;

% Parameters
parameters Nss rho delta YKratio CKratio alpha CYratio IYratio rho_a rho_z;

% Parameter calibration
rho     = 0.005;    % Discount rate (beta = 1/(1+rho))
delta   = 0.025;    % Depreciation rate (10% annual)
alpha   = 0.33;     % Capital share
rho_a   = 0.9;      % AR(1) parameter for TFP
rho_z   = 0.9;      % AR(1) parameter for investment shock

% =========================================================================
% 2) STEADY STATE RELATIONSHIPS
% =========================================================================

% since we set Z=1, state conditions are unchanged 

YKratio = (delta + rho)/alpha;
IYratio = delta * alpha / (delta + rho);
CYratio = 1 - IYratio;
CKratio = CYratio * YKratio;
Nss     = 0.33;

% =========================================================================
% 3) MODEL
% =========================================================================

model(linear);

  // (1) Labor Supply
  w = c + Nss/(1 - Nss) * n;

  // (2) Euler Equation
  c(+1) = c + z + (rho + delta)/(1 + rho) * rk(+1) - (1 - delta)/(1 + rho) * z(+1);

  // (3) Capital Accumulation
  k = (1 - delta) * k(-1) + delta*z + YKratio * y - CKratio * c;

  // (4) Production Function
  y = a + alpha * k(-1) + (1 - alpha) * n;

  // (5) Labor Demand
  w = y - n;

  // (6) Capital Demand
  rk = y - k(-1);

  // (7) Market Clearing
  y = CYratio * c + IYratio * inv;

  // (8) TFP Shock Process
  a = rho_a * a(-1) + eta_a;

  // (9) Investment Shock Process
  z = rho_z * z(-1) + eta_z;

end;

% =========================================================================
% 4) CHECKS
% =========================================================================

steady;
check;

% =========================================================================
% 5) SHOCKS
% =========================================================================

shocks;
  var eta_a; stderr 0.009;
  var eta_z; stderr 0.01;
end;

% =========================================================================
% 6) SIMULATION
% =========================================================================

stoch_simul(irf=50) a z c n w rk k y a inv;