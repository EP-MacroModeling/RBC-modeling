clc;
%close all;

% =========================================================================
% 1) VARIABLES AND PARAMETERS
% =========================================================================

% Endogenous variables
var c n w rk k a y inv;

% Exogenous variables
varexo eta_a;

% Parameters 
parameters Nss rho delta alpha rho_a YKratio CKratio CYratio IYratio;

% Parameter calibration
rho     = 0.005;    % Discount rate (beta = 1/(1+rho))
delta   = 0.025;    % Depreciation rate (10% annual)
alpha   = 0.33;     % Capital share
rho_a   = 0.9;    % AR(1) parameter for TFP

% =========================================================================
% 2) STEADY STATE RELATIONSHIPS
% =========================================================================
YKratio = (rho+delta)/alpha;
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
  c(+1) = c + (rho+delta)/(1-rho)*rk(+1);

// (3) Capital Accumulation
  k = (1 - delta) * k(-1) + YKratio * y - CKratio * c;

  // (4) Production Function
  y = a + alpha * k + (1 - alpha) * n;

  // (5) Labor Demand
  w = y - n;

  // (6) Capital Demand
  rk = y - k;

  // (7) Market Clearing
  y = CYratio * c + IYratio * inv;

  // (8) TFP Shock Process
  a = rho_a * a(-1) + eta_a;

end;

% =========================================================================
% 4) CHECKS
% =========================================================================
 
steady;
checks;

% =========================================================================
% 5) SHOCKS
% =========================================================================

shocks;
  var eta_a; stderr 0.009;
end;

% =========================================================================
% 6) SIMULATION
% =========================================================================

stoch_simul(irf=50) a c n w rk k y a inv;

