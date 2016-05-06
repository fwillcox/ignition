% ME 140 PROJECT 7: HYBRID ROCKET MOTOR
% Kristine Chen, Frankie Willcox, Eveylyn Xue, Andrew Lovett

% NOTE: All code for this project should be able to deal with effects of 
% both (1) chemical equilibrium and (2) frozen flow in nozzle assumptions

% CONSTANTS
PSI_PER_BAR = 14.5038;
At = NaN;                             % [m^2], Area of Nozzle Throat
Ae = NaN;                             % [m^2], Area of Nozzle Exit

npts = 100;
Patm =           1     * PSI_PER_BAR; % [psi], Assumption for Part 2
Pcombustor_actual = 68 * PSI_PER_BAR; % [psi], Assumption for Part 2 - Combustion Chamber Pressure
Pcombustor_lab = 10    * PSI_PER_BAR; % [psi], use for data taken in lab

% -------------------------------------------------------------------
% PART 1: Characteristic Velocity and Nozzle Temps vs. Mixture Ratio
% -------------------------------------------------------------------

% ASSUME:
% (i)   ideal polyethylene/oxygen rocket motor 
% (ii)  adiabatic combustion (dQ = 0)
% (ii)  variable Cp
% (iii) dissociation --TODO: how to account for this?

% Mixture Ratio (rmix)
rmix_lean = 10; 
rmix_rich = 1;
rmix = linspace(rmix_lean, rmix_rich,npts);

% Characteristic Velocity (cstar)
P0 = NaN;
mdot = NaN;
cstar = (P0 * At) / mdot;

% Nozzle Stagnation Temperature (T0)

% Nozzle Throat Temperature (T)


% -------------------------------------------------------------------
% PART 2: Nozzle Exit Velocity, Exit Temp, Thrust, Expansion Ratio
% -------------------------------------------------------------------
% ASSUME: 
% (i)   ideal nozzle
% (ii)  variable Cp
% (iii) dissociation

% Nozzle Exit Velocity (Ve)

% Nozzle Exit Temperature (Te)

% Thrust Coefficient (Cf)
% Note: F = Cf* mdot & cstar
Cf = F / (P0 * At);          

%  Optimal Nozzle Expanision Ratio (e)
e = Aexit / Athroat;

% ----------------------------------------
% PART 2: Improve the Hybrid Rocket Motor (using impulse/specific impulse)
% ----------------------------------------
% END GOAL: determine cstar & Cf for old & new design

