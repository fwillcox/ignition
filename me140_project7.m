% ME 140 PROJECT 7: HYBRID ROCKET MOTOR
% Kristine Chen, Frankie Willcox, Eveylyn Xue, Andrew Lovett

% NOTE: All code for this project should be able to deal with effects of 
% both (1) chemical equilibrium and (2) frozen flow in nozzle assumptions

%kkdkd

% ------
% TODOS:
% ------
% (1) ask TAs about using lab data for mdot & cstar (see note below)
% (2) ask TAs about calcuating enthalpy of formaiton, hf from LHV & HHV
% (3) check general approach used to calculate T0 and T

clear all; close all; clc;

% ---------
% CONSTANTS
% ---------
PA_PER_BAR = 100000;
KG_PER_G = 10^-3;
C_TO_K = 273.15;
At = 0.0001824;                        % [m^2], Area of Nozzle Throat
Ae = NaN;                              % [m^2], Area of Nozzle Exit

npts = 10;             
P0_combustor_actual = 68 * PA_PER_BAR; % [Pa], Assumption for Part 2 - Combustion Chamber Pressure
P0_combustor_lab = 10    * PA_PER_BAR; % [Pa], use for data taken in lab
Patm =  101300;                        % [Pa]
Tref = 25 + C_TO_K;                    % [K]

% -------------------------------------------------------------------
% PART 1: Characteristic Velocity and Nozzle Temps vs. Mixture Ratio
% -------------------------------------------------------------------

% ASSUME:
% (i)   ideal polyethylene/oxygen rocket motor 
% (ii)  adiabatic combustion (dQ = 0)
% (ii)  variable Cp
% (iii) dissociation --TODO: how to account for this?

% Initialize the gas object -- Only do this once per script
gas = IdealGasMix('gri30.xml');

% Create the composition vectors
% There are 53 species includd in GRI30, so this will be the size of our array
nsp = nSpecies(gas); 

% Index of Species
iC2H4  = speciesIndex(gas,'C2H4');
iO2  = speciesIndex(gas,'O2'); 
iCO2 = speciesIndex(gas,'CO2');
iH2  = speciesIndex(gas,'H2');
iH2O = speciesIndex(gas,'H2O');
iOH  = speciesIndex(gas,'OH');

% Mixture Ratio (rmix)
rmix_lean = 10;                                     % rmix = mdot_O2 / mdot_C2H4_fuel
rmix_rich = 1;
rmix = linspace(rmix_lean, rmix_rich,npts);

            

for i = 1:length(rmix)
   % Mol Fractions (x)
   x_r = zeros(nsp,1);
   x_r(iC2H4) = 1;                                  % ASSUME: mol_fraction_fuel = 1
   x_r(iO2)   = x_r(iC2H4) * rmix(i);
   x_p = zeros(nsp,1);

   % Set state of gas (currently just products) to Tref and Patm
   set(gas,'T',Tref,'P',Patm,'X',x_r);

   % -------------------------------------------------------------------
   % ASK TA's: HOW DO WE USE LHV AND HHV TO GET ENTHALPY OF FORMATION???
   % -------------------------------------------------------------------

   % Find the Stagnation Enthalpy of the process h_r_C2H4
   h_r = enthalpy_mole(gas);                        % [J/kmol], Stagnation Enthalpy of C2H4, since we found it at Tref
   % hf_C2H4 = NaN;                                 % [J/kmol], get somehow from LHV and HHV
   % hdiff = hf_C2H4 - h_r;
   % h0 = h_r + h_diff;                             % [J/kmol], True Stagnation Enthalpy

   % Test temperatures until h_products = h_reactants 
   h_p = -1e10;                                     % [J/kmol], start at very low enthalpy
   Ttest = Tref;                                    
   dT = 1;                                          % [K], Temperature step

   % --------------------------------------------------
   % STEP 1: Find T0 (stagnation temperature at nozzle)
   % --------------------------------------------------
   % ASSUME: adiabatic, therefore h_products = h_reactants
   while h_p < h_r
        Ttest = Ttest + dT;
        set(gas,'T',Ttest,'P',P0,'X',x_p);
        h_p = enthalpy_mole(gas);                  
   end

   T0 = Ttest;                                      % [K] T0, Stagnation Temperature of Nozzle = Adiabatic Flame Temp of Combustor

    
   % -----------------------------------------------------------------
   % STEP 2 (equilibrate): Find Tthroat (throat temperature at nozzle)
   % -----------------------------------------------------------------
   assumption = 'allowedToEquilibrate'; % 'frozenAtNozzle' 
   switch 'allowedToEquilibrate'
       case 'allowedToEquilibrate'                
            h_test = 0;
            Ttest = T0;
            while htest < h0
                Ttest = T_test - dT;
                
                % Allow the gas to reach equilibrium at current Ttest
                set(gas,'T',Ttest,'P',P0,'X',x_p);        % Set gas state to Ttest & P0
                equilibrate(gas,'HP');                    % Allow reactants to reach equilibrium (minimize g) at constant h & P
                htest = enthalphy_mol(gas) + (soundSpeed(gas)^2)/2;
            end
            Tthroat = Ttest;
   
       case 'frozenAtNozzle'
            % Freeze the state of the gas to T0 & P0
            set(gas,'T',T0,'P',P0_combustor_lab,'X',x_r); % Set gas state to T0 & P0
            equilibrate(gas,'HP');                        % Allow reactants to reach equilibrium (minimize g) at constant h & P          
            
            h_test = 0;
            Ttest = T0;
            while htest < h0
                Ttest = T_test - dT;
                set(gas,'T',Ttest,'P',P0,'X',x_p);
                htest = enthalphy_mol(gas) + (soundSpeed(gas)^2)/2;
            end
            Tthroat = Ttest;
       
       otherwise 
           warning('Unexpected assumption. Set assumption = to allowedToEquilibrate or frozenAtNozzle');
   end

% -------------------------------------------------------------------
% ASK TA's: DO WE USE LAB DATA TO CALCULATE MDOT_TOTAL SO THAT WE CAN GET CSTAR???
% -------------------------------------------------------------------
% Characteristic Velocity (cstar)
% dm = 10.7 * KG_PER_G;         % [kg] total change in mass
% tfire = NaN;                  % [s] firing period, from our lab data
% mdot_total = dm / tfire;      % [kg/s] mass flow rate of fuel & oxidizer combined
cstar = (P0 * At) / mdot_total;



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

%comment for testing 

% ----------------------------------------
% PART 3: Improve the Hybrid Rocket Motor (using impulse/specific impulse)
% ----------------------------------------
% END GOAL: determine cstar & Cf for old & new design

