%% Description
% This function simulates the Brion dynamic concentration for batch BZ reaction based on the reacant initial
% concentration (Malonic acid,Bromate ion and Cerium(III) ion). BZ model is based on the paper "Malonic acid 
% concentration as a control parameter in the kinetic analysis of the Belousov–Zhabotinsky reaction under 
% batch conditions"(DOI:https://doi.org/10.1039/B804919J)

%%
function BrionConcList = InSilicoBZReaction(ReactantConcList)

%% BZ reaction parameters 

% Input final time and number of times step
tProbeEnd = 4000; % track reaction till 4000 second
tProbeTotalStep = 4000;% set number of time step %4000

%% Perform Batch BZ reaction simulation

BrionConcList = zeros(size(ReactantConcList,1),tProbeTotalStep);

for i = 1:size(ReactantConcList,1)
    BrionConcList(i,:) = BatchBZReaction_SlavicaFun(tProbeEnd,tProbeTotalStep,ReactantConcList(i,:));
end

%% Batch BZ model (Slavica) description

% The batch BZ reaction in this code is based on the paper "Malonic acid concentration as a control parameter 
% in the kinetic analysis of the Belousov–Zhabotinsky reaction under batch conditions"
% (DOI:https://doi.org/10.1039/B804919J). In this code, it is possible to alter 3 of the initial reactant
% concentration (Malonic acid (CH2COOH),Bromate ion (BrO3-), Cerium(III) ion (Ce3+)).The output of the function 
% is the time series concentration of the Bromide ion (Br-). 

% Chemical notation in code: 
% Malonic acid - MA
% Bromide ion - Brion
% Bromate ion - BrO3ion
% Cerium(III) ion - Ce3ion
% Cerium(IV) ion - Ce4ion
% Hypobromous acid - HOBr
% Hydrogen ion - Hion
% Bromine - Br2
% Bromous acid - HBrO2
% Dibromine monoxide - Br2O
% Bromine dioxide (radical) - BrO2
% Bromomalonic acid - BrMA

function BrionConcTimeSeries = BatchBZReaction_SlavicaFun(tProbeEnd,tProbeTotalStep,ReactantConc)

%% Set parameters

% Rate constant
k1  = 2.55e9;  % (mol^-2)(dm^6)(s^-1)
k_1 = 3.18;    % s^-1
k2  = 5.93e6;  % (mol^-2)(dm^6)(s^-1)
k3  = 3.21e3;  % s^-1
k_3 = 3.22e8;  % (mol^-1)(dm^3)(s^-1)
k4  = 2.86;    % (mol^-3)(dm^9)(s^-1)
k5  = 3.49e3;  % (mol^-1)(dm^3)(s^-1)
k6  = 44.70;   % (mol^-2)(dm^6)(s^-1)
k_6 = 6.7e7;   % (mol^-1)(dm^3)(s^-1)
k7  = 3.2e4;   % (mol^-2)(dm^6)(s^-1)
k_7 = 1.12e4;  % (mol^-1)(dm^3)(s^-1)
k8  = 4.24;    % (mol^-1)(dm^3)(s^-1)
k9  = 0.36;    % (mol^-1)(dm^3)(s^-1)
k10 = 47.17;   % (mol^-1)(dm^3)(s^-1)
k11 = 4.23e-2; % (mol^-1)(dm^3)(s^-1)
k12 = 1.10e-2; % s^-1

% Initial value for products and reactants(that are not part of the changeable variables) - values are based
% on values within the paper 

Hion0   = 1.29;   % (mol)(dm^-3)
HOBr0   = 1.5e-8; % (mol)(dm^-3)
Br20    = 0;      % (mol)(dm^-3)
HBrO20  = 0;      % (mol)(dm^-3)
Br2O0   = 0;      % (mol)(dm^-3)
BrO20   = 0;      % (mol)(dm^-3)
Ce4ion0 = 0;      % (mol)(dm^-3)
BrMA0   = 0;      % (mol)(dm^-3)
Brion0  = 1.5e-5; % (mol)(dm^-3)

%% Set up ode equation

% initial concentration 
MA0      = ReactantConc(1); % mol/dm^3
BrO3ion0 = ReactantConc(2); % mol/dm^3
Ce3ion0  = ReactantConc(3); % mol/dm^3

% Probing time frame
t0Probe = linspace(0,tProbeEnd,tProbeTotalStep);

% Set up differential equation

dMA_dt      = @(MA,Brion,BrO3ion,Ce3ion,Ce4ion,HOBr,Hion,Br2,HBrO2,Br2O,BrO2,BrMA) ...
              (-k8*MA*Br2 - k9*MA*Ce4ion - k11*Br2O*MA);

dBrion_dt   = @(MA,Brion,BrO3ion,Ce3ion,Ce4ion,HOBr,Hion,Br2,HBrO2,Br2O,BrO2,BrMA) ...
              (k_1*Br2 + k8*MA*Br2 + k10*BrMA*Ce4ion - k1*Brion*HOBr*Hion - k2*HBrO2*Brion*Hion - k4*Brion*BrO3ion*(Hion^2));

dBrO3ion_dt = @(MA,Brion,BrO3ion,Ce3ion,Ce4ion,HOBr,Hion,Br2,HBrO2,Br2O,BrO2,BrMA) ...
              (k5*(HBrO2^2) + k_6*(BrO2^2) - k4*Brion*BrO3ion*(Hion^2) - k6*BrO3ion*HBrO2*Hion);

dCe3ion_dt  = @(MA,Brion,BrO3ion,Ce3ion,Ce4ion,HOBr,Hion,Br2,HBrO2,Br2O,BrO2,BrMA) ...
              (k_7*Ce4ion*HBrO2 + k9*MA*Ce4ion + k10*BrMA*Ce4ion-k7*Ce3ion*BrO2*Hion);

dCe4ion_dt  = @(MA,Brion,BrO3ion,Ce3ion,Ce4ion,HOBr,Hion,Br2,HBrO2,Br2O,BrO2,BrMA) ...
              (k7*Ce3ion*BrO2*Hion - k_7*Ce4ion*HBrO2 - k9*MA*Ce4ion - k10*BrMA*Ce4ion);

dHOBr_dt    = @(MA,Brion,BrO3ion,Ce3ion,Ce4ion,HOBr,Hion,Br2,HBrO2,Br2O,BrO2,BrMA) ...
              (k_1*Br2 + 2*k3*Br2O + k4*Brion*BrO3ion*(Hion^2) + k5*(HBrO2^2) + k11*Br2O*MA - k1*Brion*HOBr*Hion - 2*k_3*(HOBr^2));

dHion_dt    = @(MA,Brion,BrO3ion,Ce3ion,Ce4ion,HOBr,Hion,Br2,HBrO2,Br2O,BrO2,BrMA) ...
              (k_1*Br2 + k5*(HBrO2^2) + k_6*(BrO2^2) + k_7*Ce4ion*HBrO2 + k8*MA*Br2 + k9*MA*Ce4ion - k1*Brion*HOBr*Hion - k2*HBrO2*Brion*Hion - 2*k4*Brion*BrO3ion*(Hion^2) - k6*BrO3ion*HBrO2*Hion - k7*Ce3ion*BrO2*Hion);

dBr2_dt     = @(MA,Brion,BrO3ion,Ce3ion,Ce4ion,HOBr,Hion,Br2,HBrO2,Br2O,BrO2,BrMA) ...
              (k1*Brion*HOBr*Hion - k_1*Br2 - k8*MA*Br2 - k12*Br2); 

dHBrO2_dt   = @(MA,Brion,BrO3ion,Ce3ion,Ce4ion,HOBr,Hion,Br2,HBrO2,Br2O,BrO2,BrMA) ...
               (k4*Brion*BrO3ion*(Hion^2) + k_6*(BrO2^2) + k7*Ce3ion*BrO2*Hion - k2*HBrO2*Brion*Hion - 2*k5*(HBrO2^2) - k6*BrO3ion*HBrO2*Hion - k_7*Ce4ion*HBrO2);

dBr2O_dt    = @(MA,Brion,BrO3ion,Ce3ion,Ce4ion,HOBr,Hion,Br2,HBrO2,Br2O,BrO2,BrMA) ...
              (k2*HBrO2*Brion*Hion + k_3*(HOBr^2) - k3*Br2O - k11*Br2O*MA);

dBrO2_dt    = @(MA,Brion,BrO3ion,Ce3ion,Ce4ion,HOBr,Hion,Br2,HBrO2,Br2O,BrO2,BrMA) ...
              (2*k6*BrO3ion*HBrO2*Hion + k_7*Ce4ion*HBrO2 - 2*k_6*(BrO2^2) - k7*Ce3ion*BrO2*Hion);

dBrMA_dt    = @(MA,Brion,BrO3ion,Ce3ion,Ce4ion,HOBr,Hion,Br2,HBrO2,Br2O,BrO2,BrMA) ...
              (k8*MA*Br2 + k11*Br2O*MA - k10*BrMA*Ce4ion);

% Rescale initial value for each variables
InitialValues = [MA0,Brion0,BrO3ion0,Ce3ion0,Ce4ion0,HOBr0,Hion0,Br20,HBrO20,Br2O0,BrO20,BrMA0];


%% Perform ODE calculation

F = ode(ODEFcn=@(t,yt) [dMA_dt(yt(1),yt(2),yt(3),yt(4),yt(5),yt(6),yt(7),yt(8),yt(9),yt(10),yt(11),yt(12)); ...
                        dBrion_dt(yt(1),yt(2),yt(3),yt(4),yt(5),yt(6),yt(7),yt(8),yt(9),yt(10),yt(11),yt(12)); ...
                        dBrO3ion_dt(yt(1),yt(2),yt(3),yt(4),yt(5),yt(6),yt(7),yt(8),yt(9),yt(10),yt(11),yt(12)); ...
                        dCe3ion_dt(yt(1),yt(2),yt(3),yt(4),yt(5),yt(6),yt(7),yt(8),yt(9),yt(10),yt(11),yt(12)); ...
                        dCe4ion_dt(yt(1),yt(2),yt(3),yt(4),yt(5),yt(6),yt(7),yt(8),yt(9),yt(10),yt(11),yt(12)); ...
                        dHOBr_dt(yt(1),yt(2),yt(3),yt(4),yt(5),yt(6),yt(7),yt(8),yt(9),yt(10),yt(11),yt(12)); ...
                        dHion_dt(yt(1),yt(2),yt(3),yt(4),yt(5),yt(6),yt(7),yt(8),yt(9),yt(10),yt(11),yt(12)); ...
                        dBr2_dt(yt(1),yt(2),yt(3),yt(4),yt(5),yt(6),yt(7),yt(8),yt(9),yt(10),yt(11),yt(12)); ...
                        dHBrO2_dt(yt(1),yt(2),yt(3),yt(4),yt(5),yt(6),yt(7),yt(8),yt(9),yt(10),yt(11),yt(12)); ...
                        dBr2O_dt(yt(1),yt(2),yt(3),yt(4),yt(5),yt(6),yt(7),yt(8),yt(9),yt(10),yt(11),yt(12)); ...
                        dBrO2_dt(yt(1),yt(2),yt(3),yt(4),yt(5),yt(6),yt(7),yt(8),yt(9),yt(10),yt(11),yt(12)); ...
                        dBrMA_dt(yt(1),yt(2),yt(3),yt(4),yt(5),yt(6),yt(7),yt(8),yt(9),yt(10),yt(11),yt(12))], ...
        InitialValue=InitialValues,...
        Solver = "ode15s",AbsoluteTolerance=1e-9,RelativeTolerance=1e-9); 

S = solve(F,t0Probe);

%% OutputValues

BrionConcTimeSeries = log10(S.Solution(2,:)); % Brion concentration dynamic concentration

end

end