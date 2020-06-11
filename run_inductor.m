%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             Author: Paavo Rasilo & Wilmar Martinez                  %%%
%%%             Title:  Run Inductor v1.1                               %%%
%%%             Date:   24.05.2020                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================
% Runs the Simulink model for a PWM-supplied lamination core inductor
% taking into account hysteresis, eddy-current and excess losses in the
% laminated core.
%
% Note: Accounting for the excess losses (cex > 0) causes an algebraic loop
% in the Simulink model. The warning related to this algebraic loop can be
% disabled in the Simulink model parameters.
%
% Note2: It is set by default to evalaute 400Hz of fundamental frequency.
%
% Note3: For more info please follow the recent publications of the authors.
%==========================================================================

clear all; close all; clc
addpath hyst util

mu0 = 4*pi*1e-7;
nu0 = 1/mu0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DC link voltage, modulation index, phase of the duty cycle and initial
% value of the flux density. These need to be iterated to obtain the
% desired flux density.
Udc = 128.11;   
am = 0.5;
uphase = 0;
B0 = -1;
%uphase = -0.0199;
%B0 = -1;

% Fundamental and switching frequencies and deadtime
f = 400;
fs = 500e3;
deadtime = 100e-9;

% Other parameters
nb       = 1; % Number of skin-effect terms
flag_hy  = 1; % Account for hysteresis or not?

% Total number of fundamental periods to be simulated
periods = 1.1;

% Dimensions of the toroid
h = 7e-3;           % Height of iron part
do = 127e-3;        % Outer diameter
di = 102e-3;        % Inner diameter
Afe = (do-di)/2*h;  % Iron cross-section area
lfe = (do+di)/2*pi; % Length of flux path
N1 = 254;           % Number of primary turns
delta = 0;          % Air-gap length
d = 0.35e-3;        % Thickness of core lamination
sigma = 1/52e-8;    % Electrical conductivity of core
dens = 7650;        % Mass density of core
cex = 0.3137;       % Excses loss coefficient
R = 0.56;           % Resistance
Lsig = 1.8735e-06;  % Leakage inductance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate some values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulink file name
if nb == 1;     sname = 'inductor'; else
    sname = 'inductor_skineffect'; end

% Skin effect related terms
if nb > 1
    
    % Basis function values in integration points [2]
    nz = nb;
    [z,w] = lgwt(nz, -0.5, 0);
    w = w/sum(w);
    clear alfa
    for i = 1 : nb
        alfa(:,i) = cos((i-1)*z*2*pi);
    end
    
    % Dynamic coefficient matrix
    C = zeros(nb,nb);
    for m = 0 : nb-1
        for n = 0 : nb-1
            if (m == 0)  && (n == 0)
                cc = 1/12;
            elseif (m == n) && (n > 0)
                cc = 1/(2*pi^2*(m+n)^2);
            elseif (m*n == 0) && (m+n>0)
                cc = (-1)^(m+n+1)/(4*pi^2*(m+n)^2);
            else
                cc = 0;
            end
            C(m+1,n+1) = cc;
        end
    end
    C = sigma*d^2*C;
    Cinv = inv(C(2:end,2:end));
end

% Eddy-current loss resistance
if nb == 1
    Rfe = N1^2*Afe/lfe/(sigma*d^2/12);
else
    Rfe = N1^2*Afe/lfe/(sigma*d^2/12-C(1,2:end)*Cinv*C(2:end,1));
end

% Carrier frequency
fcar = fs/2;
fprintf('Simulating with hy = %d, nb = %d, switching %g kHz, deadtime %g ns\n', flag_hy, nb, fs/1e3, deadtime*1e9);

% Open Simulink file and scope for plotting
if ~bdIsLoaded(sname)
    open_system(sname);
end
open_system([sname '/Scope for everything']);

% If excess-loss coefficient is set to zero comment out the excess loss
% branch to avoid algebraic loop
blockname = [sname '/Laminated-core inductor with hysteresis, eddy currents and excess losses/Excess-loss model'];
if cex <= 0
    set_param(blockname, 'commented', 'on');
else
    set_param(blockname, 'commented', 'off');
end

% Simulation time
if deadtime == 0
    dt = 150e-9;
else
    dt = deadtime/2;
end
nper = round(1/f/dt);
ntot = periods*nper;
t = 0:dt:ntot*dt;
tmax = t(end);

% Simulate
warning off
sim(sname);
warning on

% Simulation stops before finishing
if max(tout) ~= tmax
    fprintf('Stopped by user?\n')
    return
    
    % Finished correctly
else
    
    % Take only last period
    % U and I have different sampling than B and H since they are
    % measured from Simscape with a discrete solver
    ind = (tout >= tmax-1/f);
    tout = tout(ind);
    U = interp1(U.time, U.signals.values, tout);
    I = interp1(I.time, I.signals.values, tout);
    B = B(ind);
    H = H(ind);
    phy = phy(ind);
    pcl = pcl(ind);
    pex = pex(ind);
    tout = tout-tout(1);
    
 % Centering
 B = B-(max(B)+min(B))/2;
 H = H-(max(H)+min(H))/2; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Handle results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tot(1,1)=Udc;
Tot(2,1)=max(B);
Tot(3,1)=max(H);
Tot(4,1)=phy(end);
Tot(5,1)=pcl(end);
Tot(6,1)=pex(end);
Tot(7,1)=ploop(end);

% Show Losses
fprintf('  Maximum flux density:      %g T\n', max(B));
fprintf('  Maximum flux intensity:    %g A/m\n', max(H));
fprintf('  Hysteresis loss:           %g W/kg\n', phy(end));
fprintf('  Eddy-current loss:         %g W/kg\n', pcl(end));
fprintf('  Excess loss:               %g W/kg\n', pex(end));
fprintf('  Sum:                       %g W/kg\n', phy(end)+pcl(end));
fprintf('  Total loss from B(H) loop: %g W/kg\n', ploop(end));

% Plot B(H) loop
figure;
plot(H, B);
title(sprintf('{\\itf}_{s} = %g kHz', fs/1e3),'FontSize',16,'FontName','Times New Roman');
xlabel('Magnetic field intensity {\itH}_s (A/m)','FontSize',16,'FontName','Times New Roman');
ylabel('Magnetic flux density {\itB}_0 (T)','FontSize',16,'FontName','Times New Roman');
xlim(1.1*[min(H) max(H)]);
ylim(1.1*[min(B) max(B)]);
grid on

ddt=0:(1/f*1000)/length(I):(1/f*1000);