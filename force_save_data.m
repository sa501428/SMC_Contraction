%
% Generate force using data from ACH files. Plot for vs time and save data
% as txt files. 
%

function force_save_data
global f1 iter f r k Cach Tach vCa1 vCa2 RCa cNaCa vNaCa cc vd Rd vCl vK vCa3 RK cb sc lambda cW...
    beta gamma GCa GNaCa B C D L FNaK GCl GK F Kr epsilon VM4 k4 PMV kV r2 uu Fcontract myT Calcium
global maxAch

% - Running ODEs for Muscle Contraction

% variables for ode equations part A
beta = .13;         % translation factor
gamma = 197;        % scaling factor
epsilon = .9;       % rate constant for linear IP3
lambda = 45;        % channel constant
B = 2.025;          % SR uptake rate constant
C = 55;             % CICR rate constant
D = .24;            % Ca extrusion by ATPase constant
F = .23;            % maximal influx rate
L = .025;           % leak from SR rate constant
cb = 1;             % half point SR ATPase activation
cc = .9;            % half point CICR activation
cNaCa = .5;         % half point Na Ca exchange activation
cW = 0;             % translation factor
FNaK = .0432;       % net whole cell flux
GCa = .00129;       % whole cell conductance for VOCCs
GCl = .00134;       % whole cell conductance Cl
GK = .00446;        % whole cell conductance K
GNaCa = .00316;     % whole cell conductance for Na Ca exchange
k4 = .5;            % half saturation constant IP3 degradation
Kr = 1;             % half saturation constant Ca entry
kV = -58;           % half saturation constant IP3 voltage synthesis
PMV = .8;           % max rate voltage IP3 synthesis
r2 = 8;             % hill coefficient
RCa = 8.5;          % maximum slope of VOCC activation
Rd = 250;           % slope of voltage dependance
RK = 12;            % maximum slope Ca activation
sc = 2;             % half point CICR efflux
uu = 4;              % hill coefficient
vCa1 = 100;         % reversal potential VOCCs
vCa2 = -24;         % half point VOCC activation
vCa3 = -27;         % half point Ca channel activation
vCl = -25;          % reversal potential Cl
vd = -100;          % intercept voltage dependance
vK = -104;          % reversal potential K
VM4 = 2;            % max nonlinear IP degradation
vNaCa = -40;        % reversal potential Na Ca exchange
Fcontract = [];
myT = [];
Calcium = [];

% forward rates, reverse rates, k values for ODE part B
f = [12 480 5 840 28 120 7.5 5 7.6];            % forward rates
r = [12 1200 135 45.4 .0308 4 3.75 25 22.8];    % reverse rates
k = [27 10 15 5 16 15 10];                      % catalyst rates

% initial levels of materials
%    [ACH   IP3 CSR Wi Vi   Ca   CaM  MLCK   C2CM C4CM CM_MK C2CM_MK C4CM_MK CM_BP BP      M       Mp   AMp    AM MLCP]
x20 = [10^-3 .36 1.3 .1 -38 .19 .9285 9.6506 .0015 0.0 0.3332 0.2713 0.013 2.8207 15.1793 23.9558 .0144 .0166 .0132 7.5];
x20 = x20';


% initial state vector
q0 = x20;


% load ACH data from files
Tach = dlmread('timeACH.txt');
Tach = Tach * 10^6;
Cach1 = dlmread('conACHR2b.txt');
figure(1)
linecol = ['b' 'g' 'r' 'c' 'm'];

maxAch = dlmread('maxACH.txt');

% run ODE system for each synapse distance
for mycolumn = 1:5
    Cach = Cach1(:,mycolumn);
    [t,q] = ode23(@(t,x) basicspringsrhs(t,x,mycolumn),[0,20],q0);
    
    % save data
    dlmwrite(strcat(num2str(mycolumn),'force_timeR2b.txt'), myT)
    dlmwrite(strcat(num2str(mycolumn),'force_forceR2b.txt'), Fcontract)
    dlmwrite(strcat(num2str(mycolumn),'force_CaR2b.txt'), Calcium)   
    figure(1)
    ids = find(abs(myT)>10,1);
    %Calcium(ids) = 0;
    %Fcontract(ids) = 0;
    plot(myT(1:ids),Calcium(1:ids),linecol(mycolumn))
    hold on
    
    % reset the vector which is saving the data
    myT = [];
    Fcontract = [];
    Calcium = [];
end
% label graph
legend('20 microns','30 microns','40 microns','50 microns','60 microns')
xlabel('Time (seconds)')
ylabel('Ca^{2+} Concentration (\mu M)')
hold off
return

%
% function basicspringsrhs
%
% sub-function which gets position and velocities for nodes based on the
% force equations
%
% Inputs:
%   q       matrix containing positions and velocities
%
% usage: type  basicspringsrhs in ode
%
% example:  ode23(@basicspringsrhs,[T,T+tinc],q0);
%

function dx = basicspringsrhs(t,x,j)

global f r k vCa1 vCa2 RCa cNaCa vNaCa cc vd Rd vCl vK vCa3 RK cb sc lambda cW...
    beta gamma GCa GNaCa B C D L FNaK GCl GK F Kr epsilon VM4 k4 PMV kV r2 uu Fcontract myT Calcium
global maxAch

%location of x variables, y variables, ..


% extract out node positions/velocities

%t
%x2(1)
%x2(56)

% Equations
% Acetylcholine
dx(1,1) = d_ach(t);

% Inositol 1,4,5-trisPhosphate recpetor
dx(2,1) = maxAch(j) - epsilon*x(2) - VM4 * (x(2)^uu/( x(2)^uu + k4^uu )) + PMV*(1 - x(5)^r2/( kV^r2 + x(5)^r2 ));   

% Sarcoplasmic Calcium
dx(3,1) = (B*x(6)^2/(x(6)^2 + cb^2)) ...
    - (C*(x(3)^2)*(x(6)^4)/(((x(3)^2) + ...
    (sc^2))*((x(6)^4) + (cc^4)))) ...
    - (L*x(3));

% K channel
dx(4,1) = lambda*( (( x(6) + cW )^2 /( ( x(6) + cW )^2 + beta*exp(-(x(5)-vCa3)/RK) )) ...
    - x(4) );

% Cell Potential
dx(5,1) = gamma*( - (FNaK) ...
    - (GCl*(x(5) - vCl)) ...
    - 2*(GCa * (x(5) - vCa1)/(1 + exp(-(x(5) - vCa2)/RCa))) ...
    - (GNaCa*x(6)*(x(5) - vNaCa)/(x(6) + cNaCa)) ...
    - (GK*x(4)*(x(5)-vK)) );

% Cystolic Calcium
dx(6,1) = (F*x(2)^2/(x(2)^2 + (Kr)^2)) ...
    - (GCa * (x(5) - vCa1)/(1 + exp(-(x(5) - vCa2)/RCa))) ...
    + (GNaCa*x(6)*(x(5) - vNaCa)/(x(6) + cNaCa)) ...
    - (B*x(6)^2/(x(6)^2 + cb^2)) ...
    + (C*(x(3)^2)*(x(6)^4)/(((x(3)^2) + (sc^2))*((x(6)^4) + (cc^4)))) ...
    - (D*x(6)*(1 + (x(5) - vd)/Rd)) ...
    + (L*x(3));

% Part B: CaM MLCK Ca2CaM Ca4CaM CaM_MLCK Ca2CaM_MLCK Ca4CaM_MLCK CaM_BP BP M Mp AMp AM MLCP
dx(7,1)  = r(1)*x(9) + r(3)*x(11) + r(8)*x(14) - f(1)*x(6)*x(6)*x(7) - f(3)*x(7)*x(8) - f(8)*x(7)*x(15);
dx(8,1)  = r(3)*x(11) + r(4)*x(12) + r(5)*x(13) - f(3)*x(7)*x(8) - f(4)*x(8)*x(9) - f(5)*x(8)*x(10);
dx(9,1)  = f(1)*x(6)*x(6)*x(7) + r(2)*x(10) + r(4)*x(12) + f(9)*x(6)*x(6)*x(14) - r(1)*x(9) - f(2)*x(6)*x(6)*x(9) - f(4)*x(8)*x(9) - r(9)*x(9)*x(15);
dx(10,1)  = f(2)*x(6)*x(6)*x(9) + r(5)*x(13) - r(2)*x(10) - f(5)*x(8)*x(10);
dx(11,1)  = f(3)*x(7)*x(8) + r(6)*x(12) - r(3)*x(11) - f(6)*x(6)*x(6)*x(11);
dx(12,1)  = f(4)*x(8)*x(9) + f(6)*x(6)*x(6)*x(11) + r(7)*x(13) - r(4)*x(12) - r(6)*x(12) - f(7)*x(6)*x(6)*x(12);
dx(13,1)  = f(5)*x(8)*x(10) + f(7)*x(6)*x(6)*x(12) - r(5)*x(13) - r(7)*x(13);
dx(14,1)  = f(8)*x(7)*x(15) + r(9)*x(9)*x(15) - r(8)*x(14) - f(9)*x(6)*x(6)*x(14);
dx(15,1) = r(8)*x(14) + f(9)*x(6)*x(6)*x(14) - f(8)*x(7)*x(15) - r(9)*x(9)*x(15);
dx(16,1) = -k(1)*x(13)*x(16)/(k(2) + x(16)) + k(5)*x(20)*x(17)/(k(6) + x(17)) + k(7)*x(19);
dx(17,1) = k(1)*x(13)*x(16)/(k(2) + x(16)) - k(5)*x(20)*x(17)/(k(6) + x(17)) -k(3)*x(17) + k(4)*x(18);
dx(18,1) = k(3)*x(17) - k(4)*x(18) + k(1)*x(13)*x(19)/(k(2) + x(19)) - k(5)*x(20)*x(18)/(k(6) + x(18));
dx(19,1) = -k(1)*x(13)*x(19)/(k(2) + x(19)) + k(5)*x(20)*x(18)/(k(6) + x(18)) - k(7)*x(19);
dx(20,1) = 0;

% calculate and save force and time
Fcontract = [Fcontract; (x(18) + x(19))];
myT = [myT; t];
Calcium = [Calcium; x(6)];

return



%% function d_ach
%
% sub-function which provides derivative of [ach] vs time function
%
% Inputs:
%   t       instance in time
%
% Output:
%   da      derivative at instance in time
%
% usage: type  d_ach(t)
%

function da = d_ach(t)
global Cach Tach
t = mod(t,20);
if t < 19.99925
    time = round(t*4000) + 2;
    da = (Cach(time-1) - Cach(time + 1)) / (Tach(time-1) - Tach(time + 1));
else
    da = 0;
end
return
