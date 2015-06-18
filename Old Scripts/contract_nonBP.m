% 
% cellbase2multistroke - simple cell setup
%
% two different k values (kheavy for the actin filaments, and klight for
% the myosin, crossbridges and cell membrane)
%
% multiple powerstrokes based on time (the nodes shift by a value of
% strokeinc every 3 time steps); 
%

function contract_nonBP


% global parameters
global d K M A L0 p tops bottoms centers Nmem Fx0 Nact CA N f r k CaMax Fr ff

CaMax = .72;
%tol = .01;  % tolerance for velocity when trying to reach a motionless cell
k = 1;      % spring constant
d = .9; % drag coefficient for nodes
M = 0.1; % mass of each node
ff = 1;  %frame counter for movie
Fr = zeros(1,81*3);
Tmax = 20;
p = 12;

kheavy = 150;
pins = [];

% Ellipse for membrane
Nmem = 40;  % total number of nodes along membrane
rmaj = 100;   % radius for long axis
rmin = 5;   % radius for short axis
theta = (0:2*pi/Nmem:2*pi-(2*pi)/Nmem)';    % vector of angles
xmem = rmaj*cos(theta);
ymem = rmin*sin(theta);


% Lines for upper and lower actin filaments
Nact = 15;  % number of nodes on each actin filament
ya1 = ymem(Nmem/2); xa1 = xmem(Nmem/2);     % where top actin starts
ya2 = ymem(Nmem); xa2 = xmem(Nmem);     % where bottom actin starts
xact1 = (xa1-xa1/Nact:-xa1/Nact:0)';   yact1 = zeros(length(xact1),1);
xact2 = (xa2-xa2/Nact:-xa2/Nact:0)';  yact2 = zeros(length(xact2),1);

% Line for central myosin filament
Nmy = 7; % Number of nodes on myosin; Number of crossbridges + 2
bridgespace = xact1(end)-xact1(end-1);   % space between crossbridges
xmy = transpose(-((Nmy-1)/2)*bridgespace:bridgespace:((Nmy-1)/2)*bridgespace);
ymy = zeros(length(xmy),1);


x0 = [xmem; xact1; xact2; xmy]; % concatenate x and y vectors for all cell components
y0 = [ymem; yact1; yact2; ymy];

N = length(x0); % number of nodes (N = Nmem + 2*Nact + Nmy)
u0 = zeros(N,1); % x-component of velocity of N nodes
v0 = zeros(N,1); % y-component of velocity of N nodes

% Assemble blocks within A matrix
% Membrane block
Imem = eye(Nmem);
Amem = [Imem(2:end,:); Imem(1,:)] + [Imem(:,2:end) Imem(:,1)];
% Actin block
Iact = eye(Nact-1);
Aact = [zeros(1,Nact); Iact zeros(Nact-1,1)];
Aact = Aact + Aact';
% Myosin block
Imy = eye(Nmy-1);
Amy = [zeros(1,Nmy); Imy zeros(Nmy-1,1)];
Amy = Amy + Amy';
%

% Assemble whole A matrix
A = zeros(N,N);
A(1:Nmem,1:Nmem) = Amem;    % First block membrane
A(Nmem+1:Nact+Nmem,Nmem+1:Nact+Nmem)= Aact; % Second block upper actin
A(Nmem+Nact+1:Nmem+2*Nact,Nmem+Nact+1:Nmem+2*Nact) = Aact;  % Second block upper actin
A(Nmem+2*Nact+1:end,Nmem+2*Nact+1:end) = Amy;   % Last block myosin
A(Nmem+1,Nmem/2+1) = 1;       % upper actin to membrane
A(Nmem/2+1,Nmem+1) = 1;       % "
A(Nmem+Nact+1,1) = 1;    % lower actin to membrane
A(1,Nmem+Nact+1) = 1;    % "

% Add crossbridges to A matrix
i = Nmy-3;
for j = 1:Nmy-2
    A(Nmem+2*Nact+1+j,Nmem+Nact-i) = 1;   % draw upper crossbridges
    A(Nmem+Nact-i,Nmem+2*Nact+1+j) = 1;
    A(Nmem+2*Nact+1+j,Nmem+2*Nact-j+1) = 1;   % draw lower crossbridges
    A(Nmem+2*Nact-j+1,Nmem+2*Nact+1+j) = 1;
    i = i-1;
end


% CYTOSKELETON - around edges
%mems = [1:Nmem 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
%for i = 1:Nmem
%    A(mems(i),mems(i+4)) = 1;
%    A(mems(i+4),mems(i)) = 1;
%end

% assigning cytoskeleton adjacencies
for inc = 0:Nmem-1 
    % membrane
    % basic membrane
    a1 = mod(inc,Nmem)+1;
    for sec = 1:2:round(Nmem/4)
        a2 = mod(inc+sec,Nmem)+1;
        A(a1,a2) = 1;
        A(a2,a1) = 1;
    end
end


K = k*ones(N,N);    % stiffness matrix

K(Nmem+1:Nact+Nmem,Nmem+1:Nact+Nmem)= kheavy*ones(size(Aact)); % Second block upper actin
K(Nmem+Nact+1:Nmem+2*Nact,Nmem+Nact+1:Nmem+2*Nact) = kheavy*ones(size(Aact));  % Second block upper actin
CA = [];


f = [12 480 5 840 28 120 7.5 5 7.6];            % forward rates
r = [12 1200 135 45.4 .0308 4 3.75 25 22.8];    % reverse rates
k = [27 10 15 5 16 15 10];                      % catalyst rates
x20 = [.15 .9285 9.6506 .0015 0 .3332 .2713 .013 2.8207 15.1793 23.9558 .0144 .0166 .0132 7.5];
x20 = x20';


plotsprings(x0,y0,A,0);
axis([-rmaj-30 rmaj+30 -rmin-25 rmin+25]);
pause(.1)
Fr = getframe;

% compute initial lengths of undeformed network
L0 = zeros(N,N);
for n = 1:N         % for each node
    for m = 1:N     % for each other node
        if(A(n,m)~=0)   % if connected
            L0(n,m) = sqrt((x0(n)-x0(m))^2+(y0(n)-y0(m))^2);    % find distance between nodes
        
        end
                
    end
end

% initial state vector
q0 = [x0; y0; u0; v0; x20];

% initial positions

normvel = 1; normvels = [];

tops = Nmem+Nact-(Nmy-3):Nmem+Nact;
centers = N-(Nmy-2):N-1;
bottoms = Nmem+2*Nact:-1:Nmem+2*Nact-(Nmy-3);

Fx0 = zeros(N,1);

Fcontract = 0;

% initial contractive forces from crossbridges
for i = 1:Nmy-2
    Fx0(tops(i)) = Fcontract;
    Fx0(bottoms(i)) = -Fcontract;
end

% run ode23 repeatedly until "steady state"
T(1) = 0;
Q = 1;

tinc = .25;

ii = Tmax/20;

for jj = 1:ii  %resets stepon and locations to which force is applied for each cycle 
    
tops = Nmem+Nact-(Nmy-3):Nmem+Nact;
centers = N-(Nmy-2):N-1;
bottoms = Nmem+2*Nact:-1:Nmem+2*Nact-(Nmy-3);

stepon = 1;
disp('Next cycle');

while(T(Q)<Tmax*(jj/ii))% && min(bottoms)>Nmem && min(tops)>17 && normvel>tol)

    [~,q] = ode23(@basicspringsrhs,[T(Q),T(Q)+tinc],q0);
    
 
    
    % extract out node locations and velocities
    x = q(end,0*N+1:1*N); y = q(end,1*N+1:2*N);
    u = q(end,2*N+1:3*N); v = q(end,3*N+1:4*N);
    chem = q(end,4*N+1:4*N+15);

 

    Fcontract = 700*(chem(13) + chem(14))/24;
    % new powerstroke every 3 timesteps
    if mod(T(Q),4*tinc) == 0 && stepon == 1
        % break old connections
        for i = 1:Nmy-2
            A(bottoms(i),centers(i)) = 0;
            A(centers(i),bottoms(i)) = 0;
            A(tops(i),centers(i)) = 0;
            A(centers(i),tops(i)) = 0;
            
            Fx0(bottoms(i)) = 0;
            Fx0(tops(i)) = 0;
        end
        bottoms = bottoms - 1;
        tops = tops - 1;
        % form new connections
        for i = 1:Nmy-2
            A(bottoms(i),centers(i)) = 1;
            A(centers(i),bottoms(i)) = 1;
            A(tops(i),centers(i)) = 1;
            A(centers(i),tops(i)) = 1;
            
            Fx0(bottoms(i)) = -Fcontract;
            Fx0(tops(i)) = Fcontract;
        end
    else
        Fx0(bottoms(i)) = -Fcontract;
        Fx0(tops(i)) = Fcontract;
    end
    
    if min(bottoms)==Nmem+Nact+1 || min(tops)==Nmem+1 
        Fx0 = zeros(N,1);   % turn off crossbridge forces
        stepon = 0;
    end
    
    
    
    % find maximum  magnitude of node velocity  
    normvel = max(sqrt(u.^2 + v.^2));
    normvels = [normvels;normvel];
    
    % update q0 (make vector from last row of q)
    q0 = transpose(q(end,:));
    
    % update time
    Q = Q + 1;
    T(Q) = T(Q-1)+tinc;
    ff = ff + 1;
    
    % plot nodes & springs
    plotsprings(x, y, A, T(Q));
    axis([-rmaj-30 rmaj+30 -rmin-25 rmin+25]);
    pause(.05)
    Fr(ff) = getframe;
    
end
end

figure(2)
plot(CA);
figure(3);
movie(Fr,1,12);
%movie2avi(Fr,'Smooth Muscle Cell Contraction 1 - Three Calcium Cycles');  
end

function plotsprings(x,y,A,t)




global Nmem

% plots the nodes
tt = num2str(t);
plot(x,y,'.k')
axis equal;
hold on
title(horzcat('time is ', tt, ' seconds'))
% plots the springs in between
[n1,n2] = find(A);
for n = 1:length(n1)
    id1 = n1(n);
    id2 = n2(n);
    x12 = [x(id1) x(id2)];
    y12 = [y(id1) y(id2)];
     if (id1<=Nmem && id2<=Nmem && abs(id1-id2)==1) || (id1==1 && id2==Nmem) || (id1==Nmem && id2==1)
        plot(x12,y12,'-r','linewidth',2)
    elseif id1<=Nmem && id2<=Nmem
        plot(x12,y12,'-c')
    else
        plot(x12,y12,'-b')
     end
end
hold off

end

function dq = basicspringsrhs(t,q)
% original lengths
global L0 A K M pins d p Nmem Fx0 Nact CA N f r k CaMax

% location of x variables, y variables, ..
xids = 0*N+1:1*N;
yids = 1*N+1:2*N;
uids = 2*N+1:3*N;
vids = 3*N+1:4*N;
chemids = 4*N+1:4*N+15;

% extract out node positions/velocities
x2 = q(xids);
y = q(yids);
u = q(uids);
v = q(vids);
x =  q(chemids);

% Equations
dx(1,1)  = dcalcium(t,CaMax);
dx(2,1)  = r(1)*x(4) + r(3)*x(6) + r(8)*x(9) - f(1)*x(1)*x(1)*x(2) - f(3)*x(2)*x(3) - f(8)*x(2)*x(10);
dx(3,1)  = r(3)*x(6) + r(4)*x(7) + r(5)*x(8) - f(3)*x(2)*x(3) - f(4)*x(3)*x(4) - f(5)*x(3)*x(5);
dx(4,1)  = f(1)*x(1)*x(1)*x(2) + r(2)*x(5) + r(4)*x(7) + f(9)*x(1)*x(1)*x(9) - r(1)*x(4) - f(2)*x(1)*x(1)*x(4) - f(4)*x(3)*x(4) - r(9)*x(4)*x(10);
dx(5,1)  = f(2)*x(1)*x(1)*x(4) + r(5)*x(8) - r(2)*x(5) - f(5)*x(3)*x(5);
dx(6,1)  = f(3)*x(2)*x(3) + r(6)*x(7) - r(3)*x(6) - f(6)*x(1)*x(1)*x(6);
dx(7,1)  = f(4)*x(3)*x(4) + f(6)*x(1)*x(1)*x(6) + r(7)*x(8) - r(4)*x(7) - r(6)*x(7) - f(7)*x(1)*x(1)*x(7);
dx(8,1)  = f(5)*x(3)*x(5) + f(7)*x(1)*x(1)*x(7) - r(5)*x(8) - r(7)*x(8);
dx(9,1)  = f(8)*x(2)*x(10) + r(9)*x(4)*x(10) - r(8)*x(9) - f(9)*x(1)*x(1)*x(9);
dx(10,1) = r(8)*x(9) + f(9)*x(1)*x(1)*x(9) - f(8)*x(2)*x(10) - r(9)*x(4)*x(10);
dx(11,1) = -k(1)*x(8)*x(11)/(k(2) + x(11)) + k(5)*x(15)*x(12)/(k(6) + x(12)) + k(7)*x(14);
dx(12,1) = k(1)*x(8)*x(11)/(k(2) + x(11)) - k(5)*x(15)*x(12)/(k(6) + x(12)) -k(3)*x(12) + k(4)*x(13);
dx(13,1) = k(3)*x(12) - k(4)*x(13) + k(1)*x(8)*x(14)/(k(2) + x(14)) - k(5)*x(15)*x(13)/(k(6) + x(13));
dx(14,1) = -k(1)*x(8)*x(14)/(k(2) + x(14)) + k(5)*x(15)*x(13)/(k(6) + x(13)) - k(7)*x(14);
dx(15,1) = 0;


% compute spring forces
Fx = zeros(N,1); Fy = zeros(N,1);

Px = zeros(N,1);
Py = zeros(N,1);

cellarea = polyarea(x2(1:Nmem),y(1:Nmem));
CA = [CA cellarea];

for n=1:N               % for each node
    for m=1:N           % for each other node
        if(A(n,m)~=0)   % if the nodes are attached
            % length of spring
            if (n<Nmem+2*Nact+1 && m<Nmem+2*Nact+1)
            L= sqrt( (x2(n)-x2(m))^2 + (y(n)-y(m))^2 );
            
            Fx(n) = Fx(n) + 10*K(m,n)*(L-L0(n,m))*(x2(m)-x2(n))/L;
            Fy(n) = Fy(n) + 10*K(m,n)*(L-L0(n,m))*(y(m)-y(n))/L;
            end
        end
        
        % pressure forces
        if n<=Nmem
            xdist = x2(n);
            ydist = y(n);
            dist = sqrt(xdist^2+ydist^2);
            
            % for pressure as function of area:         
            Pgen = p/cellarea;
            Pnode = Pgen/Nmem;
            Px(n) = Pnode*(xdist/dist);
            Py(n) = Pnode*(ydist/dist);

%             %  for constant pressure:
%             Px(n) = p*(xdist/dist);
%             Py(n) = p*(ydist/dist);     % Make sure one method of pressure is commented out
        end
    end
end

% compute right hand side vector --> produce vector r = q'
dq = zeros(4*N,1);
dq(xids) = u; % dx/dt = u
dq(yids) = v;  % dy/dt = v

% acceleration due to sum:
% of spring forces + pressure + applied forces + resistive damping
dq(uids) = (Fx+Fx0+Px)/M - d*u;     
dq(vids) = (Fy+Py)/M - d*v;

% fix the pinned nodes
dq(xids(pins)) = 0; dq(yids(pins)) = 0; dq(uids(pins)) = 0; dq(vids(pins)) = 0;

dq(chemids) = dx;

end

%
% function dcalcium
% 
% step wise sub-function which gives calcium levels
%
% Inputs
%   t       duration of simulation
%   max     maximum calcium levels
% 
% Outputs:
%   val     derivative of calcium concentration with respect to time at
%           instance of time input
%
% usage: type  dcalcium(t,cal)
% 
% example:  dcalcium(2,.6)
%

function val = dcalcium(t,max)

rest = .15;

t1 = 3;
t2 = 20;
t = mod(t,t2);
k = .5;

if t <= t1
    val = (max - rest)/(t1);
else
    val = -(max - rest)*k*exp(k*(-t+t1));
end

end