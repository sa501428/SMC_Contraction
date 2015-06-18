%
% mechan.m
%
% This code models the contraction of a myocyte give certain conditions and
% dimensions 
% 
% Inputs:
%   N       number of node
%   a       length of major axis
%   b       length of minor axis
%   reps    number of contraction relaxation cycles
%   xtra1   additional adjacent cytoskeletal suport (1 if used, else 0)
%   xtra2   additional cross cytoskelatal support (1 if used, else 0)
%
% usage: type  mechan(N,a,b,reps,xtra1,xtra2)  from matlab prompt
% 
% example:  mechan(100,300,50,2,1,1)
%

function mechan2(N,a,b,reps,xtra1,xtra2)

% global parameters
% K is spring constant matrix
% M is node mass
% A is adjacency matrix
% L0 is rest lengths
% d is drag coefficient
global K M A L0 d lx ly

k = 1;                  % default spring constant
maxk = round(N/10);     % spring constant as function of the # of nodes
K = k*ones(N);          % create spring constant matrix
M = 1;                  % mass of each node
A = zeros(N);           % preallocate adjacency matrix
L0 = zeros(N);          % preallocate matrix for unstretched lengths
d = .9;                 % drag coefficient
tinc = .5;              % time interval for integration
increment = 0;          % counts how many cycles have occured
Tmax = 10;               % how long a contraction/relaxation takes
lx = round(1.1*(a+b)); % for calculating plot dimensions
ly = round((.25*a+b)); % for calculating plot dimensions

x0 = zeros(N,1);        % x-coordinates of N nodes
y0 = zeros(N,1);        % y-coordinates of N nodes
u0 = zeros(N,1);        % x-component of velocity of N nodes
v0 = zeros(N,1);        % y-component of velocity of N nodes

% setting the coordinates
for inc = 0:N-1
    
    % polar equation of ellipse
    theta = 2*pi*inc/N;
    r = (a*b)/sqrt((b*cos(theta))^2+(a*sin(theta))^2);
    
    % x and y coordinates from polar data
    x0(inc+1,1) = round(r*cos(theta));
    y0(inc+1,1) = round(r*sin(theta));
end
A2 = A;
x0(N/2,1)
y0(N/2,1)
pause
% assigning membrane adjacencies
for inc = 0:N-1
    
    % basic membrane
    a1 = mod(inc,N)+1;
    %a2 = mod(inc+1,N)+1;
    a4 = mod(inc+7,N)+1;
    %A(a1,a2) = 1;
    %A(a2,a1) = 1;
    A(a1,a4) = 1;
    A(a4,a1) = 1;
    A2(a1,a4) = 1;
    A2(a4,a1) = 1;
    
    % optional additional support
    if xtra1 == 1
        more = round(N/10);
        am = mod(inc+more,N)+1;
        A(a1,am) = 1;
        A(am,a1) = 1;
        A2(a1,am) = 1;
        A2(am,a1) = 1;
        
        am = mod(inc+2*more,N)+1;
        A(a1,am) = 1;
        A(am,a1) = 1;
        A2(a1,am) = 1;
        A2(am,a1) = 1;
    end
    
    % optional additional support
    if xtra2 == 1
        more = round(N/2);
        am = mod(inc+more,N)+1;
        A(a1,am) = 1;
        A(am,a1) = 1;
        len0 = sqrt((x0(a1)-x0(am))^2+(y0(a1)-y0(am))^2)
        K(a1,am) = .3*100/len0;
        K(am,a1) = .3*100/len0;
    end
end

% This assigns AM fibers
list = [2 N/2+2; N/2 N; 4 N/2+4; N/2-2 N-2; 1 N/2+1];

% cross bridges adjacencies and make bridges stronger
for j = 1:length(list)
    A(list(j,1),list(j,2))= 1;
    A(list(j,2),list(j,1))= 1;
    len0 = sqrt((x0(list(j,1))-x0(list(j,2)))^2+(y0(list(j,1))...
              -y0(list(j,2)))^2)
    K(list(j,1),list(j,2))= 100*maxk/len0;
    K(list(j,2),list(j,1))= 100*maxk/len0;
end

% compute initial lengths of undeformed network
for n = 1:N
    for m = 1:N
        if(A(n,m)~=0)
            L0(n,m) = sqrt((x0(n)-x0(m))^2+(y0(n)-y0(m))^2);
        end
    end
end

L0save = L0;

% initial state vector
q0 = [x0; y0; u0; v0];

% initial plot
plotsprings(x0,y0,A2,'-b');
pause(1)
%L00 = L0;
% the contraction/relaxation cycles
while increment < reps
    % whole numbers cause contraction
    % resting length is shortened for contraction
    
    
    % time starts anew
    T = 0;
    miter = 0;
    miter2 = 0;

    % run ode23 until time is up
    while( T < Tmax)

        [~,q] = ode23(@basicspringsrhs,[T,T+tinc],q0);

        if mod(increment,1) == 0
            multiplier = .9;
            miter = miter + 1;
            Tmax = 10;
        else
            multiplier = 1/.95;
            miter2 = miter2 + 1;
            Tmax = 20;
        end
%         [miter, miter2]
%         increment
%         pause
        
        if miter < 10 && miter2 < 10
            if (mod(increment,1)==0)
               for j = 1:length(list)
                L0(list(j,1),list(j,2))= L0(list(j,1),list(j,2))*multiplier;
                L0(list(j,2),list(j,1))= L0(list(j,2),list(j,1))*multiplier;
               end
            else
               L0 = L0save;
            end
        end
        % else relaxation occurs
        % resting length is retored to longer length
        
        % extract out node locations
        x = q(end,1:N);
        y = q(end,N+1:2*N);
        % update q0 (make vector from last row of q)
        q0 = transpose(q(end,:));

        % update time
        T = T+tinc;
        %miter = miter + 1;

        % plot nodes & springs

        % contraction is red colored
        if mod(increment,1) == 0
            plotsprings(x, y, A2,'-r');
        % relaxation is green
        else
            plotsprings(x, y, A2,'-g');
        end

        pause(.005)
    end
    
    % half cycle completed
    increment = increment + .5;
    pause(.5)
    
end

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

function r = basicspringsrhs(~,q)

global K M A L0 d

N = length(q)/4;

% location of x variables, y variables, ..
xids = 0*N+1:1*N;
yids = 1*N+1:2*N;
uids = 2*N+1:3*N;
vids = 3*N+1:4*N;

% extract out node positions/velocities
x = q(xids);
y = q(yids);
u = q(uids);
v = q(vids);

% compute spring forces
Fx = zeros(N,1);
Fy = zeros(N,1);

for n=1:N               % for each node
    for m=1:N           % for other nodes
        if(A(n,m)~=0)   % if the nodes are attached
            
            % length of spring
            L= sqrt( (x(n)-x(m))^2 + (y(n)-y(m))^2 );

            % spring forces
            Fx(n) = Fx(n) + K(m,n)*(L-L0(n,m))*(x(m)-x(n))/L;
            Fy(n) = Fy(n) + K(m,n)*(L-L0(n,m))*(y(m)-y(n))/L;
        end
         

    end
        
end

% compute right hand side vector --> produce vector r = q'
r = zeros(4*N,1);
r(xids) = u;                % dx/dt = u
r(yids) = v;                % dy/dt = v

% acceleration due to sum of spring forces + resistive damping
r(uids) = (Fx)/M - d*u;     % du/dt = a_x
r(vids) = (Fy)/M - d*v;     % dv/dt = a_y

return

%
% function plotsprings
% 
% sub-function which plots the nodes and springs for the cell
%
% Inputs:
%   x       x coordinates
%   y       y coordinates
%   A       adjacency matrix
%   str     string controlling line colors   
%
% usage: type  plotsprings(x,y,A,str)
% 
% example:  plotsprings(x, y, A,'-g');
%

function plotsprings(x,y,A,str)

global lx ly

% plots the nodes
plot(x,y,'.k')
axis equal; axis([-lx lx -ly ly]);
hold on

% plots the springs in between
[n1,n2] = find(A);
for n = 1:length(n1)
    id1 = n1(n);
    id2 = n2(n);
    x12 = [x(id1) x(id2)];
    y12 = [y(id1) y(id2)];
    plot(x12,y12,str)
end

hold off

return