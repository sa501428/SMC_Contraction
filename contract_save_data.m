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
% usage: type  contract_save_data(N,a,b,reps,xtra1,xtra2)  from matlab prompt
%
% example:  contract_save_data(100,300,50,2,1,1)
%

function contract_save_data(NN,a,b,reps,xtra1,xtra2)
global d K M A L0 N lx ly f_data t_data saved_data


N = NN;

% global parameters
% K is spring constant matrix
% M is node mass
% A is adjacency matrix
% L0 is rest lengths
% d is drag coefficient

k = 1;                  % default spring constant
maxk = round(N/10);     % spring constant as function of the # of nodes
K = k*ones(N);          % create spring constant matrix
M = 1;                  % mass of each node
A = zeros(N);           % preallocate adjacency matrix
L0 = zeros(N);          % preallocate matrix for unstretched lengths
d = .9;                 % drag coefficient
tinc = .25;              % time interval for integration
Tmax = 20;               % how long a contraction/relaxation takes
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
        len0 = sqrt((x0(a1)-x0(am))^2+(y0(a1)-y0(am))^2);
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
        -y0(list(j,2)))^2);
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


% run cycles for the various sets of force data
for j = 1:5
    disp('run number')
    j
    f_data = dlmread(strcat(num2str(j),'force_forceR2b.txt'));
    t_data = dlmread(strcat(num2str(j),'force_timeR2b.txt'));
    
    tids = find(t_data>10);
    f_data(tids) = 0;
    
    % reset data collector
    saved_data = [];

    
    % initial state vector
    q0 = [x0; y0; u0; v0];
    
    
    %L00 = L0;
    % the contraction/relaxation cycles
    T=0;
    
    % run ode23 until time is up
    while( T < Tmax)
        
        [~,q] = ode23(@basicspringsrhs,[T,T+tinc],q0);
        
        
        % else relaxation occurs
        % resting length is retored to longer length
        
        % extract out node locations
        x = q(end,0*N+1:1*N); y = q(end,1*N+1:2*N);
        %u = q(end,2*N+1:3*N); v = q(end,3*N+1:4*N);
        %chem = q(end,4*N+1:4*N+20);
        % update q0 (make vector from last row of q)
        q0 = transpose(q(end,:));
        
        
        
        % update time
        T = T+tinc;
        %miter = miter + 1;
        
        % plot nodes & springs
        
        % contraction is red colored
        figure(1)
        plotsprings(x, y, A2,'-r',T,j);
        
    end
    
    % save data to txt file
    dlmwrite(strcat('contract_data',num2str(j),'.txt'), saved_data)
    
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

function dq = basicspringsrhs(t,q)

global L0 A K N M d L saved_data

%location of x variables, y variables, ..
xids = 0*N+1:1*N;
yids = 1*N+1:2*N;
uids = 2*N+1:3*N;
vids = 3*N+1:4*N;

% extract out node positions/velocities
x2 = q(xids);
y = q(yids);
u = q(uids);
v = q(vids);

% get major/minor axes using least squares
AE = [x2 y].^2;
aaE= (AE'*AE)\AE'*ones(length(x2),1);
data_temp = sqrt(1./aaE);

% save data to matrix
saved_data = [saved_data; t data_temp']; 


%Fcontract = 50*(x(18) + x(19))/24;
list = [2 N/2+2 N/2 N 4 N/2+4 N/2-2 N-2 1 N/2+1];

% compute spring forces
Fx = zeros(N,1);
Fy = zeros(N,1);

for j = 1:length(list)
    l_n = list(j);
    my_l = sqrt(x2(l_n)^2 + y(l_n)^2);
    Fx(l_n) = Fx(l_n) - 40*force(t)*x2(l_n)/my_l;
    Fy(l_n) = Fy(l_n) - 40*force(t)*y(l_n)/my_l;
end


for n=1:N               % for each node
    for m=1:N           % for other nodes
        if(A(n,m)~=0)   % if the nodes are attached
            
            % length of spring
            L= sqrt( (x2(n)-x2(m))^2 + (y(n)-y(m))^2 );
            
            % spring forces
            Fx(n) = Fx(n) + K(m,n)*(L-L0(n,m))*(x2(m)-x2(n))/L;
            Fy(n) = Fy(n) + K(m,n)*(L-L0(n,m))*(y(m)-y(n))/L;
        end
        
        
    end
    
end


% compute right hand side vector --> produce vector r = q'
dq = zeros(4*N,1);
dq(xids) = u;                % dx/dt = u
dq(yids) = v;                % dy/dt = v

% acceleration due to sum of spring forces + resistive damping
dq(uids) = (Fx)/M - d*u;     % du/dt = a_x
dq(vids) = (Fy)/M - d*v;     % dv/dt = a_y

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
%   j       run number
%
% usage: type  plotsprings(x,y,A,str)
%
% example:  plotsprings(x, y, A,'-g',j);
%

function plotsprings(x,y,A,str,tt,j)

global lx ly

% spline interpolation
newx = [x x(1)];
newy = [y y(1)];
slen = zeros(size(newx));
for jj=2:length(newx)
  slen(jj) = sqrt((newx(jj)-newx(jj-1))^2 + ...
                  (newy(jj)-newy(jj-1))^2) + slen(jj-1);
end
ds = slen(end)/1000;
slenlong = 0:ds:slen(end);
xx = spline(slen,newx,slenlong);
yy = spline(slen,newy,slenlong);
plot(xx,yy,'k','LineWidth',2)
axis equal; axis([-lx lx -ly ly]);
hold on

plot([x(1) x(51)],[y(1) y(51)],'b','LineWidth',2)
plot([x(3) x(50)],[y(3) y(50)],'b','LineWidth',2)
plot([x(4) x(49)],[y(4) y(49)],'b','LineWidth',2)
plot([x(100) x(53)],[y(100) y(53)],'b','LineWidth',2)
plot([x(98) x(53)],[y(98) y(53)],'b','LineWidth',2)

% plots the springs in between
[n1,n2] = find(A);
for n = 1:length(n1)
    id1 = n1(n);
    id2 = n2(n);
    x12 = [x(id1) x(id2)];
    y12 = [y(id1) y(id2)];
    plot(x12,y12,str)
end
title(strcat('Time=',num2str(tt),' seconds, E_{40} case'),'fontsize',14)
xlabel('ISMC Length (\mu m)','fontsize',14)
ylabel('ISMC Height (\mu m)','fontsize',14)
set(gca,'XTickLabel',{'-150','-100','-50','0','50','100','150'})
set(gca,'YTickLabel',{'-10','-5','0','5','10'})
if (j==3)
pause
end
hold off

return


%
%
%

function f = force(t)
global f_data t_data

indx = find(t_data >= mod(t,20),1);
f = max(f_data(indx),0);

return