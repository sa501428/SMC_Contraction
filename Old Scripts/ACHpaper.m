% function to do 2D diffusion with reaction

clear all
close all

n = 100;
m = 20;
width = .02:.01:.06;

% space discretization
dx = .1/n;
dy = .1/n;
dt = .1;
ttime = 80;
Diffusion = 4.0e-4;
cfl = Diffusion*dt/dx^2;
TAch = zeros(800,5);
concen = 90;
k = 4.0e-3;


for lala=1:5
    x = 0:dx:width(lala);
    y = 0:dy:.1;

    % create mesh of grid points
    [X,Y] = meshgrid(x,y);
    % the acetylcholine level
    nint = size(X,1)-2;
    mint = size(X,2)-2;
    A2 = zeros(nint*mint,mint*nint);
    Xint = X(2:end-1,2:end-1);
    Yint = Y(2:end-1,2:end-1);
    
    
% set up coefficient matrix
for i=2:nint-1
    for j=2:mint-1
        sp = (i-1)*mint+j;
        sp2 = i*mint+j;
        sp0 = (i-2)*mint+j;
        A2(sp,sp) = 1+4*cfl;
        A2(sp,sp+1) = -cfl;
        A2(sp,sp-1) = -cfl;
        A2(sp,sp2) = -cfl;
        A2(sp,sp0) = -cfl;
    end
end

% boundary conditions 
% bottom row
for j=2:mint-1
    sp = j;
    sp2 = mint+j;
    A2(sp,sp) = 1+3*cfl;
    A2(sp,sp+1) = -cfl;
    A2(sp,sp-1) = -cfl;
    A2(sp,sp2) = -cfl;
end

% top row
for j=2:mint-1
    sp = (nint-1)*mint+j;
    sp0 = (nint-2)*mint+j;
    A2(sp,sp) = 1+3*cfl;
    A2(sp,sp+1) = -cfl;
    A2(sp,sp-1) = -cfl;
    A2(sp,sp0) = -cfl;
end

% left side
for i=2:nint-1
    sp = (i-1)*mint+1;
    sp0 = (i-2)*mint+1;
    sp2 = i*mint+1;
    A2(sp,sp) = 1+3*cfl;
    A2(sp,sp+1) = -cfl;
    A2(sp,sp0) = -cfl;
    A2(sp,sp2) = -cfl;
end


% right side
for i=2:nint-1
    sp = (i-1)*mint+mint;
    sp0 = (i-2)*mint+mint;
    sp2 = i*mint+mint;
    A2(sp,sp) = 1+3*cfl+k*dt/dx;
    A2(sp,sp-1) = -cfl;
    A2(sp,sp0) = -cfl;
    A2(sp,sp2) = -cfl;
end

% corners
A2(1,1) = 1+2*cfl;
A2(1,2) = -cfl;
A2(1,mint+1) = -cfl;
A2(mint,mint) = 1+2*cfl+k*dt/dx;
A2(mint,mint-1) = -cfl;
A2(mint,2*mint) = -cfl;
A2((nint-1)*mint+1,(nint-1)*mint+1)=1+2*cfl;
A2((nint-1)*mint+1,(nint-1)*mint+2)=-cfl;
A2((nint-1)*mint+1,(nint-2)*mint+1)=-cfl;
A2(nint*mint,nint*mint)=1+2*cfl+k*dt/dx;
A2(nint*mint,nint*mint-1)=-cfl;
A2(nint*mint,(nint-1)*mint)=-cfl;

% initial conditions for Ach
b2 = zeros(nint*mint,1);


Ach2 = reshape(b2,mint,nint)';
sA2 = sparse(A2);
tc = 1;

for t=0:dt:ttime/10
    figure(1)
    clf
    newb2 = sA2\b2;
    Ach2 = reshape(newb2,mint,nint)';   
    TAch(tc,lala) = mean(Ach2(1:end,end));
    contourf(Xint,Yint,Ach2)
    caxis([0 concen])
    colorbar
    axis equal
    axis off
    b2=newb2;
    t
    tc = tc+1;
end

w = 4;
h = 20;
for i=ceil(nint/2)-h/2:ceil(nint/2)+h/2
    for j=1:w
        sp = (i-1)*mint+j;
        b2(sp) = concen;
    end
end

for t=ttime/10+dt:dt:ttime
    newb2 = sA2\b2;
    Ach2 = reshape(newb2,mint,nint)';   
    TAch(tc,lala) = mean(Ach2(1:end,end));
    contourf(Xint,Yint,Ach2)
    caxis([0 concen])
    colorbar
    axis equal
    axis off
    b2=newb2;
    t
    tc = tc+1;
end
end

figure(2)
plot(0:dt:ttime,TAch(:,1),'k','LineWidth',2)
hold on
plot(0:dt*20:ttime,TAch(1:20:end,2),'ro')
plot(0:dt:ttime,TAch(:,3),'g--')
plot(0:dt*20:ttime,TAch(1:20:end,4),'b*')
plot(0:dt*20:ttime,TAch(1:20:end,5),'mx')
plot(0:dt:ttime,TAch(:,2),'r')
plot(0:dt:ttime,TAch(:,4),'b')
plot(0:dt:ttime,TAch(:,5),'m')
xlabel('Time (\mu s)','fontsize',14)
ylabel('ACH (\mu M)','fontsize',14)
legend('NE','E_{30}','E_{40}','E_{50}','E_{60}','Location','northeast')

maxA = [max(TAch(:,1)),max(TAch(:,2)),max(TAch(:,3)),max(TAch(:,4)),max(TAch(:,5))];
synap = [20,30,40,50,60];

figure(3)
plot(synap,maxA,'ko','MarkerSize',5)
hold on
plot(synap,maxA,'k','LineWidth',2)
xlabel('Synapse Width (nm)','fontsize',14)
ylabel('Peak, Avg. ACh concentration (\mu M)','fontsize',14)

figure(4)
bar(maxA)
set(gca,'XTickLabel',{'20','30','40','50','60'},'fontsize',14);
xlabel('Synapse Width (nm)','fontsize',14)
ylabel('Peak, Avg. ACh Concentration (\mu M)','fontsize',14)

newtime = 0:0.001:80;
oldtime = 0:0.1:80;
y1 = interp1(oldtime,TAch(:,1),newtime);
y2 = interp1(oldtime,TAch(:,2),newtime);
y3 = interp1(oldtime,TAch(:,3),newtime);
y4 = interp1(oldtime,TAch(:,4),newtime);
y5 = interp1(oldtime,TAch(:,5),newtime);

YAch = [y1',y2',y3',y4',y5'];
size(YAch)

dlmwrite('conACHR2b.txt', YAch);

