% 
clear all
close all

% Tach = dlmread('timeACH.txt');
% Tach = Tach * 10^6;
% Cach1 = dlmread('conACH.txt');
% Cach2 = dlmread('conACHMuh.txt');
% Cach3 = dlmread('conACHR1.txt');
% Cach4 = dlmread('conACHR2.txt');
% Cach5 = dlmread('conACHR3.txt');
% Cach6 = dlmread('conACHR4.txt');
% Cach7 = dlmread('conACHR5.txt');

q = 2000;
linecol = ['k' 'ro' 'g--' 'b*' 'mx'];

figure(1)
clf
f_data1 = dlmread('1force_forceR2b.txt');
t_data1 = dlmread('1force_timeR2b.txt');
f_data2 = dlmread('2force_forceR2b.txt');
t_data2 = dlmread('2force_timeR2b.txt');
f_data3 = dlmread('3force_forceR2b.txt');
t_data3 = dlmread('3force_timeR2b.txt');
f_data4 = dlmread('4force_forceR2b.txt');
t_data4 = dlmread('4force_timeR2b.txt');
f_data5 = dlmread('5force_forceR2b.txt');
t_data5 = dlmread('5force_timeR2b.txt');
ids1 = find(t_data1>9.9,1)
ids2 = find(t_data2>9.9,1)
ids3 = find(t_data3>9.9,1)
ids4 = find(t_data4>9.9,1)
ids5 = find(t_data5>9.9,1)

fs1 = f_data1(ids1);
fe = 0;
ni1 = ids1/q;
m1 = -fs1/10;
dt = 10.0/ni1;

plot(t_data1(1:ids1),f_data1(1:ids1),'k','LineWidth',2)
pause
hold on
plot(t_data2(1:q:ids2),f_data2(1:q:ids2),'ro')
plot(t_data3(1:ids3),f_data3(1:ids3),'g--','LineWidth',2)
plot(t_data4(1:q:ids4),f_data4(1:q:ids4),'b*')
plot(t_data5(1:q:ids5),f_data5(1:q:ids5),'mx')
plot(t_data2(1:ids2),f_data2(1:ids2),'r')
plot(t_data4(1:ids4),f_data4(1:ids4),'b')
plot(t_data5(1:ids5),f_data5(1:ids5),'m')
fs1 = f_data1(ids1);
fe = 0;
ni1 = ids1/q;
m1 = -fs1/10;
dt = 10.0/ni1;
nt = 10.0/dt;
for k=1:nt
    t1(k) = 10 + dt*(k-1);
    f1(k) = m1*(t1(k)-10)+fs1;
end
plot(t1,f1,'k','LineWidth',2)
pause

fs2 = f_data2(ids2);
ni2 = ids2/q;
m2 = -fs2/10;
dt = 10.0/ni2;
nt = 10.0/dt;
for k=1:nt
    t2(k) = 10 + dt*(k-1);
    f2(k) = m2*(t2(k)-10)+fs2;
end
plot(t2,f2,'ro')

fs3 = f_data3(ids3);
ni3 = ids3/q;
m3 = -fs3/10;
dt = 10.0/ni3;
nt = 10.0/dt;
for k=1:nt
    t3(k) = 10 + dt*(k-1);
    f3(k) = m3*(t3(k)-10)+fs3;
end
plot(t3,f3,'g--','LineWidth',2)
fs4 = f_data4(ids4);
ni4 = ids4/q;
m4 = -fs4/10;
dt = 10.0/ni4;
nt = 10.0/dt;
for k=1:nt
    t4(k) = 10 + dt*(k-1);
    f4(k) = m4*(t4(k)-10)+fs4;
end
plot(t4,f4,'b*')
fs5 = f_data5(ids5);
ni5 = ids5/q;
m5 = -fs5/10;
dt = 10.0/ni5;
nt = 10.0/dt;
for k=1:nt
    t5(k) = 10 + dt*(k-1);
    f5(k) = m5*(t5(k)-10)+fs5;
end
plot(t5,f5,'mx')

xlabel('Time (s)','fontsize',14)
ylabel('Force (milliNewtons)','fontsize',14)
legend('NE','E_{30}','E_{40}','E_{50}','E_{60}','Location','northeast')

q = 4000;
clear t1 t2 t3 t4 t5
clear f1 f2 f3 f4 f5
figure(2)
clf
f_data1 = dlmread('1force_CaR2b.txt');
t_data1 = dlmread('1force_timeR2b.txt');
f_data2 = dlmread('2force_CaR2b.txt');
t_data2 = dlmread('2force_timeR2b.txt');
f_data3 = dlmread('3force_CaR2b.txt');
t_data3 = dlmread('3force_timeR2b.txt');
f_data4 = dlmread('4force_CaR2b.txt');
t_data4 = dlmread('4force_timeR2b.txt');
f_data5 = dlmread('5force_CaR2b.txt');
t_data5 = dlmread('5force_timeR2b.txt');
ids1 = find(t_data1>9.9,1)
ids2 = find(t_data2>9.9,1)
ids3 = find(t_data3>9.9,1)
ids4 = find(t_data4>9.9,1)
ids5 = find(t_data5>9.9,1)

fs1 = f_data1(ids1);
fe = 0;
ni1 = ids1/1000;
m1 = -fs1/10;
dt = 10.0/ni1;

plot(t_data1,f_data1,'k','LineWidth',2)
hold on
plot(t_data2(1:q:end),f_data2(1:q:end),'ro')
plot(t_data3,f_data3,'g--','LineWidth',2)
plot(t_data4(1:q:end),f_data4(1:q:end),'b*')
plot(t_data5(1:q:end),f_data5(1:q:end),'mx')
plot(t_data2,f_data2,'r')
plot(t_data4,f_data4,'b')
plot(t_data5,f_data5,'m')

fs1 = f_data1(ids1);
fe = 0;
ni1 = ids1/q;
m1 = (0.2-fs1)/10;
dt = 10.0/ni1;
nt = 10.0/dt;
for k=1:nt
    t1(k) = 10 + dt*(k-1);
    f1(k) = m1*(t1(k)-10)+fs1;
end
plot(t1,f1,'k','LineWidth',2)

fs2 = f_data2(ids2);
ni2 = ids2/q;
m2 = (0.2-fs2)/10;
dt = 10.0/ni2;
nt = 10.0/dt;
for k=1:nt
    t2(k) = 10 + dt*(k-1);
    f2(k) = m2*(t2(k)-10)+fs2;
end
plot(t2,f2,'ro')

fs3 = f_data3(ids3);
ni3 = ids3/q;
m3 = (0.2-fs3)/10;
dt = 10.0/ni3;
nt = 10.0/dt;
for k=1:nt
    t3(k) = 10 + dt*(k-1);
    f3(k) = m3*(t3(k)-10)+fs3;
end
plot(t3,f3,'g--','LineWidth',2)
fs4 = f_data4(ids4);
ni4 = ids4/q;
m4 = (0.2-fs4)/10;
dt = 10.0/ni4;
nt = 10.0/dt;
for k=1:nt
    t4(k) = 10 + dt*(k-1);
    f4(k) = m4*(t4(k)-10)+fs4;
end
plot(t4,f4,'b*')
fs5 = f_data5(ids5);
ni5 = ids5/q;
m5 = (0.2-fs5)/10;
dt = 10.0/ni5;
nt = 10.0/dt;
for k=1:nt
    t5(k) = 10 + dt*(k-1);
    f5(k) = m5*(t5(k)-10)+fs5;
end
plot(t5,f5,'mx')

xlabel('Time (s)','fontsize',14)
ylabel('Calcium (\mu M)','fontsize',14)
legend('NE','E_{30}','E_{40}','E_{50}','E_{60}','Location','northeast')

