%
% Generates graphs of major and minor axes vs time using data from
% contract_save_data.m
%

function contractGraph

% import each of the data files manually
j = 1;
data1 = dlmread(strcat('contract_data',num2str(j),'.txt'));

j = 2;
data2 = dlmread(strcat('contract_data',num2str(j),'.txt'));

j = 3;
data3 = dlmread(strcat('contract_data',num2str(j),'.txt'));

j = 4;
data4 = dlmread(strcat('contract_data',num2str(j),'.txt'));

j = 5;
data5 = dlmread(strcat('contract_data',num2str(j),'.txt'));
q = 50;
% plot major axis vs time
figure(1)
clf
plot(data1(:,1),data1(:,2),'k','LineWidth',2)
hold on
plot(data2(1:q:end,1),data2(1:q:end,2),'ro')
plot(data3(:,1),data3(:,2),'g--','LineWidth',2)
plot(data4(1:q:end,1),data4(1:q:end,2),'b*')
plot(data5(1:q:end,1),data5(1:q:end,2),'mx')
plot(data2(:,1),data2(:,2),'r')
plot(data4(:,1),data4(:,2),'b')
plot(data5(:,1),data5(:,2),'m')
xlabel('Time (seconds)','fontsize',14)
ylabel('Major Axis Length (\mu m)','fontsize',14)
legend('NE','E_{30}','E_{40}','E_{50}','E_{60}','Location','southeast')

mj(1) = min(data1(:,2))
mj(2) = min(data2(:,2))
mj(3) = min(data3(:,2))
mj(4) = min(data4(:,2))
mj(5) = min(data5(:,2))

figure(2)
clf
bar(mj/3)
set(gca,'XTickLabel',{'20','30','40','50','60'},'fontsize',14);
xlabel('Synapse Width (nm)','fontsize',14)
ylabel('Maximum Percent (%) Contraction','fontsize',14) 

mj/3


% plot minor axis vs time
% figure
% plot(data1(:,1),data1(:,3),data2(:,1),data2(:,3),data3(:,1),data3(:,3),data4(:,1),data4(:,3),data5(:,1),data5(:,3))
% legend('20 microns','30 microns','40 microns','50 microns','60 microns')
% xlabel('Time (seconds)')
% ylabel('Minor Axis')
% title('Minor Axis Length at Further Distances')

return