%
% Run simulation to generate data for plotting ACH concentration as
% function of time and saving data as txt files. 
%

function ACH_save_data
global Xsim psf Cach Tach

% ACH simulation preallocation
Xsim = 100;                 % length of 2d plane for ach diffusion
Tsim = 80000;               % duration of simulation (20 sec)
maxsim = 5;
alpha = .1;
% gaussian blur for diffusion factor for 2d simulation
psf = fspecial('gaussian',3,.4);
    
% preallocate vectors hold ACH/time data
Tach = zeros(Tsim,1);
Cach = zeros(Tsim,maxsim);

for count = 1:Tsim
    Tach(count) = count*0.25*10^(-9);           % saving time at instance
end

for simnum = 1:maxsim
    Ysim = (simnum+1)*10;       % width of 2d plane for ach diffusion
    sim = zeros(Xsim, Ysim);    % creating simulation grid

    % ACH released at neuron end (1000 microMolar concentration)
    for xint = round(5*Xsim/11):round(6*Xsim/11)
            sim(xint,2) = 1000;    
    end
    %amount = 0;
    count = 8001;
    % running simulation
    while count < 80000
        if count > 8000
            sim = fixer(sim,Ysim,alpha);       % diffusion time step
            amount = (sim(50,Ysim-1) + sim(51,Ysim-1))/2; % concentration
            Cach(count,simnum) = amount; % saving concentration at instance
        end
        count = count + 1;
    end
end

% save data
dlmwrite('timeACH.txt', Tach)
dlmwrite('conACH.txt', Cach)

% plot data
plot(Tach,Cach(:,1),Tach,Cach(:,2),Tach,Cach(:,3),Tach,Cach(:,4),Tach,Cach(:,5))
legend('20 microns','30 microns','40 microns','50 microns','60 microns')
xlabel('Time (seconds)')
ylabel('Concentration (microMolar)')
title('Lower Concentrations at Further Distances')

return

function sim = fixer(sim0,Ysim,alpha)
% call global variables
global Xsim

% preliminary diffusion step without any assumed boundaries
%sim = imfilter(sim0,psf);

% extracting coefficients from psf
%c1 = psf(1,1);              % coeff for corners
%c2 = psf(1,2);              % coeff for adjacencies
%c3 = psf(2,2)+c2+2*c1;      % coeff for center

%% running correction for 2 boundaries along cell and neuron membranes 
sim = sim0;
for x = 2:Xsim-1
    for y = 2:Ysim-1
    % for bottom row (neuron end)
    sim(x,y) = sim0(x,y) + alpha*(-4*sim0(x,y)+sim0(x+1,y)+sim0(x-1,y)+sim0(x,y+1)+sim0(x,y-1));
    end
end

for x = round(4*Xsim/10):round(6*Xsim/10)
    % for bottom row (neuron end)
    sim(x,2) = sim0(x,2) + alpha*(-3*sim0(x,2)+sim0(x+1,2)+sim0(x-1,2)+sim0(x,3));
    % for top row (myocyte membrane)
    sim(x,Ysim-1) = sim0(x,Ysim-1) + alpha*(-3*sim0(x,Ysim-1)+sim0(x+1,Ysim-1)+sim0(x-1,Ysim-1)+sim0(x,Ysim-2));
end


return