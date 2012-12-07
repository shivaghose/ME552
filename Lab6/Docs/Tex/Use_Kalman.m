%This is a script to run various data filering algorithms. The script
%processes 1D signal data. Additionally we assume delta_t, the time b/w
%various measurements, to be constant. The folloowing filers are applied:
%   1)  Kalman Filter
%   2)  Low Pass Filter
%   3)  Moving Avg Filter
%Each filter's attributes can be set at each block. 


%Initialize:
str_data_fileName = 'Data/Sensitivity_Calibration.xls';
data = xlsread(str_data_fileName,'C:C');
x_init = data(1,1);
%data = xlsread('Z_Vbias_Data.xls','B:B');


%Start Kalman Filtering
p_init = 1;
steadyState_Data_Variance = 0.0019; % calculated from the static case
reduceLag_factor = 100;
r =  1/(steadyState_Data_Variance * reduceLag_factor);
q = 0.001;
Simple_Kalman_Init;

data_KF = zeros(size(data,1),1);
data_KF(1,1) = x_init;
for i = 2: size(data,1)
    [data_KF(i,1), p] = Simple_Kalman(data(i,1), data_KF(i-1,1), p, r, q);
end


%Low Pass Filter
tau = 0.05;
delta_t = 0.001;
freq = delta_t/tau;
data_LPF = filter(freq, [1 freq-1], data,x_init);


%Moving Average Filter
%Z_data_MAvgF = smooth(Z_data, 100); %Not causal
MAvgWindowSize = 100;
data_MAvgF = zeros(size(data,1),1);
data_MAvgF(1,1) = x_init;
for i = 1:size(data,1)
    if(i <= MAvgWindowSize)        
        data_MAvgF(i,1) = sum(data(i:-1:1,1))/i;
    else
        data_MAvgF(i,1) = sum(data(i:-1:i-MAvgWindowSize,1))/(MAvgWindowSize+1);
    end
end


%Plot results
t = 1 : size(data,1);
t = t';

figure
hold on
plot(t,data,'-r', t, data_LPF, '-g', t, data_MAvgF, '-b', t,data_KF,'-k')
str_KF_Spec  = ['Kalman Spec: ' 'r = ' num2str(r) ', q = ' num2str(q)...
    ', p_{init} = ' num2str(p_init)];
str_LPF_Spec = ['Low Pass Filter Spec: \tau = ' num2str(tau) ', \Delta T = ' num2str(delta_t)];
str_MAvgF_Spec = ['Moving Average Filter Spec: Window Size = ' num2str(MAvgWindowSize)];
figTitle = char(str_KF_Spec, str_LPF_Spec, str_MAvgF_Spec);
title(figTitle,'FontWeight','bold')
legend('Original Data','Low Pass Filter','Moving Avg Filter', 'Kalman Filter')
hold off




%Initial Variance - data(1,1) - data(1100,1)
init_index = 1100;
init_Var_data = var(data(1:init_index,1));
init_Var_KF = var(data_KF(1:init_index,1));
init_Var_LPF = var(data_LPF(1:init_index,1));
init_Var_MAvgF = var(data_MAvgF(1:init_index,1));

%Steady State Variance - Low Pass Filter should have caught up by now
ss_index = 32200;
lastIndex = size(data,1);
ss_Var_data = var(data(ss_index:lastIndex,1));
ss_Var_KF = var(data_KF(ss_index:lastIndex,1));
ss_Var_LPF = var(data_LPF(ss_index:lastIndex,1));
ss_Var_MAvgF = var(data_MAvgF(ss_index:lastIndex,1));


figure
hold on
subplot(4,1,1)
plot(t,data,'-r')
figTitle = 'Raw Data';
title(figTitle,'FontWeight','bold')
legend('Original Data')

subplot(4,1,2)
plot(t,data_KF,'-k')
figTitle = char(str_KF_Spec, ['Steady State Data Variance: ' num2str(ss_Var_KF)]);
title(figTitle,'FontWeight','bold')
legend('Kalman Filtered Data')

subplot(4,1,3)
plot(t,data_LPF,'-k')
figTitle = char(str_LPF_Spec, ['Steady State Data Variance: ' num2str(ss_Var_LPF)]);
title(figTitle,'FontWeight','bold')
legend('Low Pass Filtered Data')

subplot(4,1,4)
plot(t,data_MAvgF,'-k')
figTitle = char(str_MAvgF_Spec, ['Steady State Data Variance: ' num2str(ss_Var_MAvgF)]);
title(figTitle,'FontWeight','bold')
legend('Moving Avg Filtered Data')
hold off