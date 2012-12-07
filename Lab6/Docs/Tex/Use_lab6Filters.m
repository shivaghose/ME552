%This is a script to run various data filering algorithms. The script
%processes 1D signal data. Additionally we assume delta_t, the time b/w
%various measurements, to be constant. The folloowing filers are applied:
%   1)  Kalman Filter
%   2)  Low Pass Filter
%   3)  Moving Avg Filter
%Each filter's attributes can be set at each block. 


%Initialize:
str_data_fileName = 'Sensitivity_Calibration.xls';
data = xlsread(str_data_fileName,'C:C');

%data = xlsread('Z_Vbias_Data.xls','B:B');


%Start Kalman Filtering
p_init = 1;
steadyState_Data_Variance = 0.0019; % calculated from the static case
reduceLag_factor = 10;
r =  1/(steadyState_Data_Variance * reduceLag_factor);
q = 0.001;
KF_args = [p_init, steadyState_Data_Variance, reduceLag_factor, q];
data_KF = lab6Filters(data, 'KF', KF_args);



%Low Pass Filter
tau = 0.05;
delta_t = 0.001;
LPF_args = [tau, delta_t];
data_LPF = lab6Filters(data, 'LPF', LPF_args);


%Moving Average Filter
%Z_data_MAvgF = smooth(Z_data, 100); %Not causal
MAvgWindowSize = 100;
MAvgF_args = [MAvgWindowSize];
data_MAvgF = lab6Filters(data, 'MAvgF', MAvgF_args);


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