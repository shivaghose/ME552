clc
clear all


str_data_fileName = 'Data/IMU_HorizontalFullSwing.xlsx';

%Get data sim time
time = xlsread(str_data_fileName,'A:A');
tStart = time(1);
tEnd = time(end) - tStart;
%delta_t = time(2) - tStart; %Sampling time
time = time - tStart;
delta_t_array = zeros(length(time),1);
for i = 2:length(time)
    delta_t_array(i,1) = time(i,1) - time(i-1,1);
end
delta_t_array(1,1) = delta_t_array(2,1);

delta_t = mean(delta_t_array);


%Import sensor data
V_Y_data = xlsread(str_data_fileName,'D:D');
V_Z_data = xlsread(str_data_fileName,'F:F');
V_Gyro_data = xlsread(str_data_fileName,'H:H');

%Import encoder data
encoder_theta = xlsread(str_data_fileName,'B:B');
encoder_theta = deg2rad(encoder_theta);
encoder_theta = mod2pi(encoder_theta);

[ZeroRate_Gyro_mean, ZeroRate_Gyro_var]  = ZeroRate(encoder_theta, V_Gyro_data);

t_interp = 0:0.0001:tEnd;
V_Z_data_interp = interpolate(time, V_Z_data, t_interp);
V_Y_data_interp  = interpolate(time, V_Y_data, t_interp);
V_Gyro_data_interp = interpolate(time, V_Gyro_data, t_interp);


%Create simulink I/P Data structs
V_Z = [t_interp' , V_Z_data_interp];
V_Y = [t_interp' , V_Y_data_interp];
V_Gyro = [t_interp' , V_Gyro_data_interp];

%Output is stored in theta_out
%We now need to compare theta_out with the encoder

%Set tau for the filters here:

tau = 1.1;

sim('SensorFusion3000');


%Static offset
[static_offset offset_index]= SensorFusionOffset(encoder_theta, theta_out);

theta_out_corrected = theta_out - static_offset;

figure
hold on
plot(time, encoder_theta, '--r', tout, theta_out_corrected, '-b', tout, theta_acc_out, '-.k', tout,theta_gyro_out, ':m')
title(['Angle Measurement Using Sensor Fusion with \tau = ' num2str(tau)],'FontWeight','bold')    
legend('Encoder','Sensor Fusion','Accelerometers Alone', 'Gyro Alone (Not Mod2Pi-ed)')
xlabel('Time (seconds)','FontSize',16)
ylabel('Arm Position (\theta) (rad)','FontSize',16)
hold off

% %Plot Error
% length_diff = length(encoder_theta)- length(theta_out_corrected);
% if(length_diff >= 0) 
%     error = encoder_theta(1:end - length_diff) - theta_out_corrected;
%     angles_abs = encoder_theta(1:end - length_diff);
% else
%     error = encoder_theta - theta_out_corrected(1:end - length_diff);
%     angles_abs = encoder_theta;
% end

% figure
% hold on
% plot(angles_abs,error, '-b')
% title('Error Magnitude v.s. Encoder Position','FontWeight','bold')    
% legend('Error Magnitude')
% xlabel('Absolute Angle)','FontSize',16)
% ylabel('Error Magnitude (\theta_{encoder} - \theta_{sensor fusion}) (rad)','FontSize',16)
% hold off

lpf_tf = tf([1],[tau 1]);
hpf_tf = tf([tau 0], [tau 1]);
figure
hold on
%bode(hpf_tf, lpf_tf);
bodeplot(hpf_tf,'r-',lpf_tf,'b--');
title(['High Pass Filter and Low Pass Filter response for with \tau = ' num2str(tau)],'FontWeight','bold')    
legend('High Pass Filter','Low Pass Filter')
hold off