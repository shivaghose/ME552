% plot our data used to obtain friction data

VAL = xlsread('steadystateVelsV.xls');

DriverVs = VAL(:,2);
AngularVels = VAL(:,3);

% have numbers in terms of angular velocity, need motor torque constant to
% really use this data
% slope of Voltage / angular Velocity (radians / s)
Sforwards = 5.428888e-04;
Sreverse = 6.728054e-04;
% the y intercept, which represents the voltage needed to overcome friction
I0forwards = 0.463786;
I0reverse = -0.6071722;

VelsForwards = 0:1:700;
VForwards = ones(length(VelsForwards),1);

VelsReverse = -700:1:0;
VReverse = ones(length(VelsReverse),1);

for i = 1:length(VelsForwards)
    VForwards(i) = I0forwards + Sforwards * VelsForwards(i);
    VReverse(i) = I0reverse + Sreverse * VelsReverse(i);
end

figure
hold on;
plot(AngularVels,DriverVs,'ko');
plot(VelsForwards,VForwards,'b');
plot(VelsReverse,VReverse,'r');
title('Driver input Voltage vs Motor Angular Velocity');
xlabel('Angular Velocity (radians / s)')
ylabel('Driver Voltage (V)');
hold off

