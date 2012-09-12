% script for reading in excel files and generating the graphs for it

filename = 'Set1_g.xls';

Num = xlsread(filename);

cntr = Num(:,1);
samplingPeriod = 0.001;

Vout = Num(:,2);

% now compute our inputs to graph along side
A = [3; 6; 0; 0];
f = [4; 4; 0; 0];
offset = [0; 0; 0; 0];
phase = [0; pi/2; 0; 0];

T = cntr * samplingPeriod;

V = zeros(length(Vout),length(A));

for i = 1 : size(V,2)
    for t = 1 : size(V,1)
        V(t,i) = A(i) * sin(2*pi*f(i) * T(t) + phase(i)) + offset(i);
    end
end

a = [1; 2; -3; -4];
Sum = sinSuperPosition(a, A, f, phase, offset ,T);

figure;
hold on;
plot(T,V(:,1),'b--');
plot(T,V(:,2),'r--');
%plot(T,V(:,3),'b--');
%plot(T,V(:,4),'r--');
plot(T,Sum,'k--');
plot(T,Vout,'k');
legend('V1','V2','Vout','Expected Vout'); % ,'V3','V4'
axis([0 T(end) (min([min(V) min(Vout) min(Sum)]) - 0.1) (max([max(V) max(Vout) max(Sum)])+ 0.1)]);
xlabel('Time (seconds)');
ylabel('Voltage (V)');

title('Problem 2g Set1')
hold off