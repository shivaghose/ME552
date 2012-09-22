[A, B] = meshgrid(10:100:5000, 1:0.2:5);

I = 0.1;
x = 0.002;

F = zeros(size(A));
for a = 1:size(A,1)
    for b = 1:size(A,2)
        F(a,b) = estimateForceSimple(A(a,b),B(a,b),I,x);
    end
end

figure
surf(A,B,F);
title('Model contour');
xlabel('A');
ylabel('B');
zlabel('F');
