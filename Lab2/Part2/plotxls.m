function [] = plotxls(filename)
% plots 2nd column vs the first column of an excel file
VAL = xlsread(filename);

plot(VAL(:,1),VAL(:,2));
title(filename)
xlabel('X')
ylabel('Y')

end

