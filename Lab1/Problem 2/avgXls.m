function [avg] = avgXls(filename)
% averages out the columns of the given file

NUM = xlsread(filename);

meanV = mean(NUM);

avg = meanV(2);


end

