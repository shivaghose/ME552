function [ filterData ] = lab6Filters( data, filterType, filterArgs )
%LAB6FILTERS Various filters for ME552 Lab 6 implemented as a function
%[ filterData ] = lab6Filters( data, filterType, filterArgs )
%   Inputs: 
%       data - Should be Nx1 Vector
%
%
%       filterType - can be 'KF', 'LPF' or 'MAvgF'
%
%
%       filterArgs:
%        Kalman Filter:  [p_init, steadyState_Data_Variance, reduceLag_factor, q]
%           Recommended Values:
%               p_init = 1;
%               steadyState_Data_Variance = From data
%               reduceLag_factor = 100
%               q = 0.001
%        
%        Low Pass Filter: [tau, delta_t]
%
%        High Pass Filter: [tau, delta_t]
%        
%        Moving Avg. Filter: [MAvgWindowSize]


filterTypeValid = strcmpi(filterType,'KF') + strcmpi(filterType,'LPF') +...
    strcmpi(filterType,'MAvgF') + strcmpi(filterType,'HPF') + ...
    strcmpi(filterType,'LPF_d') + strcmpi(filterType,'HPF_d');
assert(filterTypeValid>0, 'Error: filterType should be "KF", "LPF", "MAvgF" or "HPF"');

%Stad Initialization for all functions:
x_init = data(1,1);
%Kalman Filter
if strcmpi(filterType,'KF')
%     disp('Using Kalman Filter:');
    p_init = filterArgs(1);
    steadyState_Data_Variance = filterArgs(2);
    reduceLag_factor = 10;
    r =  1/(steadyState_Data_Variance * reduceLag_factor);
    q = 0.001;
    
%     disp(['Initializing r to ' num2str(r)]);
%     disp(['Initializing q to ' num2str(q)]);
%     disp(['Initializing p_init to ' num2str(p_init)]);
%     disp(['Initializing x_init to ' num2str(x_init)]);

    p = p_init + q;    
    
    filterData = zeros(size(data,1),1);
    filterData(1,1) = x_init;
    for i = 2: size(data,1)
        [filterData(i,1), p] = Simple_Kalman(data(i,1), filterData(i-1,1), p, r, q);
    end
    return;
end


%Low Pass Filter
if strcmpi(filterType,'LPF')
%     disp('Using Low Pass Filter');
    tau = filterArgs(1);
    delta_t = filterArgs(2);
    alpha = delta_t/tau;
    filterData = filter(alpha, [1 alpha-1], data,x_init);
    return;
end
if strcmpi(filterType,'LPF_d')
%     disp('Using Low Pass Filter');
    tau = filterArgs(1);    
    filterData = zeros(size(data,1),1);
    filterData(1,1) = x_init;
    for i = 2:size(data,1)
        delta_t = data(i,2) - data(i-1,2);
        alpha = delta_t/(delta_t + tau);
        filterData(i,1) = alpha*data(i,1) + (1 - alpha)*filterData(i-1,1);
    end
    return;
end
    

%High Pass Filter
if strcmpi(filterType,'HPF')
%     disp('Using Low Pass Filter');
    tau = filterArgs(1);
    delta_t = filterArgs(2);
    freq = delta_t/tau;
    filterData = filter([1-freq freq-1],[1 freq-1], data,x_init);
    return;
end
if strcmpi(filterType,'HPF_d')
    tau = filterArgs(1);        
    filterData = zeros(size(data,1),1);
    filterData(1,1) = x_init;
    for i = 2:size(data,1)
        delta_t = data(i,2) - data(i-1,2);
        alpha = tau / (tau + delta_t);
        filterData(i,1) = alpha * (filterData(i-1,1) + data(i,1) - data(i-1,1));
    end
end


if strcmpi(filterType,'MAvgF')
%     disp('Using Moving Average Filter');
    MAvgWindowSize = filterArgs(1);
    filterData = zeros(size(data,1),1);
    filterData(1,1) = x_init;
    for i = 1:size(data,1)
        if(i <= MAvgWindowSize)        
            filterData(i,1) = sum(data(i:-1:1,1))/i;
        else
            filterData(i,1) = sum(data(i:-1:i-MAvgWindowSize,1))/(MAvgWindowSize+1);
        end
    end
    return;
end

end