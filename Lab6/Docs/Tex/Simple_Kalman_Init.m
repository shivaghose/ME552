%   Initializes the values for the Simple_Kalman function.

if(~exist('r','var'))
    r = 1;     
end
if(~exist('q','var'))
    q = 1;    
end
if(~exist('p_init','var'))
    p_init = 1;    
end
if(~exist('x_init','var'))
    x_init = 1;    
end

disp(['Initializing r to ' num2str(r)]);
disp(['Initializing q to ' num2str(q)]);
disp(['Initializing p_init to ' num2str(p_init)]);
disp(['Initializing x_init to ' num2str(x_init)]);

p = p_init + q;
x = x_init;