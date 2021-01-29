function [t,u] = RK(t,u0,f,method)
%This function is an explicit Runge-Kutta ODE solver that solves systems 
% of ODES using either the 'classic fourth-order' or 'explicit midpoint'
% method depending on user preference.

%Inputs:
% t     : Time interval (i.e. t = (0:0.1:10))
% u0    : Initial conditions (i.e u0 = [1,3.07])
% f     : Function handle for ODEs to be solved (i.e. f = @(t,u)([1 +
% (u(1)^2)*u(2) - (3+1)*u(1); 3*u(1) - (u(1)^2)*u(2)]) for the Brusselator
% system).
% method: Runge-Kutta method to be used (i.e. 'classic fourth-order')

%Outputs:
% t    : Time interval (same as input)
% u    : Solution in matrix form (time steps by dimension of problem)

%Select the method to be used.
switch method
    case 'classic fourth-order'
        a = [0,0,0,0;0.5,0,0,0;0,0.5,0,0;0,0,1,0];
        b = [1/6,1/3,1/3,1/6];
        c = [0,0.5,0.5,1];
    case 'explicit midpoint'
        a = [0,0;0.5,0];
        b = [0,1];
        c = [0,0.5];
    otherwise
        fprintf('Error: Please define the Runge Kutta method. \n')
end

%SOLVING THE EQUATIONS.

%Pre-define the solution matrix and specify initial conditions
u = zeros(length(u0),length(t));
u(:,1) = u0;

%Iterate over each time step explicitly
for n = 1:length(t)-1
    
    %This function carries out the iterative step of a general form
    % of the Runge-Kutta method with inputs: (time step, initial time,
    % intial condition, function, coefficient matrices).
    
    h = t(n+1)-t(n);            %length of time step
    dim = length(u0);           %dimension of the ODE problem
    S = length(b);              %order of the RK method (2nd/4th)
    k = zeros(dim,S);           %matrix for other k values
    k(:,1) = h*f(t(n),u(:,n));  %definition of k1
    
    %calculate the coefficients k
    for i = 2:S
        temp = zeros(dim,1);
        for j = 1:i-1
            temp = temp + a(i,j)*k(:,j);
        end
        k(:,i) = h*f(t(n) + c(i)*h,u(:,n) + temp);
    end
    
    %calculate the final solution  
    u(:,n+1) = u(:,n) + sum(b.*k,2);
    
end
u = u';  %transpose solution matrix
end



