function [t_fine,u_fine,err,k] = parareal(f,n,tspan,u0,N,Ndg,Ndf,epsilon)
%This function will implement the parareal alogrithm for a system of first
% order ODEs. This version only updates time sub-intervals that have not
% converged.

%Inputs:
% f:           Function handle for function to be solved (i.e. f = @(t,u)([u(1);u(2)])
% n:           Dimensions of ODE system (i.e. n = 2)
% tspan:       Time interval over which to integrate (i.e. [0,12])
% u0:          Initial conditions at tspan(1) (i.e. [0,1])
% N:           Number of 'proccesors' (temporal sub-intervals) (i.e. N = 40)
% Ndg:         Number of coarse time steps (i.e. Ndg = 40)
% Ndf:         Number of fine times steps (i.e. Ndf = 4000)
% epsilon:     Error tolerance (i.e. 10^(-10))

%Outputs:
% t_fine:      Time interval at the fine resolution
% u_fine:      ODE solution at the fine resolution
% err:         Error at each time step at each k
% k:           Iterations taken until convergence

%INITIALISATION
fprintf('Initialising parareal... \n')

L = tspan(2) - tspan(1);           %length of interval
L_sub = L/N;                       %length of sub-interval
dt = L/Ndg;                        %coarse time step
t_sub = (tspan(1):L_sub:tspan(2)); %sub-interval vector
t_top = t_sub(2:end);              %extra array for parfor loop
dtt = L/Ndf;                       %fine time step
t_fine = (tspan(1):dtt:tspan(2));  %fine time interval

if mod(Ndg,N)~=0 || mod(Ndf,Ndg)~=0  %number time steps must be multiples of each other
    fprintf("Ndf must be a multiple of Ndg and Ndg must be a multiple of N")
    return
end

%The following solution matrices are structured: ith time step by kth iteration
u_fine = zeros(length(t_fine),n*(N+1)); %final solution matrix at fine resolution
u = zeros(N+1,n*(N+1));      %stores most accurate up to date solutions (at sub-interval time steps)
uG = zeros(N+1,n*(N+1));     %stores coarse solutions (at sub-interval time steps)
uF = zeros(N+1,n*(N+1));     %stores fine solutions (at sub-interval time steps)
err = zeros(N,N+1);          %error at each time slice
I = 1;                       %counts how many intervals have converged

%Define the initial condition to be exact at start of each iteration
u(1,:) = repmat(u0,1,N+1);
uG(1,:) = u(1,:);
uF(1,:) = u(1,:);
u_fine(1,:) = u(1,:);


%PARAREAL ALGORITHM

%Step 1 (k = 0): Use G (coarse 'classic fourth-order') to find 'rough' initial conditions
fprintf('Stage 1: Inital coarse solve...')
[~,temp] = RK(t_sub(1):dt:t_sub(end),u0,f,'classic fourth-order');
uG(:,1:n) = temp(1:round(L_sub/dt):end,:);
u(:,1:n) = uG(:,1:n);
clear temp
fprintf('done. \n')


%Step 2 (k>0): Use F (finer 'classic fourth-order') in parallel using these initial
%conditions and update using the predictor-corrector formula.
for k = 1:N
    
    % Give an indication as to where we are in the code for the console
    if k == 1
        fprintf('Stage 2: Parareal iteration number (out of %.0f): 1 ',N)
    elseif k == N
        fprintf('%.0f.',N)
    else
        fprintf('%.0f ',k)
    end
    
    %Now we use the previously found (or updated) initial conditions and
    %solve with the fine solver (in parallel)
    dim_indices = (n*(k-1)+1:n*k);        %current indices
    dim_indices_next = ((n*k)+1:n*(k+1)); %next indices
    fine_trajecs = cell(N,1);
    parfor i = I:N
        [~,u_f] = RK((t_sub(i):dtt:t_top(i)),u(i,dim_indices),f,'classic fourth-order');
        uF(i+1,dim_indices) = u_f(end,:);                 %save the solution from final time step in interval
        fine_trajecs{i,1} = u_f(2:end,:);                 %temporarily store solution in cell array
    end
    u_fine(((Ndf/N)*(I-1))+2:end,dim_indices) = vertcat(fine_trajecs{:,1});  %concatenate stored solutions to final output
    
    
    %Predictor-corrector step (for each sub-interval serially)
    for i = I:N
        %First need to find uG for next iteration step using coarse solver
        [~,u_temp] = RK((t_sub(i):dt:t_sub(i+1)),u(i,dim_indices_next),f,'classic fourth-order');
        uG(i+1,dim_indices_next) = u_temp(end,:);
        
        %Do the predic-correc step and save final solution value
        u(i+1,dim_indices_next) = uF(i+1,dim_indices) + uG(i+1,dim_indices_next) - uG(i+1,dim_indices);
    end
    
    %error catch
    if sum(isnan(uG),'all') ~= 0
        error("NaN values in initial coarse solve - pick more coarse time steps")
    end
    
    %Error checking
    err(k,:) = vecnorm( u(:,dim_indices_next) - u(:,dim_indices), inf, 2);   %error at each time step
    
    II = I;
    for p = II:N
        if p == II
            %we know I is now converged so we copy solutions to all future k
            u(p+1,(n*(k+1))+1:end) = repmat(u(p+1,dim_indices_next),1,N-k);
            uG(p+1,(n*(k+1))+1:end) = repmat(uG(p+1,dim_indices_next),1,N-k);
            uF(p+1,(n*k)+1:end) = repmat(uF(p+1,dim_indices),1,N-k+1);
            u_fine(((Ndf/N)*(p-1))+2:((Ndf/N)*p)+1,(n*k)+1:end) = repmat(u_fine(((Ndf/N)*(p-1))+2:((Ndf/N)*p)+1,dim_indices),1,N-k+1);
            I = I + 1;
        elseif p > II
            if err(k,p) < epsilon
                %if further intervals beyond I are converged we can copy solutions to all future k
                u(p+1,(n*(k+1))+1:end) = repmat(u(p+1,dim_indices_next),1,N-k);
                uG(p+1,(n*(k+1))+1:end) = repmat(uG(p+1,dim_indices_next),1,N-k);
                uF(p+1,(n*k)+1:end) = repmat(uF(p+1,dim_indices),1,N-k+1);
                u_fine(((Ndf/N)*(p-1))+2:((Ndf/N)*p)+1,(n*k)+1:end) = repmat(u_fine(((Ndf/N)*(p-1))+2:((Ndf/N)*p)+1,dim_indices),1,N-k+1);
                I = I + 1;
            else
                break
            end
        end
    end
    
    %break the parareal iteration if all steps converged
    if I == N + 1
        break
    end
end

%output the fine solution up to the kth cycle.
u_fine = u_fine(:,1:n*k);

fprintf(' \n')
fprintf('Step 3: Parareal complete. \n')

end
