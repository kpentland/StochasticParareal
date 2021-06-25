function [t,u,err,k] = parareal(f,tspan,u0,N,Ng,Nf,epsilon)
%This function will implement the parareal alogrithm for a system of first
% order ODEs. This version only updates time sub-intervals that have not
% converged.

%Inputs:
% f:           Function handle for function to be solved (i.e. f = @(t,u)([u(1);u(2)])
% tspan:       Time interval over which to integrate (i.e. [0,12])
% u0:          Initial conditions at tspan(1) (i.e. [0,1])
% N:           Number of 'proccesors' (temporal sub-intervals) (i.e. N = 40)
% Ng:          Number of coarse time steps (i.e. Ng = 40)
% Nf:          Number of fine times steps (i.e. Nf = 4000)
% epsilon:     Error tolerance (i.e. 10^(-10))

%Outputs:
% t:           Vector of time sub-intervals (at which solutions located)
% u:           Parareal solution to ODE system on the mesh given by 't'
% err:         Error at each time sub-interval at each k
% k:           Iterations taken until convergence

%INITIALISATION
fprintf('Initialising parareal... \n')

n = length(u0);                    %dimension of the ODE system
L = tspan(2) - tspan(1);           %length of interval
L_sub = L/N;                       %length of sub-interval
dT = L/Ng;                         %coarse time step
dt = L/Nf;                         %fine time step
t = (tspan(1):L_sub:tspan(2));     %time sub-interval vector  
t_shift = t(2:end);                %shifted time vector for parfor loop below

if mod(Ng,N)~=0 || mod(Nf,Ng)~=0  %number time steps must be multiples of each other
    fprintf("Nf must be a multiple of Ng and Ng must be a multiple of N - change time steps")
    return
end

% solution storage matrices (time step x system dimension*iteration)
u = zeros(N+1,n*(N+1));      %stores parareal solutions (at sub-interval time steps)
uG = zeros(N+1,n*(N+1));     %stores coarse solutions (at sub-interval time steps)
uF = zeros(N+1,n*(N+1));     %stores fine solutions (at sub-interval time steps)

err = zeros(N+1,N);          %error at each time slice
I = 1;                       %counts how many intervals have converged

% define the initial condition to be exact at start of each iteration
u(1,:) = repmat(u0,1,N+1);
uG(1,:) = u(1,:);
uF(1,:) = u(1,:);


%PARAREAL ALGORITHM

%Step 1 (k = 0): Use G (coarse solver) to find 'rough' initial conditions
fprintf('Stage 1: Inital coarse solve...')
[~,temp] = RK(t(1):dT:t(end),u0,f,'classic fourth-order');
uG(:,1:n) = temp(1:round(L_sub/dT):end,:); clear temp;
u(:,1:n) = uG(:,1:n);
fprintf('done. \n')


%Step 2 (k>0): Use F (fine solver) in parallel using current best initial
% conditions and update using the predictor-corrector formula.
for k = 1:N
    
    % give an indication as to which iteration we're at for the console
    if k == 1
        fprintf('Stage 2: Parareal iteration number (out of %.0f): 1 ',N)
    elseif k == N
        fprintf('%.0f.',N)
    else
        fprintf('%.0f ',k)
    end
    
    % use the previously found (or updated) initial conditions and
    % solve with the fine solver (in parallel)
    dim_indices = (n*(k-1)+1:n*k);        %current indices
    dim_indices_next = ((n*k)+1:n*(k+1)); %next indices
    parfor i = I:N
        [~,u_f] = RK((t(i):dt:t_shift(i)),u(i,dim_indices),f,'classic fourth-order');
        uF(i+1,dim_indices) = u_f(end,:);    %save the solution from final time step
    end
    
    
    %Predictor-corrector step (for each sub-interval serially)
    for i = I:N
        %First need to find uG for next iteration step using coarse solver
        [~,u_temp] = RK((t(i):dT:t(i+1)),u(i,dim_indices_next),f,'classic fourth-order');
        uG(i+1,dim_indices_next) = u_temp(end,:);
        
        %Do the predictor-corrector step and save final solution value
        u(i+1,dim_indices_next) = uF(i+1,dim_indices) + uG(i+1,dim_indices_next) - uG(i+1,dim_indices);
    end
    
    %error catch
    if sum(isnan(uG),'all') ~= 0
        error("NaN values in initial coarse solve - increase Ng!")
    end
    
    %%CONVERGENCE CHECKS
    %if an error is small enough up to a certain time interval, all solutions are saved and the
    %next k considers only unconverged chunks
    err(:,k) = vecnorm( u(:,dim_indices_next) - u(:,dim_indices), inf, 2)';   %error at each time step
    
    II = I;
    for p = II:N
        if p == II
            %we know I is now converged so we copy solutions to all future k
            u(p+1,(n*(k+1))+1:end) = repmat(u(p+1,dim_indices_next),1,N-k);
            uG(p+1,(n*(k+1))+1:end) = repmat(uG(p+1,dim_indices_next),1,N-k);
            uF(p+1,(n*k)+1:end) = repmat(uF(p+1,dim_indices),1,N-k+1);
            I = I + 1;
        elseif p > II
            if err(p,k) < epsilon
                %if further intervals beyond I are converged we can copy solutions to all future k
                u(p+1,(n*(k+1))+1:end) = repmat(u(p+1,dim_indices_next),1,N-k);
                uG(p+1,(n*(k+1))+1:end) = repmat(uG(p+1,dim_indices_next),1,N-k);
                uF(p+1,(n*k)+1:end) = repmat(uF(p+1,dim_indices),1,N-k+1);
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

%output the matrix containing the solutions after 1,2,3...k iterations
u = u(:,(n+1):n*(k+1));   

fprintf(' \n')
fprintf('Step 3: Parareal complete. \n')

end
