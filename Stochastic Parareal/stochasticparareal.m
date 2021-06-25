function [t,U,ERR,K,UG,UF] = stochasticparareal(f,tspan,u0,N,Ng,Nf,sample_rule,epsilon,m,sims)
%This function implements the stochastic parareal alogrithm to solve a
% system of first order ODEs in a time-parallel fashion. Sampling rules are
% based on multivariate normal or t-copula probability distributions.

%Inputs:
% f:           Function handle for function to be solved (i.e. f = @(t,u)([u(1);u(2)])
% tspan:       Time interval over which to integrate (i.e. [0,12])
% u0:          Initial conditions at tspan(1) (i.e. [0,1])
% N:           Number of temporal sub-intervals (i.e. N = 40)
% Ng:          Number of coarse time steps (i.e. Ndg = 40)
% Nf:          Number of fine times steps (i.e. Ndf = 4000)
% sample_rule: Sampling rule (i.e. 1,2 3 or 4)
% epsilon:     Error tolerance (i.e. 10^(-10))
% m:           Number of samples to take at each sub-interval (i.e. m = 10)
% sims:        Number of independent simulations of stochastic parareal (i.e. sims = 5)

%Outputs:
% t:           Vector of time sub-intervals (at which solutions located)
% U:           Cell array of ODE solutions for each sim
% ERR:         Cell array of errors between consecutive iterations for each sim
% K:           Vector of iterations taken until convergence for each sim
% UG:          Cell array of coarse solutions for each sim
% UF:          Cell array of fine solutions for each sim

%INITIALISATION
fprintf('Initialising... \n')

n = length(u0);                    %dimension of the ODE system
L = tspan(2) - tspan(1);           %length of interval
L_sub = L/N;                       %length of sub-interval
dT = L/Ng;                         %coarse time step
dt = L/Nf;                         %fine time step
t = (tspan(1):L_sub:tspan(2));     %sub-interval vector
t_shift = t(2:end);                %extra array for parfor loop

if mod(Ng,N)~=0 || mod(Nf,Ng)~=0  %number time steps must be multiples of each other
    error("Ndf must be a multiple of Ndg and Ndg must be a multiple of N")
end

% solution storage cell arrays for independent simulations
U = cell(sims,1);                  %parareal solutions
ERR = cell(sims,1);                %errors
K = zeros(sims,1);                 %iterations at which parareal converges
UG = cell(sims,1);                 %coarse solutions
UF = cell(sims,1);                 %fine solutions


% STOCHASTIC PARAREAL
fprintf('Step 1: Begin stochastic parareal... \n')

% loop if running more than one indpendent simulation of stochastic parareal.
for s = 1:sims
    
    % give an indication as to which simulation we're at for the console
    if s == 1
        fprintf('Step 2: Independent simulation number (out of %.0f): 1 ',sims)
    elseif s == sims
        fprintf('%.0f.',sims)
    else
        fprintf('%.0f ',s)
    end
    
    % solution matrices (structure: ith time step by kth iteration) (with
    %initial condition pre-defined across all k = 0)
    u = zeros(N+1,n*(N+1));                 u(1,:) = repmat(u0,1,N+1); %stores most accurate solutions (at sub-intervals at each k)
    uG = zeros(N+1,n*(N+1));                uG(1,:) = u(1,:);          %stores coarse solutions (at sub-intervals at each k)
    uF = zeros(N+1,n*(N+1));                uF(1,:) = u(1,:);          %stores fine solutions (at sub-intervals at each k)
    
    err = zeros(N+1,N);                                                %stores continuity errors (at sub-inervals at each k)
    I = 1;                                                             %index tracking which time intervals have converged
    k = 1;                                                             %iteration number
    
    dim_indices = (n*(k-1)+1:n*k);                                     %defines current indices
    dim_indices_next = ((n*k)+1:n*(k+1));                              %defines next indices
    
    %INITIAL COARSE SOLVE (k = 0)
    %Use G (coarse 'classic fourth-order') to find approximate initial conditions
    [~,temp] = RK(t(1):dT:t(end),u0,f,'classic fourth-order');        %solve the ODE
    uG(:,dim_indices) = temp(1:round(L_sub/dT):end,:);                %save solutions
    u(:,dim_indices) = uG(:,dim_indices);                             %save most up to date solution
    clear temp
    
    
    %INITIAL FINE SOLVE (k = 0)
    % use F with previously found initial conditions to solve in parallel
    parfor i = 1:N
        [~,temp] = RK((t(i):dt:t_shift(i)),uG(i,dim_indices),f,'classic fourth-order');    %solve ODE at fine resolution in each interval
        uF(i+1,dim_indices) = temp(end,:);                                                    %save the solution at sub-interval boundary
    end
    clear fine_trajecs temp
    
    %PREDICTOR-CORRECTOR STEP (for each sub-interval serially)
    for i = 1:N
        %First need to find uG for next iteration using coarse solver
        [~,temp] = RK((t(i):dT:t(i+1)),u(i,dim_indices_next),f,'classic fourth-order');
        uG(i+1,dim_indices_next) = temp(end,:);
        
        %Do the predictor-corrector step and save solution
        u(i+1,dim_indices_next) = uF(i+1,dim_indices) + uG(i+1,dim_indices_next) - uG(i+1,dim_indices);
    end
    clear temp
    
    %error catch (due to explicit solver - if used)
    if sum(isnan(uG),'all') ~= 0
        error('NaN values in uG - select higher number of coarse steps')
        %K(s,1) = NaN;
        %break
    end
    
    %%CONVERGENCE CHECKS
    %if an error is small enough up to a certain time interval, all solutions are saved and the
    %next k considers only unconverged chunks
    err(:,k) = vecnorm( u(:,dim_indices_next) - u(:,dim_indices), inf, 2);   % infinity norm at each time step
    
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
    
    
    %PARAREAL ITERATIONS
    for k = 2:N
        
        dim_indices = (n*(k-2)+1:n*(k-1));               %defines current indices
        dim_indices_next = ((n*(k-1))+1:n*k);            %defines next indices
        
        %CALCULATE CORRELATIONS BETWEEN DIMENSIONS (if n > 1)
        COR = repmat(diag(ones(1,n)),[1 1 N+1]);         %store correlation matrices for each time step in 3D array
        if n > 1
            if k > 2
                for i = I:N
                    temp = corr( fine_trajecs_end(:,:,i),'rows','complete' );
                    %ensure positive semi-definiteness of the matrix
                    temp = ((temp + temp')/2) - eye(n);
                    temp = temp - (temp > 0)*200*eps + (temp < 0)*200*eps;
                    COR(:,:,i) = temp + eye(n);
                end
            end
        end
        
        %INITIAL VALUE SAMPLING AND PROPAGATION
        %sample initial values and propagate in parallel.
        chunks = N - I;                                 %number of unconverged time chunks - 1
        sampled_initial_values = zeros((chunks*m)+1,n); %stores the m*chunks sampled initial values + the converged trajec in the Ith sub-interval
        propagated_F_trajecs = zeros((chunks*m)+1,n);        %stores m*chunks of propagated 'sampled initial values' + the converged trajec in the Ith sub-interval        
        parfor II = 1:(chunks*m)+1
            
            %in the first unconverged interval we have a correct initial
            %condition --> only need a single F run.
            if II == (chunks*m) + 1
                %integrate using F on appropriate interval and store in cell array
                [~,temp] = RK((t(I):dt:t_shift(I)),u(I,dim_indices_next),f,'classic fourth-order');
                sampled_initial_values(II,:) = u(I,dim_indices_next);
                propagated_F_trajecs(II,:) = temp(end,:);
            else
                %find time step over which to integrate over
                index_n = mod(II-1,chunks) + I + 1;
                
                %ensure one of the samples in each chunk is the predictor-corrector value we just found.
                if II <= chunks
                    %integrate using F on appropriate interval and store in cell array
                    [~,temp] = RK((t(index_n):dt:t_shift(index_n)),u(index_n,dim_indices_next),f,'classic fourth-order');
                    sampled_initial_values(II,:) = u(index_n,dim_indices_next);
                    propagated_F_trajecs(II,:) = temp(end,:);
                    
                    %otherwise we sample from the appropriate distribution
                else
                    if n == 1
                        if sample_rule == 1      %univariate normal
                            sample = normrnd(uF(index_n,dim_indices),abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices)));
                        elseif sample_rule == 2  %univariate normal
                            sample = normrnd(u(index_n,dim_indices_next),abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices)));
                        elseif sample_rule == 3  %univariate uniform
                            sample = random('uniform',uF(index_n,dim_indices) - sqrt(3)*abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices)),uF(index_n,dim_indices) + sqrt(3)*abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices)) );
                        elseif sample_rule == 4  %univariate uniform
                            sample = random('uniform',u(index_n,dim_indices_next) - sqrt(3)*abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices)),u(index_n,dim_indices_next) + sqrt(3)*abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices)) );
                        end
                    elseif n > 1
                        if sample_rule == 1      %mutivariate normal
                            sigma = abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices));
                            cov = repmat(sigma,n,1).*repmat(sigma',1,n).*COR(:,:,index_n);
                            sample = mvnrnd(uF(index_n,dim_indices),cov);
                        elseif sample_rule == 2  %multivariate normal
                            sigma = abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices));
                            cov = repmat(sigma,n,1).*repmat(sigma',1,n).*COR(:,:,index_n);
                            sample = mvnrnd(u(index_n,dim_indices_next),cov);
                        elseif sample_rule == 3  %multivariate t copula (with nu = 1)
                            sigma = abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices));
                            temp = copularnd('t',COR(:,:,index_n),1,1);
                            sample = 2*sqrt(3)*sigma.*temp + uF(index_n,dim_indices) - sqrt(3)*sigma;
                        elseif sample_rule == 4  %multivariate t copula (with nu = 1)
                            sigma = abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices));
                            temp = copularnd('t',COR(:,:,index_n),1,1);
                            sample = 2*sqrt(3)*sigma.*temp + u(index_n,dim_indices_next) - sqrt(3)*sigma;
                        end
                    end
                    
                    %integrate using F on appropriate interval and store in cell array
                    [~,temp] = RK((t(index_n):dt:t_shift(index_n)),sample,f,'classic fourth-order');
                    sampled_initial_values(II,:) = sample;
                    propagated_F_trajecs(II,:) = temp(end,:);
                end
            end
        end
        
        %LOCATING THE OPTIMAL INITAL VALUES
        time_index = mod((1:chunks*m)-1,chunks) + I + 1;      %locates correct interval indices from parfor loop
        opt_init_vals = NaN(N+1,n);                           %store optimal initial values
        fine_trajecs_end = zeros(m,n,N+1);                    %store the fine propagations of each sample (for correlation calculation in next k)
        for i = I:N 
            trajec_indices = find(time_index == i);           %find trajectory indices
            
            %in the I(th) interval, the trajectory is already optimal having come from a converged initial value
            if i == I
                uF(i+1,dim_indices_next) = propagated_F_trajecs(end,:);
                
                %for later unconverged intervals, select trajectory that is closest to the fine trajectory in the previous interval
            else
                %calculate errors in each dimension for each sampled trajec
                diffs =  abs( sampled_initial_values(trajec_indices,:) - uF(i,dim_indices_next) );
                fine_trajecs_end(:,:,i) = propagated_F_trajecs(trajec_indices,:); %store propagations for correlations later on
                
                %pick smallest 2-norm
                [~,min_index] = min(vecnorm(diffs,2,2));
                               
                %store solutions
                opt_init_vals(i,:) = sampled_initial_values(trajec_indices(min_index),:);
                uF(i+1,dim_indices_next) = propagated_F_trajecs(trajec_indices(min_index),:);
            end
        end
              
        %EXTRA G RUNS: need to run G from the new optimal initial values for
        %correction in the next step.
        for i = I+1:N
            [~,temp] = RK((t(i):dT:t_shift(i)),opt_init_vals(i,:),f,'classic fourth-order');
            uG(i+1,dim_indices_next) = temp(end,:);
        end
        clear temp
        
        dim_indices = (n*(k-1)+1:n*k);           %re-defines current indices
        dim_indices_next = ((n*k)+1:n*(k+1));    %re-defines next indices
        
        %PREDICTOR-CORRECTOR STEP (for each sub-interval serially)
        for i = I:N
            %First need to find uG for next iteration using coarse solver
            [~,temp] = RK((t(i):dT:t(i+1)),u(i,dim_indices_next),f,'classic fourth-order');
            uG(i+1,dim_indices_next) = temp(end,:);
            
            %Do the predictor-corrector step and save solution
            u(i+1,dim_indices_next) = uF(i+1,dim_indices) + uG(i+1,dim_indices_next) - uG(i+1,dim_indices);
        end
        
        %%CONVERGENCE CHECKS
        %if an error is small enough up to a certain time interval, all solutions are saved and the
        %next k considers only unconverged chunks
        err(:,k) = vecnorm( u(:,dim_indices_next) - u(:,dim_indices), inf, 2);   %error at each time step
        
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
        
        %break the stochastic parareal iteration if all steps converged
        if I == N + 1
            break
        end
    end
    
    %save solutions for this independent simulation and do next
    if s == 1
        U = u(:,(n+1):n*(k+1));
        ERR = err;
        K(s,1) = k;
        UG = uG(:,1:n*(k+1));
        UF = uF(:,1:n*(k+1));    
    else
        U{s,1} = u(:,(n+1):n*(k+1));
        ERR{s,1} = err;
        K(s,1) = k;
        UG{s,1} = uG(:,1:n*(k+1));
        UF{s,1} = uF(:,1:n*(k+1));
    end    
end

fprintf(' \n')
fprintf('Step 3: Stochastic Parareal complete. \n')

end
