function [t_fine,U,ERROR,K,UG,UF] = stochasticparareal(f,n,tspan,u0,N,Ndg,Ndf,sample_rule,epsilon,m,sims)
%This function implements the stochastic parareal alogrithm to solve a
% system of first order ODEs in a time-parallel fashion. Sampling rules are
% based on multivariate normal or t-copula probability distributions.

%Inputs:
% f:           Function handle for function to be solved (i.e. f = @(t,u)([u(1);u(2)])
% n:           Dimension of ODE system (i.e. n = 2)
% tspan:       Time interval over which to integrate (i.e. [0,12])
% u0:          Initial conditions at tspan(1) (i.e. [0,1])
% N:           Number of temporal sub-intervals (i.e. N = 40)
% Ndg:         Number of coarse time steps (i.e. Ndg = 40)
% Ndf:         Number of fine times steps (i.e. Ndf = 4000)
% sample_rule: Sampling rule (i.e. 1,2 3 or 4)
% epsilon:     Error tolerance (i.e. 10^(-10))
% m:           Number of samples to take at each sub-interval (i.e. m = 10)
% sims:        Number of independent simulations of stochastic parareal (i.e. sims = 5)

%Outputs:
% t_fine: Time interval at the fine resolution
% U:      Cell array of ODE solutions at the fine resolution for each sim
% ERROR:  Cell array of errors between consecutive iterations for each sim
% K:      Vector of iterations taken until convergence for each sim
% UG:     Cell array of coarse solutions for each sim
% UF:     Cell array of fine solution (at sub-interval boundaries) for each sim


%INITIALISATION
fprintf('Initialising... \n')

if mod(Ndg,N)~=0 || mod(Ndf,Ndg)~=0  %number time steps must be multiples of each other
    error("Ndf must be a multiple of Ndg and Ndg must be a multiple of N")
end

L = tspan(2) - tspan(1);           %length of interval
L_sub = L/N;                       %length of sub-interval
dt = L/Ndg;                        %coarse time step
dtt = L/Ndf;                       %fine time step

t_sub = (tspan(1):L_sub:tspan(2)); %sub-interval vector
t_top = t_sub(2:end);              %extra array for parfor loop
t_fine = (tspan(1):dtt:tspan(2));  %fine time interval

U = cell(sims,1);                  %store solutions from independent simulations
ERROR = cell(sims,1);              %store errors
K = zeros(sims,1);                 %store iterations at which parareal converges
UG = cell(sims,1);                 %store coarse solutions from independent simulations
UF = cell(sims,1);                 %store fine solutions from independent simulations


% STOCHASTIC PARAREAL
fprintf('Step 1: Begin stochastic parareal... \n')

%Only used if running more than one indpendent simulation of stochastic parareal.
for s = 1:sims
    
    % Give an indication as to where we are in the code for the console
    if s == 1
        fprintf('Step 2: Independent simulation number (out of %.0f): 1 ',sims)
    elseif s == sims
        fprintf('%.0f.',sims)
    else
        fprintf('%.0f ',s)
    end
    
    %solution matrices (structure: ith time step by kth iteration) (with
    %initial condition pre-defined across all k = 0)
    u = zeros(N+1,n*(N+1));                 u(1,:) = repmat(u0,1,N+1); %stores most accurate solutions (at sub-intervals at each k)
    uG = zeros(N+1,n*(N+1));                uG(1,:) = u(1,:);          %stores coarse solutions (at sub-intervals at each k)
    uF = zeros(N+1,n*(N+1));                uF(1,:) = u(1,:);          %stores fine solutions (at sub-intervals at each k)
    u_fine = zeros(length(t_fine),n*(N+1)); u_fine(1,:) = u(1,:);      %final solution matrix (finest resolution at each k)
    err = zeros(N+1,N);                                                %stores continuity errors (at sub-inervals at each k)
    I = 1;                                                             %index tracking which time intervals have converged
    k = 1;                                                             %iteration number
    
    dim_indices = (n*(k-1)+1:n*k);                                     %defines current indices
    dim_indices_next = ((n*k)+1:n*(k+1));                              %defines next indices
    
    %INITIAL COARSE SOLVE (k = 0)
    %Use G (coarse 'classic fourth-order') to find approximate initial conditions
    [~,ug] = RK(t_sub(1):dt:t_sub(end),u0,f,'classic fourth-order');   %solve the ODE
    uG(:,dim_indices) = ug(1:round(L_sub/dt):end,:);                   %save solutions
    u(:,dim_indices) = uG(:,dim_indices);                              %save most up to date solution
    clear ug
    
    
    %INITIAL FINE SOLVE (k = 0)
    %Use F with previously found initial conditions to solve in parallel
    fine_trajecs = zeros(Ndf/N,n,N);
    parfor i = 1:N
        [~,uf] = RK((t_sub(i):dtt:t_top(i)),uG(i,dim_indices),f,'classic fourth-order');    %solve ODE at fine resolution in each interval
        uF(i+1,dim_indices) = uf(end,:);                                                    %save the solution at sub-interval boundary
        fine_trajecs(:,:,i) = uf(2:end,:);                                                  %store solution at fine resolution in 3D array
    end
    u_fine(2:end,dim_indices) = reshape(permute(fine_trajecs, [2 1 3]), n, [])';            %concatenate stored trajectories to save
    clear fine_trajecs
    
    %PREDICTOR-CORRECTOR STEP (for each sub-interval serially)
    for i = 1:N
        %First need to find uG for next iteration using coarse solver
        [~,ug] = RK((t_sub(i):dt:t_sub(i+1)),u(i,dim_indices_next),f,'classic fourth-order');
        uG(i+1,dim_indices_next) = ug(end,:);
        
        %Do the predictor-corrector step and save solution
        u(i+1,dim_indices_next) = uF(i+1,dim_indices) + uG(i+1,dim_indices_next) - uG(i+1,dim_indices);
    end
    clear ug
    
    %error catch (due to explicit solver)
    if sum(isnan(uG),'all') ~= 0
        error('NaN values in uG - select higher number of coarse steps')
        %K(s,1) = NaN;
        %break
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
            u_fine(((Ndf/N)*(p-1))+2:((Ndf/N)*p)+1,(n*k)+1:end) = repmat(u_fine(((Ndf/N)*(p-1))+2:((Ndf/N)*p)+1,dim_indices),1,N-k+1);
            I = I + 1;
        elseif p > II
            if err(p,k) < epsilon
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
                    temp = temp - (temp > 0)*100*eps + (temp < 0)*100*eps;
                    COR(:,:,i) = temp + eye(n);
                end
            end
        end
        
        %INITIAL VALUE SAMPLING AND PROPAGATION
        %sample initial values and propagate in parallel.
        chunks = N - I;                            %number of unconverged time chunks - 1
        sampledF_trajecs = cell((chunks*m)+1,1);   %stores m*chunks trajectories of propagated 'sampled initial values' + the converged trajec in the Ith sub-interval        
        parfor II = 1:(chunks*m)+1
            
            %in the first unconverged interval we have a correct initial
            %condition --> only need a single F run.
            if II == (chunks*m) + 1
                %integrate using F on appropriate interval and store in cell array
                [~,uf] = RK((t_sub(I):dtt:t_top(I)),u(I,dim_indices_next),f,'classic fourth-order');
                sampledF_trajecs{II,:} = uf;
            else
                %find time step over which to integrate over
                index_n = mod(II-1,chunks) + I + 1;
                
                %ensure one of the samples in each chunk is the predictor-corrector value we just found.
                if II <= chunks
                    %integrate using F on appropriate interval and store in cell array
                    [~,uf] = RK((t_sub(index_n):dtt:t_top(index_n)),u(index_n,dim_indices_next),f,'classic fourth-order');
                    sampledF_trajecs{II,:} = uf;
                    
                    %otherwise we sample from the appropriate distribution
                else
                    if n == 1
                        if sample_rule == 1      %univariate normal
                            sample_init_vals = mvnrnd(uF(index_n,dim_indices),corr2cov(abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices)),COR(:,:,index_n)));
                        elseif sample_rule == 2  %univariate normal
                            sample_init_vals = mvnrnd(u(index_n,dim_indices_next),corr2cov(abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices)),COR(:,:,index_n)));
                        elseif sample_rule == 3  %univariate uniform
                            sample_init_vals = random('uniform',uF(index_n,dim_indices) - sqrt(3)*abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices)),uF(index_n,dim_indices) + sqrt(3)*abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices)) );
                        elseif sample_rule == 4  %univariate uniform
                            sample_init_vals = random('uniform',u(index_n,dim_indices_next) - sqrt(3)*abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices)),u(index_n,dim_indices_next) + sqrt(3)*abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices)) );
                        end
                    elseif n > 1
                        if sample_rule == 1      %mutivariate normal
                            sigma = abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices));
                            cov = repmat(sigma,n,1).*repmat(sigma',1,n).*COR(:,:,index_n);
                            sample_init_vals = mvnrnd(uF(index_n,dim_indices),cov);
                        elseif sample_rule == 2  %multivariate normal
                            sigma = abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices));
                            cov = repmat(sigma,n,1).*repmat(sigma',1,n).*COR(:,:,index_n);
                            sample_init_vals = mvnrnd(u(index_n,dim_indices_next),cov);
                        elseif sample_rule == 3  %multivariate t copula (with nu = 1)
                            sigma = abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices));
                            temp = copularnd('t',COR(:,:,index_n),1,1);
                            sample_init_vals = 2*sqrt(3)*sigma.*temp + uF(index_n,dim_indices) - sqrt(3)*sigma;
                        elseif sample_rule == 4  %multivariate t copula (with nu = 1)
                            sigma = abs(uG(index_n,dim_indices_next) - uG(index_n,dim_indices));
                            temp = copularnd('t',COR(:,:,index_n),1,1);
                            sample_init_vals = 2*sqrt(3)*sigma.*temp + u(index_n,dim_indices_next) - sqrt(3)*sigma;
                        end
                    end
                    
                    %integrate using F on appropriate interval and store in cell array
                    [~,uf] = RK((t_sub(index_n):dtt:t_top(index_n)),sample_init_vals,f,'classic fourth-order');
                    sampledF_trajecs{II,:} = uf;
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
                uF(i+1,dim_indices_next) = sampledF_trajecs{end,1}(end,:);
                u_fine((Ndf/N)*(i-1)+1:((Ndf/N)*i)+1,dim_indices_next) = sampledF_trajecs{end,1};
                
                %for later unconverged intervals, select trajectory that is closest to the fine trajectory in the previous interval
            else
                %calculate errors in each dimension for each sampled trajec
                diffs = zeros(m,n);
                for j = 1:m
                    diffs(j,:) = abs( sampledF_trajecs{trajec_indices(j),1}(1,:) - uF(i,dim_indices_next) );
                    fine_trajecs_end(j,:,i) = sampledF_trajecs{trajec_indices(j),1}(end,:);
                end
                
                %pick smallest 2-norm
                [~,min_index] = min(vecnorm(diffs,2,2));
                
                %select the optimal trajectory across dimensions
                optimal_trajecF = sampledF_trajecs{trajec_indices(min_index),1};
                
                %store solutions
                opt_init_vals(i,:) = optimal_trajecF(1,:);
                uF(i+1,dim_indices_next) = optimal_trajecF(end,:);
                u_fine((Ndf/N)*(i-1)+2:(((Ndf/N)*i)+1),dim_indices_next) = optimal_trajecF(2:end,:);
            end
        end
              
        %EXTRA G RUNS: need to run G from the new optimal initial values for
        %correction in the next step.
        for i = I+1:N
            [~,ug] = RK((t_sub(i):dt:t_top(i)),opt_init_vals(i,:),f,'classic fourth-order');
            uG(i+1,dim_indices_next) = ug(end,:);
        end
        
        dim_indices = (n*(k-1)+1:n*k);           %re-defines current indices
        dim_indices_next = ((n*k)+1:n*(k+1));    %re-defines next indices
        
        %PREDICTOR-CORRECTOR STEP (for each sub-interval serially)
        for i = I:N
            %First need to find uG for next iteration using coarse solver
            [~,ug] = RK((t_sub(i):dt:t_sub(i+1)),u(i,dim_indices_next),f,'classic fourth-order');
            uG(i+1,dim_indices_next) = ug(end,:);
            
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
                u_fine(((Ndf/N)*(p-1))+2:((Ndf/N)*p)+1,(n*k)+1:end) = repmat(u_fine(((Ndf/N)*(p-1))+2:((Ndf/N)*p)+1,dim_indices),1,N-k+1);
                I = I + 1;
            elseif p > II
                if err(p,k) < epsilon
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
        
        %break the stochastic parareal iteration if all steps converged
        if I == N + 1
            break
        end
    end
    
    %save solutions for this independent simulation and do next (iff s > 1)
    U{s,1} = u_fine(:,1:dim_indices_next(end));
    ERROR{s,1} = err;
    K(s,1) = k;
    UG{s,1} = uG(:,1:dim_indices_next(end));
    UF{s,1} = uF(:,1:dim_indices_next(end));    
end

fprintf(' \n')
fprintf('Step 3: Stochastic Parareal complete. \n')

end
