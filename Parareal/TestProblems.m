%% Test problems
%Below are example test problems using the parareal algorithm.
%To run each test, simply run each section below independently.

%% Test Problem 1: One-dimensional nonlinear ODE

clear; close all; clc

%Inputs:
f = @(t,u)( sin(u)*cos(u) - 2*u + exp(-t/100)*sin(5*t) + log(1+t)*cos(t) );     %function handle for ODE
tspan = [0,100];                        %time interval
u0 = 1;                                 %intial conditions
N = 40;                                 %no. of time sub-intervals steps
Ng = 80;                                %no. of coarse steps (in each sub-interval)
Nf = Ng*100;                            %no. of fine steps (in each sub-interval)
epsilon = 10^(-10);                     %error tolerance 

%solve with parareal
[t,U,err,K] = parareal(f,tspan,u0,N,Ng,Nf,epsilon);

%solve the ODE in parallel with fine solver using initial conditions from parareal
dt = (tspan(2)-tspan(1))/Nf;    t_fine = (tspan(1):dt:tspan(2));

n = length(u0);                       %dimension of system
dim_indices = (n*(K-1)+1:n*K);        %final solution indices
U_fine = zeros(Nf+1,length(u0));
fine_trajecs = cell(N,1);
parfor i = 1:N
    [~,temp] = RK((t(i):dt:t(i+1)),U(i,dim_indices),f,'classic fourth-order');
    if i < N
    fine_trajecs{i,1} = temp(1:end-1,:);                 
    else
    fine_trajecs{i,1} = temp;                 
    end
end
U_fine(1:end,:) = vertcat(fine_trajecs{:,1});
    
    
%solve using the fine solver serially (for comparison)
[~,u_fine] = RK((tspan(1):dt:tspan(2)),u0,f,'classic fourth-order');

%plot both solutions together
figure(1)
hold on
plot(t_fine,u_fine,'k','LineWidth',2)
plot(t_fine(1:10:Nf+1),U_fine(1:10:Nf+1),'ob')
xlabel('$t$','Interpreter','latex'); ylabel('$u_1(t)$','Interpreter','latex');
grid on; box on;
xlim([0 18]); ylim([-2 2]);
xticks((0:3:18))
legend({'$\mathcal{F}$ solution','$\mathcal{P}$ solution'},'Interpreter','latex','location','northwest')


%% Test Problem 2: Two-dimensional Brusselator system

clear; close all; clc

%Inputs:
f = @(t,u)([1 + (u(1)^2)*u(2) - (3+1)*u(1); 3*u(1) - (u(1)^2)*u(2)]);     %function handle for ODE
tspan = [0,15.3];                      %time interval
u0 = [1,3.07];                         %intial conditions
N = 25;                                %no. of time sub-intervals steps
Ng = 25;                               %no. of coarse steps (in each sub-interval)
Nf = Ng*100;                           %no. of fine steps (in each sub-interval)
epsilon = 10^(-6);                     %error tolerance 

%solve with parareal
[t,U,err,K] = parareal(f,tspan,u0,N,Ng,Nf,epsilon);

%solve the ODE in parallel with fine solver using initial conditions from parareal
dt = (tspan(2)-tspan(1))/Nf;    t_fine = (tspan(1):dt:tspan(2));

n = length(u0);                       %dimension of system
dim_indices = (n*(K-1)+1:n*K);        %final solution indices
U_fine = zeros(Nf+1,length(u0));
fine_trajecs = cell(N,1);
parfor i = 1:N
    [~,temp] = RK((t(i):dt:t(i+1)),U(i,dim_indices),f,'classic fourth-order');
    if i < N
    fine_trajecs{i,1} = temp(1:end-1,:);                 
    else
    fine_trajecs{i,1} = temp;                 
    end
end
U_fine(1:end,:) = vertcat(fine_trajecs{:,1});

%solve using the fine solver serially (for comparison)
[~,u_fine] = RK((tspan(1):dt:tspan(2)),u0,f,'classic fourth-order');

%plot both solutions together in 2D phase space
figure(1)
hold on
plot(u_fine(:,1),u_fine(:,2),'k','LineWidth',2)  
plot(U_fine((1:50:Nf+1),end-1),U_fine((1:50:Nf+1),end),'ob')
xlabel('$u_1$', 'Interpreter','latex'); ylabel('$u_2$', 'Interpreter','latex');
axis([0 4 0 5])
xticks((0:0.5:4))
yticks((0:0.5:5))
legend({'$\mathcal{F}$ solution','$\mathcal{P}$ solution'},'Interpreter','latex','location','northeast')
box on; grid on
hold off


%% Test Problem 3: Three-dimensional Lorenz system

clear; close all; clc

%Inputs:
f = @(t,u)([10*(u(2) - u(1)); 28*u(1) - u(2) - u(1)*u(3); u(1)*u(2) - (8/3)*u(3)]);    %function handle for ODE
tspan = [0,18];                        %time interval
u0 = [-15,-15,20];                     %intial conditions
N = 50;                                %no. of time sub-intervals steps
Ng = 250;                              %no. of coarse steps (in each sub-interval)
Nf = Ng*75;                            %no. of fine steps (in each sub-interval)
epsilon = 10^(-8);                     %error tolerance 

%solve with parareal
[t,U,err,K] = parareal(f,tspan,u0,N,Ng,Nf,epsilon);

%solve the ODE in parallel with fine solver using initial conditions from parareal
dt = (tspan(2)-tspan(1))/Nf;    t_fine = (tspan(1):dt:tspan(2));

n = length(u0);                       %dimension of system
dim_indices = (n*(K-1)+1:n*K);        %final solution indices
U_fine = zeros(Nf+1,length(u0));
fine_trajecs = cell(N,1);
parfor i = 1:N
    [~,temp] = RK((t(i):dt:t(i+1)),U(i,dim_indices),f,'classic fourth-order');
    if i < N
    fine_trajecs{i,1} = temp(1:end-1,:);                 
    else
    fine_trajecs{i,1} = temp;                 
    end
end
U_fine(1:end,:) = vertcat(fine_trajecs{:,1});

%solve using the fine solver serially (for comparison)
[~,u_fine] = RK((tspan(1):dt:tspan(2)),u0,f,'classic fourth-order');


%plot the solution in 3D phase space
figure(1)
hold on
plot3(u_fine(:,1),u_fine(:,2),u_fine(:,3),'k','LineWidth',2)
plot3(U_fine((1:(Nf/Ng):Nf+1),end-2),U_fine((1:(Nf/Ng):Nf+1),end-1),U_fine((1:(Nf/Ng):Nf+1),end),'ob')
xlabel('$u_1$', 'Interpreter','latex'); 
ylabel('$u_2$', 'Interpreter','latex'); 
zlabel('$u_3$', 'Interpreter','latex');
view(45,10)
legend({'$\mathcal{F}$ solution','$\mathcal{P}$ solution'},'Interpreter','latex','location','northeast')
box on; grid on;
hold off

%% Test Problem 1: One-dimensional Bernoulli equation

clear; close all; clc

%Inputs:
f = @(t,u)(2*(u/(1+t)) - (t*u)^2);     %function handle for ODE
tspan = [0,10];                        %time interval
u0 = 2;                                %intial conditions
N = 20;                                %no. of time sub-intervals steps
Ng = 20;                               %no. of coarse steps (in each sub-interval)
Nf = Ng*100;                           %no. of fine steps (in each sub-interval)
epsilon = 10^(-10);                    %error tolerance 

%solve with parareal
[t,U,err,K] = parareal(f,tspan,u0,N,Ng,Nf,epsilon);

%solve the ODE in parallel with fine solver using initial conditions from parareal
dt = (tspan(2)-tspan(1))/Nf;    t_fine = (tspan(1):dt:tspan(2));

n = length(u0);                       %dimension of system
dim_indices = (n*(K-1)+1:n*K);        %final solution indices
U_fine = zeros(Nf+1,length(u0));
fine_trajecs = cell(N,1);
parfor i = 1:N
    [~,temp] = RK((t(i):dt:t(i+1)),U(i,dim_indices),f,'classic fourth-order');
    if i < N
    fine_trajecs{i,1} = temp(1:end-1,:);                 
    else
    fine_trajecs{i,1} = temp;                 
    end
end
U_fine(1:end,:) = vertcat(fine_trajecs{:,1});

%solve using the fine solver serially (for comparison)
[~,u_fine] = RK((tspan(1):dt:tspan(2)),u0,f,'classic fourth-order');

%plot both solutions together
figure(1)
hold on
plot(t_fine,u_fine,'k','LineWidth',2)              %fine solver solution
plot(t_fine(1:40:Nf+1),U_fine((1:40:Nf+1),end),'ob') %parareal solution
xlabel('$t$','Interpreter','latex'); ylabel('$u_1(t)$','Interpreter','latex');
grid on; box on;
legend('northeast',{'$\mathcal{F}$ solution','$\mathcal{P}$ solution'},'Interpreter','latex')
hold off


%% Test Problem 5: Two-dimensional 'square limit cycle' system 

clear; close all; clc

%Inputs:
f = @(t,u)([-sin(u(1))*(0.1*cos(u(1)) + cos(u(2))); -sin(u(2))*(0.1*cos(u(2)) - cos(u(1)))  ]);     %function handle for ODE
tspan = [0,60];                        %time interval
u0 = [1.5,1.5];                        %intial conditions
N = 30;                                %no. of time sub-intervals steps
Ng = 30;                               %no. of coarse steps (in each sub-interval)
Nf = Ng*100;                           %no. of fine steps (in each sub-interval)
epsilon = 10^(-8);                     %error tolerance 

%solve with parareal
[t,U,err,K] = parareal(f,tspan,u0,N,Ng,Nf,epsilon);

%solve the ODE in parallel with fine solver using initial conditions from parareal
dt = (tspan(2)-tspan(1))/Nf;    t_fine = (tspan(1):dt:tspan(2));

n = length(u0);                       %dimension of system
dim_indices = (n*(K-1)+1:n*K);        %final solution indices
U_fine = zeros(Nf+1,length(u0));
fine_trajecs = cell(N,1);
parfor i = 1:N
    [~,temp] = RK((t(i):dt:t(i+1)),U(i,dim_indices),f,'classic fourth-order');
    if i < N
    fine_trajecs{i,1} = temp(1:end-1,:);                 
    else
    fine_trajecs{i,1} = temp;                 
    end
end
U_fine(1:end,:) = vertcat(fine_trajecs{:,1});

%solve using the fine solver serially (for comparison)
[~,u_fine] = RK((tspan(1):dt:tspan(2)),u0,f,'classic fourth-order');

%plot both solutions together
figure(1)
hold on
plot(u_fine(:,1),u_fine(:,2),'k','LineWidth',2)
plot(U_fine((1:60:Nf+1),end-1),U_fine((1:60:Nf+1),end),'ob')
axis([-0.75 pi+0.75 -0.75 pi+0.75])
xlabel('$u_1$', 'Interpreter','latex'); ylabel('$u_2$', 'Interpreter','latex');
xticks((-0.5:0.5:4.5))
yticks((-0.5:0.5:4.5))
legend({'$\mathcal{F}$ solution','$\mathcal{P}$ solution'},'Interpreter','latex','location','northeast')
box on; grid on
hold off
