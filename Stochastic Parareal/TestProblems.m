%% Test problems
%Below are example test problems using the stochastic parareal algorithm.
%To run each test, simply run each section below independently.

%% Test Problem 1: Bernoulli equation

clear; close all; clc

%Inputs:
f = @(t,u)(2*(u/(1+t)) - (t*u)^2);     %function handle for ODE
n = 1;                                 %dimension of system
tspan = [0,10];                        %time interval
u0 = 2;                                %intial conditions
N = 20;                                %no. of time sub-intervals steps
Ndg = 20;                              %no. of coarse steps (in each sub-interval)
Ndf = Ndg*100;                         %no. of fine steps (in each sub-interval)
sample_rule = 1;                       %sampling rule to employ
epsilon = 10^(-10);                    %error tolerance 
m = 30;                                %no. of samples at each sub-interval
sims = 1;                              %no. of independent sims

%solve with stochastic parareal
[t_fine,U,~,K,~,~] = stochasticparareal(f,n,tspan,u0,N,Ndg,Ndf,sample_rule,epsilon,m,sims);
U = U{1,1};  %selects solution from cell array

%solve using the fine solver serially
[~,u] = RK(t_fine,u0,f,'classic fourth-order');

%plot both solutions together
figure(1)
hold on
plot(t_fine,u,'k','LineWidth',2)              %fine solver solution
plot(t_fine(1:40:Ndf),U((1:40:Ndf),end),'ob') %stochastic parareal solution
xlabel('$t$','Interpreter','latex'); ylabel('$u_1(t)$','Interpreter','latex');
grid on; box on;
legend('northeast',{'$\mathcal{F}$ solution','$\mathcal{P}_s$ solution'},'Interpreter','latex')
hold off


%% Test Problem 2: Brusselator system

clear; close all; clc

%Inputs:
f = @(t,u)([1 + (u(1)^2)*u(2) - (3+1)*u(1); 3*u(1) - (u(1)^2)*u(2)]);     %function handle for ODE
n = 2;                                 %dimension of system
tspan = [0,15.3];                      %time interval
u0 = [1,3.07];                         %intial conditions
N = 25;                                %no. of time sub-intervals steps
Ndg = 25;                              %no. of coarse steps (in each sub-interval)
Ndf = Ndg*100;                         %no. of fine steps (in each sub-interval)
sample_rule = 1;                       %sampling rule to employ
epsilon = 10^(-6);                     %error tolerance 
m = 10;                                %no. of samples at each sub-interval
sims = 1;                              %no. of independent sims

%solve with stochastic parareal
[t_fine,U,~,K,~,~] = stochasticparareal(f,n,tspan,u0,N,Ndg,Ndf,sample_rule,epsilon,m,sims);
U = U{1,1};  %selects solution from cell array

%solve using the fine solver serially
[~,u] = RK(t_fine,u0,f,'classic fourth-order');

%plot both solutions together
figure(1)
hold on
plot(u(:,1),u(:,2),'k','LineWidth',2)
plot(U((1:50:Ndf+1),end-1),U((1:50:Ndf+1),end),'ob')
xlabel('$u_1$', 'Interpreter','latex'); ylabel('$u_2$', 'Interpreter','latex');
axis([0 4 0 5])
xticks((0:0.5:4))
yticks((0:0.5:5))
legend({'$\mathcal{F}$ solution','$\mathcal{P}_s$ solution'},'Interpreter','latex','location','northeast')
box on; grid on
hold off


%% Test Problem 3: Lorenz system

clear; close all; clc

%Inputs:
f = @(t,u)([10*(u(2) - u(1)); 28*u(1) - u(2) - u(1)*u(3); u(1)*u(2) - (8/3)*u(3)]);    %function handle for ODE
n = 3;                                 %dimension of system
tspan = [0,18];                        %time interval
u0 = [-15,-15,20];                     %intial conditions
N = 50;                                %no. of time sub-intervals steps
Ndg = 250;                             %no. of coarse steps (in each sub-interval)
Ndf = Ndg*75;                          %no. of fine steps (in each sub-interval)
sample_rule = 2;                       %sampling rule to employ
epsilon = 10^(-8);                     %error tolerance 
m = 10;                                %no. of samples at each sub-interval
sims = 1;                              %no. of independent sims

%solve with stochastic parareal
[t_fine,U,~,K,~,~] = stochasticparareal(f,n,tspan,u0,N,Ndg,Ndf,sample_rule,epsilon,m,sims);
U = U{1,1};  %selects solution from cell array

%solve using the fine solver serially
[~,u] = RK(t_fine,u0,f,'classic fourth-order');

%plot the solution in 3D phase space
figure(1)
hold on
plot3(u(:,1),u(:,2),u(:,3),'k','LineWidth',2)
plot3(U((1:(Ndf/Ndg):Ndf+1),end-2),U((1:(Ndf/Ndg):Ndf+1),end-1),U((1:(Ndf/Ndg):Ndf+1),end),'ob')
xlabel('$u_1$', 'Interpreter','latex'); 
ylabel('$u_2$', 'Interpreter','latex'); 
zlabel('$u_3$', 'Interpreter','latex');
view(45,10)
legend({'$\mathcal{F}$ solution','$\mathcal{P}_s$ solution'},'Interpreter','latex','location','northeast')
box on; grid on;
hold off


%% Test Problem 4: Nonlinear ODE

clear; close all; clc

%Inputs:
f = @(t,u)( sin(u)*cos(u) - 2*u + exp(-t/100)*sin(5*t) + log(1+t)*cos(t) );     %function handle for ODE
n = 1;                                 %dimension of system
tspan = [0,100];                        %time interval
u0 = 1;                                %intial conditions
N = 40;                                %no. of time sub-intervals steps
Ndg = 80;                              %no. of coarse steps (in each sub-interval)
Ndf = Ndg*100;                         %no. of fine steps (in each sub-interval)
sample_rule = 1;                       %sampling rule to employ
epsilon = 10^(-10);                    %error tolerance 
m = 2;                                %no. of samples at each sub-interval
sims = 1;                              %no. of independent sims

%solve with stochastic parareal
[t_fine,U,~,K,~,~] = stochasticparareal(f,n,tspan,u0,N,Ndg,Ndf,sample_rule,epsilon,m,sims);
U = U{1,1};  %selects solution from cell array

%solve using the fine solver serially
[~,u] = RK(t_fine,u0,f,'classic fourth-order');

%plot both solutions together
figure(1)
hold on
plot(t_fine,u,'k','LineWidth',2)
plot(t_fine(1:10:Ndf),U((1:10:Ndf),end),'ob')
xlabel('$t$','Interpreter','latex'); ylabel('$u_1(t)$','Interpreter','latex');
grid on; box on;
xlim([0 18]); ylim([-2 2]);
xticks((0:3:18))
legend({'$\mathcal{F}$ solution','$\mathcal{P}_s$ solution'},'Interpreter','latex','location','northwest')


%% Test Problem 5: Square limit cycle equations 

clear; close all; clc

%Inputs:
f = @(t,u)([-sin(u(1))*(0.1*cos(u(1)) + cos(u(2))); -sin(u(2))*(0.1*cos(u(2)) - cos(u(1)))  ]);     %function handle for ODE
n = 2;                                 %dimension of system
tspan = [0,60];                        %time interval
u0 = [1.5,1.5];                                %intial conditions
N = 30;                                %no. of time sub-intervals steps
Ndg = 30;                              %no. of coarse steps (in each sub-interval)
Ndf = Ndg*100;                         %no. of fine steps (in each sub-interval)
sample_rule = 2;                       %sampling rule to employ
epsilon = 10^(-8);                    %error tolerance 
m = 3;                                %no. of samples at each sub-interval
sims = 1;                              %no. of independent sims

%solve with stochastic parareal
[t_fine,U,~,K,~,~] = stochasticparareal(f,n,tspan,u0,N,Ndg,Ndf,sample_rule,epsilon,m,sims);
U = U{1,1};  %selects solution from cell array

%solve using the fine solver serially
[~,u] = RK(t_fine,u0,f,'classic fourth-order');

%plot both solutions together
figure(1)
hold on
plot(u(:,1),u(:,2),'k','LineWidth',2)
plot(U((1:60:Ndf+1),end-1),U((1:60:Ndf+1),end),'ob')
axis([-0.75 pi+0.75 -0.75 pi+0.75])
xlabel('$u_1$', 'Interpreter','latex'); ylabel('$u_2$', 'Interpreter','latex');
xticks((-0.5:0.5:4.5))
yticks((-0.5:0.5:4.5))
legend({'$\mathcal{F}$ solution','$\mathcal{P}_s$ solution'},'Interpreter','latex','location','northeast')
box on; grid on
hold off
