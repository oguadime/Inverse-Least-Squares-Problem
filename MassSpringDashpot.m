function lab1
% (used with permission from: Veronika Vasylkivska 11/2008)

% set up the values of m,c,k
global m c k
m = 2;
c = 0.5;
k = 1.5;

% set up the parameters of integration
global T x0 v0
T = 20;
x0 = 10;
v0 = 25;

% parameters of real solution
nu = sqrt(4*k*m - c^2)/(2*m);
A = x0;
B = (v0 + c*x0/(2*m))/nu;
r = c/(2*m);

X = @ (t) exp(-r*t).*(A*cos(nu*t) + B*sin(nu*t));

%% (a) solve the differential equation and plot it together 
% with analytic solution

K = k/m;
C = c/m;
disp(['True Parameters: K = ' num2str(k/m) ', C = ' num2str(c/m)]);
[te,Xe] = ode45(@ (t,y) dy(t,y,K,C),[0 T],[x0 v0]);
%Determine maximum time step used by integrator:
dt=te(2:end)-te(1:end-1);
dtmax=max(dt);

figure
plot(te,Xe(:,1),'-o',te,X(te));
xlabel('Time, t');
ylabel('x(t)');
set(gca, 'LineWidth', 1.2);
legend('x_h(t)','x(t)');
title(['Numerical and analytic solutions: m = ' num2str(m) ', c = ' ...
       num2str(c) ', k = ' num2str(k), '(max time step = ', num2str(dtmax),')']);

% plot a difference of two solutions
figure
plot(te,(Xe(:,1)-X(te)),'-o');
xlabel('Time, t');
ylabel('x_h(t) - x(t)');
set(gca, 'LineWidth', 1.2);
title(['Numerical residual using: m = ' num2str(m) ', c = ' num2str(c) ', k = ' num2str(k)]);



%% (b) generate the "simulated" data using the analytic solution
M = 100;
global h Sim_data
h = T/(M - 1);
t = 0:h:T; % equally-spaced time steps
Sim_data = X(t)';
disp([t' Sim_data]);


%% (c) compare a numerical solution method to data using the least squares objective functional

[t,Xh] = ode45(@ (t,y) dy(t,y,K,C),t,[x0 v0]);
disp(['LSOF for the numerical solution with true parameters is ' num2str(LSOF(Xh(:,1),Sim_data))]);

%% (d)
disp(' ')
disp('Optimization using lsqnonlin');

% array D will keep all necessary information
D = zeros(6,6);

% set up the "initial" guesses for C and K
D(:,1) = [K; K*1.1; K*0.9; K*2.3; K*5; K*20];
D(:,2) = [C; C*0.89; C*1.3; C*4; C*5; C*12];

options=optimset('lsqnonlin');
newopts=optimset(options,'Display','off');

% Perform optimization for different initial guesses
for j=1:size(D,1)
    D(j,3) = 0.5*sum(myfun(D(j,1:2)).^2);
    [D(j,4:5),resnorm,residual,exitflag,output] = lsqnonlin(@myfun,D(j,1:2),[],[],newopts);
    D(j,6) = 0.5*resnorm;
    D(j,7) = output.funcCount;
end;
format short g
disp('  InitGuess K  InitGuess C     InitCost     FinEst K     FinEst C      FinCost     FunEvals');
disp(D);

%% (e)-(f)
% choose one of the initial estimates randomly
row = 4;
disp(['Estimates with noise (initial estimate: K = ' num2str(D(row,1)) ', C = ' num2str(D(row,2)),')']);
P = zeros(6,6);
P(:,1) = [0; 0.01; .02; 0.05; 0.1; 0.2]*10; % setting noise levels
S = Sim_data;
res=zeros(6,1);

for j = 1:size(P,1)
    Sim_data = S + P(j,1)*randn(M,1); % add noise to data
    P(j,2) = 0.5*sum(myfun(D(row,1:2)).^2);
    [P(j,3:4),res(j),residual,exitflag,output] = lsqnonlin(@myfun,D(row,1:2),[],[],newopts);
    P(j,5) = 0.5*res(j); 
    P(j,6) = output.funcCount;
    pause(0.1);
end;

% displaying the data on the estimation of parameters with noised data
disp('   NoiseLevel     InitCost     FinEst K     FinEst C      FinCost     FunEvals');    
disp(P);

figure;
plot(P(:,1),P(:,5),'-o');
xlabel('Noise Level');
ylabel('min LSOF');
set(gca, 'LineWidth', 1.2);
title('The Minimum of the Objective Function Value');

figure;
loglog(P(2:end,1),P(2:end,5),'-o');
axis('tight')
xlabel('Noise Level');
ylabel('min LSOF');
set(gca, 'LineWidth', 1.2);
title('The Minimum of the Objective Function Value');

%% (g)
% first we need to build a matrix e(q_0)
st = 0.01;
epsilon = zeros(M,2);
interval = zeros(size(P,1),4);
leng = zeros(size(P,1),2);

for j = 1:size(P,1)
    % find the appropriate solutions to approximate the partial derivatives
    % with respect to K and C
    [t,X_p_K] = ode45(@ (t,y) dy(t,y,P(j,3)+st/2,P(j,4)),t,[x0 v0]);
    [t,X_m_K] = ode45(@ (t,y) dy(t,y,P(j,3)-st/2,P(j,4)),t,[x0 v0]);
    [t,X_p_C] = ode45(@ (t,y) dy(t,y,P(j,3),P(j,4)+st/2),t,[x0 v0]);
    [t,X_m_C] = ode45(@ (t,y) dy(t,y,P(j,3),P(j,4)-st/2),t,[x0 v0]);
    
    [t,Xh] = ode45(@ (t,y) dy(t,y,P(j,3),P(j,4)),t,[x0 v0]);
    
    % calculate the sample average and variance
    sigmaS=res(j)/(M-2);

    % building the covariance matrix
    epsilon(:,1) = (X_p_K(:,1) - X_m_K(:,1))/st;
    epsilon(:,2) = (X_p_C(:,1) - X_m_C(:,1))/st;
    CovMatrix = sigmaS*inv(epsilon'*epsilon);
    
    % sigma for K and C
    sigma_K = sqrt(abs(CovMatrix(1,1)));
    sigma_C = sqrt(abs(CovMatrix(2,2)));
    
    % compute the lengths of the confidence intervals and intervals themself
    leng(j,1) = tinv(0.975,M-2)*sigma_K;
    interval(j,1) = P(j,3) - leng(j,1);
    interval(j,2) = P(j,3) + leng(j,1);
    
    leng(j,2) = tinv(0.975,M-2)*sigma_C;
    interval(j,3) = P(j,4) - leng(j,2);
    interval(j,4) = P(j,4) + leng(j,2);
end;

% displaying the confidence intervals
disp('Confidence Intervals');
disp('   NoiseLevel        Est K        Est C     LLimit K     HLimit K     LLimit C     HLimit C');
disp([P(:,[1 3 4]),interval]);

% plotting the estimates for K and C together with confidence intervals

figure;
errorbar(P(:,1),P(:,3),leng(:,1),'-o');
xlabel('Noise Level');
ylabel('K');
set(gca, 'LineWidth', 1.2);
title('The Confidence Interval for K');

figure;
errorbar(P(:,1),P(:,4),leng(:,2),'-o');
xlabel('Noise Level');
ylabel('C');
set(gca, 'LineWidth', 1.2);
title('The Confidence Interval for C');

%% section of subfunctions

function F = dy(t,y,K,C)
    F = zeros(2,1);   
    F(1) = y(2);
    F(2) = -K*y(1) - C*y(2);
    
function x = myfun(q)
    global h T Sim_data x0 v0
    [t,Xq] = ode45(@ (t,y) dy(t,y,q(1,1),q(1,2)),0:h:T,[x0 v0]);
    x = Xq(:,1) - Sim_data;
    
function value = LSOF(Xh,X)
    value = .5*sum((Xh - X).^2);
