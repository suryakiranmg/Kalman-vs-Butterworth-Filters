clc; clear all;close all;
%Given Parameters for Digital Filter Implementation
fs = 100;
fc = 5;
%Given Data
load noisy.data
yn = noisy(:,1); 
xn = noisy(:,2); %True Data
un = noisy(:,3);
Length = 1000;

%%%%%%####-----Digital Filter Implementation------######%%%
% Butterworth Coefts
[b,a] = butter(2,fc/(fs/2),'low')
% Digital Filter
filtered_data = filter(b,a,yn);
t = 1:1:Length;
% Plots Digital Filter Output Vs True Data
figure('Name','Butterworth - Filtered Data vs True Data','NumberTitle','off');
hold on;grid on;
xlabel('seq no.'); ylabel('magnitude')
plot(t,filtered_data,'-r','LineWidth',2);
plot(t,xn,'-b','LineWidth',2);
legend('Filtered Data','True Data');
hold off; 
%Delay from digital filter
delay_timesteps = finddelay(xn,filtered_data);
digital_filter_delay_in_secs = delay_timesteps * 0.01

%%%%%%-----Kalman Filter Implementation------%%%
Q = 0.01;
R = 1;
A = 0.9; B = 0.5; C = 1;

X(1) = 0;
P(1) = 1;  % Initial choice for error covarance
%Memory Allocation
K_plot = zeros(1, Length);
P_plot = zeros(1, Length);
P_plot(1) = 1;
%Discrete Kalman Filter Equations
for n = 2:Length 
X(n) = A * X(n-1) + B * un(n-1);
P_est = ( A* P(n-1) *A' ) + Q;
K = P_est / (P_est + R);
X(n) = X(n) + K * (yn(n) - X(n)) ;
P(n) = ( 1 - K ) * P_est ;
% K and P for plots
K_plot(n) = K ;
P_plot(n) = P(n);
end
%Plots for performance comparison
figure('Name','Kalman - Filtered Data vs True Data','NumberTitle','off');
hold on;grid on;
xlabel('seq no.'); ylabel('magnitude')
plot(t,X,'-r','LineWidth',2);
plot(t,xn,'-b','LineWidth',2);
legend('Kalman Filtered Data','True Data');
hold off;
%Filter delay
delay_timesteps = finddelay(xn,X);
Kalman_filter_delay_in_secs = delay_timesteps * 0.01
%Kalman Gain and Error Plots
figure('Name','Kalman Gain (K) vs TimeSteps','NumberTitle','off');
plot(t,K_plot,'-r','LineWidth',2);grid on;
figure('Name','Posterior Error Covariance Matrix (P) vs TimeSteps','NumberTitle','off');
plot(t,P_plot,'-b','LineWidth',2);grid on;


