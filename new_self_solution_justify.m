clc;
clear;              
m = 0.013;
zeta = 0.05;
omega_n = 2 * pi * 778;
M = m * eye(2);
C = 2 * m * zeta * omega_n * eye(2);
K = m * omega_n^2 * eye(2);
l = 2;
mega_omega =5000;
dt = 1/10/mega_omega;         
t1 = 0;
t2 = 0.5;
t = t1:dt:t2;      
x = zeros(2,length(t));
y = zeros(2,length(t));
u = zeros(2,length(t));
HXX = zeros(1,length(t));
HXY = zeros(1,length(t));
HYX = zeros(1,length(t));
HYY = zeros(1,length(t));
G = zeros(2+4*l,length(t));
mu_x = zeros(2+4*l,length(t));
mu_y = zeros(2+4*l,length(t));
F = zeros(2,length(t));
x(:,1) = 0;        
b = 0.00095;
N=2;
tau = 60 / (N * mega_omega);
a = 0.0025;
R = 0.0025;
y(:,1) = 2*pi*R*mega_omega;  
phi_e = 0;
phi_a = acos(1 - a / R);
Kt = 6 * 10^8;
Kn = 2 * 10^8;
lambda_x = eye(2);
lambda_y = eye(2);
omega = 2*pi*N*mega_omega/60;


for i = 1:tau/dt
    H = Hl(i*dt,mega_omega,N,phi_e,phi_a,Kt,Kn);
    HXX(1,i) = H(1,1);
    HXY(1,i) = H(1,2);
    HYX(1,i) = H(2,1);
    HYY(1,i) = H(2,2);
    F(:,i) = b*H*x(:,i);
    [M1,N1,K1,L1] = mnkl1(M,C,K,F(:,i),u(:,i),omega,i*dt,l,x(1,i),x(2,i),y(1,i),y(2,i));
    [M2,N2,K2,L2] = mnkl1(M,C,K,F(:,i),u(:,i),omega,(i+0.5)*dt,l,x(1,i)+0.5*dt*K1(1,1),x(2,i)+0.5*dt*K1(2,1),y(1,i)+0.5*dt*L1(1,1),y(2,i)+0.5*dt*L1(2,1));
    [M3,N3,K3,L3] = mnkl1(M,C,K,F(:,i),u(:,i),omega,(i+0.5)*dt,l,x(1,i)+0.5*dt*K2(1,1),x(2,i)+0.5*dt*K2(2,1),y(1,i)+0.5*dt*L2(1,1),y(2,i)+0.5*dt*L2(2,1));
    [M4,N4,K4,L4] = mnkl1(M,C,K,F(:,i),u(:,i),omega,(i+1)*dt,l,x(1,i)+dt*K3(1,1),x(2,i)+dt*K3(2,1),y(1,i)+dt*L3(1,1),y(2,i)+dt*L3(2,1));
    mu_x(:,i+1) = mu_x(:,i) + 1/6*dt*(M1+2*M2+2*M3+M4);
    mu_y(:,i+1) = mu_y(:,i) + 1/6*dt*(N1+2*N2+2*N3+N4);
    u(:,i+1) = -[mu_x(:,i+1).'* g1(i*dt,omega,l,x(1,i),x(2,i));mu_y(:,i+1).'* g1(i*dt,omega,l,x(1,i),x(2,i))];
    G(:,i) = g1(i*dt,omega,l,x(1,i),x(2,i));
    x(:,i+1) = x(:,i) + 1/6*dt*(K1+2*K2+2*K3+K4);
    y(:,i+1) = y(:,i) + 1/6*dt*(L1+2*L2+2*L3+L4);
end

for i = tau/dt+1:length(t)-1
    H = Hl(i*dt,mega_omega,N,phi_e,phi_a,Kt,Kn);
    HXX(1,i) = H(1,1);
    HXY(1,i) = H(1,2);
    HYX(1,i) = H(2,1);
    HYY(1,i) = H(2,2);
    F(:,i) = b*H*(x(:,i-tau/dt)-x(:,i)); 
    [M1,N1,K1,L1] = mnkl2(M,C,K,F(:,i),u(:,i),omega,i*dt,l,x(1,i-tau/dt),x(1,i),x(2,i-tau/dt),x(2,i),y(1,i),y(2,i));
    [M2,N2,K2,L2] = mnkl2(M,C,K,F(:,i),u(:,i),omega,(i+0.5)*dt,l,x(1,i-tau/dt),x(1,i)+0.5*dt*K1(1,1),x(2,i-tau/dt),x(2,i)+0.5*dt*K1(2,1),y(1,i)+0.5*dt*L1(1,1),y(2,i)+0.5*dt*L1(2,1));
    [M3,N3,K3,L3] = mnkl2(M,C,K,F(:,i),u(:,i),omega,(i+0.5)*dt,l,x(1,i-tau/dt),x(1,i)+0.5*dt*K2(1,1),x(2,i-tau/dt),x(2,i)+0.5*dt*K2(2,1),y(1,i)+0.5*dt*L2(1,1),y(2,i)+0.5*dt*L2(2,1));
    [M4,N4,K4,L4] = mnkl2(M,C,K,F(:,i),u(:,i),omega,(i+1)*dt,l,x(1,i-tau/dt),x(1,i)+dt*K3(1,1),x(2,i-tau/dt),x(2,i)+dt*K3(2,1),y(1,i)+dt*L3(1,1),y(2,i)+dt*L3(2,1));
    mu_x(:,i+1) = mu_x(:,i) + 1/6*dt*(M1+2*M2+2*M3+M4);
    mu_y(:,i+1) = mu_y(:,i) + 1/6*dt*(N1+2*N2+2*N3+N4);
    u(:,i+1) = -[mu_x(:,i+1).'* g2(i*dt,omega,l,x(1,i-tau/dt),x(1,i),x(2,i-tau/dt),x(2,i));mu_y(:,i+1).'* g2(i*dt,omega,l,x(1,i-tau/dt),x(1,i),x(2,i-tau/dt),x(2,i))];
    G(:,i) = g2(i*dt,omega,l,x(1,i-tau/dt),x(1,i),x(2,i-tau/dt),x(2,i));
    x(:,i+1) = x(:,i) + 1/6*dt*(K1+2*K2+2*K3+K4);
    y(:,i+1) = y(:,i) + 1/6*dt*(L1+2*L2+2*L3+L4);
end
subplot(2,2,1)
plot(t,x)
xlabel('Time t');
ylabel('Displacement x');
subplot(2,2,2)
plot(t,y)
xlabel('Time t');
ylabel('Velocity x');
subplot(2,2,3)
plot(t,F)
xlabel('Time t')
ylabel('Force F')
subplot(2,2,4)
plot(t,u)
xlabel('Time t')
ylabel('Force u')

figure
subplot(2,2,1)
plot(t,HXX)
xlabel('Time t')
ylabel('Hxx')
subplot(2,2,2)
plot(t,HXY)
xlabel('Time t')
ylabel('Hxy')
subplot(2,2,3)
plot(t,HYX)
xlabel('Time t')
ylabel('Hyx')
subplot(2,2,4)
plot(t,HYY)
xlabel('Time t')
ylabel('Hyy')
figure
plot(t,G)
xlabel('Time t')
ylabel('g')

