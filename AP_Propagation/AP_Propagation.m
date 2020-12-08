% Mathematical Modeling of Action Potential with Transmission Equations and Hodgkin-Huxley Model
% BENG 221 Problem Solving Report

dt=0.0001;  %Time-step
dx=1;       %Unit step length along axon
s=dt/dx;
t=0:dt:10;  %Total time- 10ms
x=0:dx:10;  %Total length- 10units

% Initial values of rate constants of gating variables- alpha and beta
alpha_n0=0.01*(0+10)/(-1+exp((0+10)/10));
beta_n0=0.125*exp(0/80);
alpha_m0=0.1*(0+25)/(-1+exp((0+25)/10));
beta_m0=4*exp(0/18);
alpha_h0=0.07*exp(0/20); 
beta_h0=(1+exp((0+30)/10))^-1;

% Empty matrices to store the values of n, m, h, v, phi(Vx), psi(Vt)
% phi(Vx)= dV/dx, psi(Vt)= dV/dt
% Values are stored across length of axon(rows) and time(columns)
n=zeros(length(x),length(t));
m=zeros(length(x),length(t));
h=zeros(length(x),length(t));
v=zeros(length(x),length(t));
psi=zeros(length(x),length(t));
phi=zeros(length(x),length(t));

% Values of gating variables n, m, h at t=0
n(:,1)=alpha_n0/(alpha_n0+beta_n0);
m(:,1)=alpha_m0/(alpha_m0+beta_m0);
h(:,1)=alpha_h0/(alpha_h0+beta_h0);

%Setting initial stimulus voltage

%Impulse of -15mV can cross the threshold and is sufficient to produce an
%action potential
%-15mV impulse

for i=1:1:30000
    v(1,i)= -15;   %Impulse at x=0 across all time points = -15mV
end


%{
%Below threshold - unable to produce an action potential
for i=1:1:30000
    v(1,i)= -5;  %Impulse at x=0 across all time points = -5mV
end
%}
for i=1:1:length(x)-1
    phi(i,1)=1.5;  % phi at t=0 along axon
    v(i+1,1)=v(i,1)+dx*phi(i,1); % V at t=0 along axon
end

%Finite differece equations to solve the dfferential equations
for j=1:1:length(t)-1
    for i=1:1:length(x)-1
        n(i,j+1)=n(i,j) + dt*N(n(i,j),v(i,j));
        m(i,j+1)=m(i,j) + dt*M(m(i,j),v(i,j));
        h(i,j+1)=h(i,j) + dt*H(h(i,j),v(i,j));
        v(i,j+1)=v(i,j) + dt*psi(i,j);
        phi(i+1,j)=phi(i,j) + s*(psi(i+1,j)-psi(i,j));
        psi(i,j+1)=psi(i,j) + s*(phi(i+1,j)-phi(i,j)) - dt*F(psi(i,j),v(i,j),n(i,j),m(i,j),h(i,j));
     end
end

% ?? Is the last difference equation for psi dimensionally correct? 
% ?? ‘s’ term should be ‘1/s’ since ? is (dV/dx) and s has been defined as d?/dx.
% ?? Introduce a smaller step-size length along the axon so that 
% 1/s= dx/d?; will be a small.

figure(1)
plot(t,v(1,:))
set(gca, 'YDir','reverse')

figure(2)
plot(t,v(1,:),t,v(2,:),t,v(3,:),t,v(4,:),t,v(5,:),t,v(6,:),t,v(7,:),t,v(8,:),t,v(9,:),t,v(10,:))
legend('1','2','3','4','5','6','7','8','9','10')
set(gca, 'YDir','reverse')
