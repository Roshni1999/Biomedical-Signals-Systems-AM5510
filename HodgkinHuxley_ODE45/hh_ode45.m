% Assignment4: Hodgkin-Huxley Model
% Solved using MATLAB ode45
% Roshni (BE17B009)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Units
% t                   time                    ms
% V                   membrane potantial      mV
% INa,IK,Il           ionic current           uA/cm^2
% Cm                  capacitance             uF/cm^2
% gNa, gK, gl         conductance             uS/cm^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants specified
global Cm I ENa EK El gbarNa gbarK gbarl 
Cm = 0.01;       % Membrane Capcitance uF/cm^2
dt = 0.04;       % Time Step ms
t = 0:dt:25;     % Time Array ms
I = 0.2;         % External Current Applied
ENa = 55.17;     % mv Na reversal potential
EK = -72.14;     % mv K reversal potential
El = -49.42;     % mv Leakage reversal potential
gbarNa = 1.2;    % uS/cm^2 Na conductance
gbarK = 0.36;    % uS/cm^2 K conductance
gbarl = 0.003;   % uS/cm^2 Leakage conductance

V0 = -60;                       % Initial Membrane voltage
m0 = am(V0)/(am(V0)+bmm(V0));   % Initial m-value
n0 = an(V0)/(an(V0)+bn(V0));    % Initial n-value
h0 = ah(V0)/(ah(V0)+bh(V0));    % Initial h-value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain V, m, n, h by solving ODE's
statevar_i = [V0, m0, h0, n0];            % Initial Conditions
[time,statevars] = ode45(@dydt_ode45, t ,statevar_i) ; 

% Values obtained by solving ODE's
V = statevars(:, 1); 
m = statevars(:, 2); 
h = statevars(:, 3); 
n = statevars(:, 4); 

gNa = gbarNa*m.^3.*h;
gK = gbarK*n.^4;
gl = gbarl;

INa = gNa.*(V-ENa);
IK = gK.*(V-EK);
Il = gl.*(V-El);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results

figure(1)
plot(time, V)
xlabel('time (ms)')
ylabel('V_m (mV)')

figure(2)
hold on
plot(t,INa,'b')
plot(t,IK,'r')
plot(t,Il,'g')
legend('I_N_a','I_K','I_L')
xlabel('time (ms)')
ylabel('Ionic current')
hold off

figure(3)
hold on
plot(t,gNa,'b')
plot(t,gK,'r')
legend('g_N_a','g_K')
xlabel('time (ms)')
ylabel('Conductances')

figure(4)
hold on
plot(t,m,'b')
plot(t,h,'r')
plot(t,n,'g')
legend('m','h','n')
xlabel('time (ms)')
ylabel('Gating variables')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




