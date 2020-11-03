function deriv = dydt_ode45(time,statevar)

% Assignment4: Hodgkin-Huxley Model
% Solved using MATLAB ode45
% Roshni (BE17B009)

global Cm I ENa EK El gbarNa gbarK gbarl 

%Initial conditions
V = statevar(1); m = statevar(2); h = statevar(3); n = statevar(4); 

gNa = gbarNa*m^3*h;
gK = gbarK*n^4;
gl = gbarl;

INa = gNa*(V-ENa);
IK = gK*(V-EK);
Il = gl*(V-El);

%ODE's
dm = am(V)*(1-m) - bmm(V)*m ;
dh = ah(V)*(1-h) - bh(V)*h ;
dn = an(V)*(1-n) - bn(V)*n ; 
dV = (I-(INa+IK+Il))/Cm;

deriv = [dV;dm;dh;dn];

return
