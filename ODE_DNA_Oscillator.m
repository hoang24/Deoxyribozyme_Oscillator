% deoxyribozyme oscillator parameters
h = 0.7849; % fraction of the reactor chamber that is well-mixed
b = 5e-7; % reaction rate constant for the gate-substrate reaction (nM/s)
e = 8.8750e-11; % efflux rate (nL/s)
V = 7.54e-9; % volume of the reactor (nL)
Sm_base = 5.45e-6; % base value for substrate-influx rate (nmol/s)

% initial conditions
S1o = 0; S2o = 0; S3o = 0; % initial concentration of substrates S1, S2, S3, respectively (nM)
P1o = 1000; P2o = 0; P3o = 0; % initial concentration of products P1, P2, P3, respectively (nM)
G = 2500; % uniform gate concentration (nM)
G1 = G; G2 = G; G3 = G;

% influx rates
Sm1 = Sm_base; Sm2 = Sm_base; Sm3 = Sm_base;

% timing parameters
ts = 0; % start time 
T = 500; % maximum elapsed time

conds = [S1o, S2o, S3o, P1o, P2o, P3o];
tspan = [ts, T];

% concentrations: S1 = x(1), S2 = x(2), S3 = x(3), P1 = x(4), P2 = x(5), P3 = x(6); 
eqns = @(t,x) [ (Sm1   /V) - h*b*x(1)*(G1 - x(6)) - (e/V)*x(1); ...
                         (Sm2   /V) - h*b*x(2)*(G2 - x(4)) - (e/V)*x(2); ...
                         (Sm3   /V) - h*b*x(3)*(G3 - x(5)) - (e/V)*x(3); ...
                                      h*b*x(1)*(G1 - x(6)) - (e/V)*x(4); ...
                                      h*b*x(2)*(G2 - x(4)) - (e/V)*x(5); ...
                                      h*b*x(3)*(G3 - x(5)) - (e/V)*x(6)];
[t,x] = ode23t(eqns, tspan, conds);
S1 = x(:,1); S2 = x(:,2); S3 = x(:,3); % concentration of substrates over time
P1 = x(:,4); P2 = x(:,5); P3 = x(:,6); % concentration of products over time

summaryTable = [t, S1, S2, S3, P1, P2, P3];
