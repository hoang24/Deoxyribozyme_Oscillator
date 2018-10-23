% deoxyribozyme oscillator parameters
h = 0.7849; % fraction of the reactor chamber that is well-mixed
b = 5e-7; % reaction rate constant for the gate-substrate reaction (nM/s)
e = 8.8750e-11; % efflux rate (nL/s)
V = 7.54e-9; % volume of the reactor (nL)
Sm_base = 5.45e-6; % base value for substrate-influx rate (nmol/s)

% initial conditions
load('ODE_non-perturb.mat');

% initial concentration of substrates S1, S2, S3, respectively (nM)
S1o = summaryTable(end,2); 
S2o = summaryTable(end,3); 
S3o = summaryTable(end,4); 

% initial concentration of products P1, P2, P3, respectively (nM)
P1o = summaryTable(end,5); 
P2o = summaryTable(end,6); 
P3o = summaryTable(end,7); 

G = 2500; % uniform gate concentration (nM)
G1 = G; G2 = G; G3 = G;

S1 = S1o; S2 = S2o; S3 = S3o; P1 = P1o; P2 = P2o; P3 = P3o;

% influx rates
Sm3 = Sm_base;

% timing parameters
ts = 500; % start time 
T = 1000; % maximum elapsed time
thold = 50;
period = (T - ts) / thold;

t = ts;
for i = 1 : period
    Sm1 = Sm_base * rand; Sm2 = Sm_base * rand; 
    tspan = [ts + thold*(i-1), ts + thold*i];
    conds = [S1(end), S2(end), S3(end), P1(end), P2(end), P3(end)];
    % concentrations: S1 = x(1), S2 = x(2), S3 = x(3), P1 = x(4), P2 = x(5), P3 = x(6); 
    eqns = @(t,x) [ (Sm1   /V) - h*b*x(1)*(G1 - x(6)) - (e/V)*x(1); ...
                    (Sm2   /V) - h*b*x(2)*(G2 - x(4)) - (e/V)*x(2); ...
                    (Sm3   /V) - h*b*x(3)*(G3 - x(5)) - (e/V)*x(3); ...
                                 h*b*x(1)*(G1 - x(6)) - (e/V)*x(4); ...
                                 h*b*x(2)*(G2 - x(4)) - (e/V)*x(5); ...
                                 h*b*x(3)*(G3 - x(5)) - (e/V)*x(6)];
    [tPerturb,x] = ode23t(eqns, tspan, conds);
      
    tPerturb = tPerturb(2:end,:);
    x = x(2:end,:);
    t = [t; tPerturb];
    S1 = [S1; x(:,1)]; S2 = [S2; x(:,2)]; S3 = [S3; x(:,3)]; % concentration of substrates over time
    P1 = [P1; x(:,4)]; P2 = [P2; x(:,5)]; P3 = [P3; x(:,6)]; % concentration of products over time 
end
summaryTable = [t, S1, S2, S3, P1, P2, P3];