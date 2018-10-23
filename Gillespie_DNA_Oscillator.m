% Reactions of the DNA enzyme oscillator CRN
% (1)  S1  --(hb(G1-P3))-->  P1
% (2)  S2  --(hb(G2-P1))-->  P2
% (3)  S3  --(hb(G3-P2))-->  P3
% (4)  0   --(Sm1/V)-->      S1
% (5)  0   --(Sm2/V)-->      S2
% (6)  0   --(Sm3/V)-->      S3
% (7)  S1  --(e/V)-->        0
% (8)  S2  --(e/V)-->        0
% (9)  S3  --(e/V)-->        0
% (10) P1  --(e/V)-->        0
% (11) P2  --(e/V)-->        0
% (12) P3  --(e/V)-->        0

% update 1 nM in a Gillespie iteration
% V = 7.54 nL; CM = 1 nM 
% --> #mol = CM * V = 7.54 amol 
% --> #molecule = NA * #mol = 4540588 molecules
% --> update 4 540 588 molecules in a Gillespie iteration

clear
close all
clc

fileName = datestr(datetime('now'));
tStart = tic;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
h = 0.7849; % fraction of the reactor chamber that is well-mixed
b = 5e-7; % reaction rate constant for the gate-substrate reaction (nM/s)
e = 8.8750e-11; % efflux rate (L/s)
V = 7.54e-9; % volume of the reactor (L)
Sm_base = 5.45e-6; % base value for substrate-influx rate (nmol/s)

S1 = 0; S2 = 0; S3 = 0; % initial concentration of substrates S1, S2, S3, respectively (nM)
P1 = 1000; P2 = 0; P3 = 0; % initial concentration of products P1, P2, P3, respectively (nM)
G = 2500; % uniform gate concentration (nM)
G1 = G; G2 = G; G3 = G;

t = 0; % start time 
T = 500; % maximum elapsed time

maxSize = 2400 * T; % 2400: calibrated value
summaryTable = zeros(maxSize, 7); 
summaryTable(1, 1:7) = [t, S1, S2, S3, P1, P2, P3]; % table summarizing concentration of different species at different time

numIteration = 0; % number of iteration

dCon = 1; % the change in concentration in every iteration

P = zeros(maxSize, 12); 

% Run Gillespie algorithm
while (t < T)
   numIteration = numIteration + 1;
   % random influx rate for S1 and S2 and S3
   Sm1 = Sm_base * 1;
   Sm2 = Sm_base * 1;
   Sm3 = Sm_base * 1;
    
   % Reaction rate for each reaction
   r1 =   h * b * S1 * (G1 - P3);
   r2 =   h * b * S2 * (G2 - P1);
   r3 =   h * b * S3 * (G3 - P2);
   r4 = Sm1 / V;
   r5 = Sm2 / V;
   r6 = Sm3 / V;
   r7 = e * S1 / V;
   r8 = e * S2 / V;
   r9 = e * S3 / V;
   r10 = e * P1 / V;
   r11 = e * P2 / V;
   r12 = e * P3 / V;
   R = [ r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12 ];
   
   % Total reaction rate
   rTot = sum(R);
   
   % Probability for each reaction to happen
   %P = zeros(1, length(R));
   for i = 1 : length(R)
       P(numIteration, i) = R(i) / rTot;
   end  
   
   % increase time step by dt
   mu = 1 / rTot; % mean of the exponential distribution
   dt = exprnd(mu);
   t = t + dt;
   
   % update species concentration
   x = sum(rand >= cumsum([0, P(numIteration,:)]));
   if (x == 1)
       S1 = S1 - dCon;
       P1 = P1 + dCon;
   elseif (x == 2)
       S2 = S2 - dCon;
       P2 = P2 + dCon;
   elseif (x == 3) 
       S3 = S3 - dCon;
       P3 = P3 + dCon;
   elseif (x == 4)
       S1 = S1 + dCon;
   elseif (x == 5)
       S2 = S2 + dCon;
   elseif (x == 6)
       S3 = S3 + dCon;
   elseif (x == 7)
       S1 = S1 - dCon;
   elseif (x == 8)
       S2 = S2 - dCon;
   elseif (x == 9)
       S3 = S3 - dCon;
   elseif (x == 10)
       P1 = P1 - dCon;
   elseif (x == 11)
       P2 = P2 - dCon;
   elseif (x == 12)
       P3 = P3 - dCon;
   end
   
   %summaryTable = [summaryTable ; t, S1, S2, S3, P1, P2, P3];
   summaryTable(numIteration + 1, 1:7) = [t, S1, S2, S3, P1, P2, P3];
end

summaryTable = summaryTable(1:(numIteration + 1), :);
P = P(1:numIteration, :);
tElapsed = toc(tStart);

save(fileName, 'summaryTable', 'numIteration', 'tElapsed');
