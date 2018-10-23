% readout driver
tic
plotData_2(summaryTable);
iterations = 100000;
coeff = 0.001;
bit = 10; % change me
[HamVec, HamDist] = HammingDistance(inputs);
[actualOut, error, weights] = readOut(summaryTable, HamVec, bit, iterations, coeff);
toc