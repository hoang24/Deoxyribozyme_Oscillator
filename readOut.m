function [actualOut, error, weights] = readOut(summaryTable, target, bit, iterations, coeff)
    % This function implements the readout layer into the deoxyribozyme
    % oscillator CRN with Gillespie algorithm
    % Parameters:
        % summaryTable = table of Gillespie time and concentrations
        % tPerturb_vec = vector of the perturbed time
        % desiredOut = Hamming vector
        % bit = number of bits/ number of perturbations
    % Returns:
    
    % decompress table of concentrations into analog signals
    t  = summaryTable(:, 1);
    S1 = summaryTable(:, 2);
    S2 = summaryTable(:, 3);
    S3 = summaryTable(:, 4);
    P1 = summaryTable(:, 5);
    P2 = summaryTable(:, 6);
    P3 = summaryTable(:, 7);
    
    % Get random indices of the species
    iS1 = randi(length(S1), [bit, 1]);
    iS2 = randi(length(S2), [bit, 1]);
    iS3 = randi(length(S3), [bit, 1]);
    %
    iP1 = randi(length(P1), [bit, 1]);
    iP2 = randi(length(P2), [bit, 1]);
    iP3 = randi(length(P3), [bit, 1]);
    
    % Vector of concentration of the substrates and products at specific times (analog)
    vS1 = S1(iS1);
    vS2 = S2(iS2);
    vS3 = S3(iS3);
    %
    vP1 = P1(iP1);
    vP2 = P2(iP2);
    vP3 = P3(iP3);
    
    % Digital input vector of scaled concentration of the substrates and products at specific times
    inS1 = vS1 / max(S1);
    inS2 = vS2 / max(S2);
    inS3 = vS3 / max(S3);
    %
    inP1 = vP1 / max(P1);
    inP2 = vP2 / max(P2);
    inP3 = vP3 / max(P3);
    
    % Bias
    bias = 1;
    
    % Weights
    wBias = rand;
    wS1 = rand;
    wS2 = rand;
    wS3 = rand;
    wP1 = rand;
    wP2 = rand;
    wP3 = rand;
    
    % Reallocation
    actualOut = zeros(bit, iterations);
    error = zeros(bit, iterations);
    weights = [wS1, wS2, wS3, wP1, wP2, wP3, wBias];
    
    for i = 1:iterations
        for j = 1:bit
            net = inS1(j)*wS1 + inS2(j)*wS2 + inS3(j)*wS3 ...
                + inP1(j)*wP1 + inP2(j)*wP2 + inP3(j)*wP3 + bias*wBias;
            actualOut(j,i) = 1./(1+exp(-net));
            error(j,i) = target(j) - actualOut(j,i);
            
            wS1 = wS1 + coeff * inS1(j) * error(j,i); % inputS1 weight update
            wS2 = wS2 + coeff * inS2(j) * error(j,i); % inputS2 weight update
            wS3 = wS3 + coeff * inS3(j) * error(j,i); % inputS3 weight update
            
            wP1 = wP1 + coeff * inP1(j) * error(j,i); % inputP1 weight update
            wP2 = wP2 + coeff * inP2(j) * error(j,i); % inputP2 weight update
            wP3 = wP3 + coeff * inP3(j) * error(j,i); % inputP3 weight update
            
            wBias = wBias + coeff * bias * error(j,i); % bias weight update
            
            weights = [weights; wS1, wS2, wS3, wP1, wP2, wP3, wBias];
        end
    end
    error = abs(error);

end