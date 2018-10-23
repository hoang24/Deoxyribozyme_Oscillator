function [HamVec, HamDist] = HammingDistance (A_inputs)
    % decompress table of inputs into analog inputs
    A_Sm1 = A_inputs(:, 1);
    A_Sm2 = A_inputs(:, 2); 
    
    % convert analog signals to digital inputs
    D_Sm1 = Analog_to_Digital(A_Sm1);
    D_Sm2 = Analog_to_Digital(A_Sm2);
    
    HamVec = (D_Sm1 == D_Sm2); % Hamming vector
    HamDist = sum(HamVec); % Hamming Distance
end