function [D_vec] = Analog_to_Digital (A_vec)
    % this function convert the analog inputs into digital inputs
    % Parameters: A_vec  = vector of analog inputs/signals (influx rates/concentrations)
    % Return: D_vec = vector of digital inputs/signals (bit streams)
    baseVal = max(A_vec);
    D_vec = round(A_vec/baseVal);
end