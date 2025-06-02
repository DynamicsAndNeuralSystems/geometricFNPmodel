function [a_array1, b_array1] = generate_connectome_edr(a_array,b_array, lambda, topology, isperiodic)
%GENERATE_CONNECTOME_EDR Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 5 || isempty(isperiodic)
        isperiodic = 1;
    end

    if ~all(size(a_array) == size(b_array))
        error('a_array and b_array must be same size');
    end
    num_fnps = size(a_array, 2);

    a_array1 = a_array; b_array1 = b_array;
    % Create connectome constrained by lambda parameter
    for j = 1:num_fnps
        a = a_array1(:, j); b = b_array1(:, j);
        test = 0;
        while (test == 0)
            distx = abs(a(1) - b(1));
            disty = abs(a(2) - b(2));
            if isperiodic == 1 %% Account for periodic conditions
                distx = min(distx, topology.L - distx);
                disty = min(disty, topology.L - disty);
            end
            dist = norm([distx; disty], 2);
            test = 1;
            if rand() > (exp(-dist*(100*lambda)))
                test = 0;
                a = topology.L * rand(2, 1);
                b = topology.L * rand(2, 1);                    
            end
        end
        a_array1(:, j) = a; b_array1(:, j) = b;
    end
end

