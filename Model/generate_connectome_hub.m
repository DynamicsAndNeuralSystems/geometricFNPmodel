function [a_array1, b_array1] = generate_connectome_hub(a_array,b_array, lambda, topology, hub_centres, hublength)
%GENERATE_CONNECTOME_EDR Summary of this function goes here
%   Detailed explanation goes here

    if ~all(size(a_array) == size(b_array))
        error('a_array and b_array must be same size');
    end
    num_hubs = size(hub_centres, 1);
    num_fnps = size(a_array, 2);

    a_array1 = a_array; b_array1 = b_array;

    % Create connectome

    for j = 1:num_fnps

        a = a_array1(:, j); b = b_array1(:, j);
        test = 0;
        reject_rv = rand();
        while (test == 0)
            test = 1;
            if reject_rv < lambda
                % Check if LRC is not incident with any hub region
                % Calculate L-inf distance of Ri from hub centre
                dist_a = zeros(1, num_hubs);
                for iter0 = 1:num_hubs
                    distx = abs(a(1) - hub_centres(iter0, 1));
                    disty = abs(a(2) - hub_centres(iter0, 2));
                    dist_a(iter0) = norm([distx; disty], Inf);
                end
                % Calculate L-inf distance of Rf from hub centre
                dist_b = zeros(1, num_hubs);
                for iter0 = 1:num_hubs
                    distx = abs(b(1) - hub_centres(iter0, 1));
                    disty = abs(b(2) - hub_centres(iter0, 2));
                    dist_b(iter0) = norm([distx; disty], Inf);
                end
                % Rejection 1: a and b are both inside a hub
                test1 = any(dist_a <= hublength/2) & any(dist_b <= hublength/2);
                % Rejection 2: a and b are both outside a hub
                test2 = all(dist_a > hublength/2) & all(dist_b > hublength/2);
                if test1 || test2
                    test = 0;
                    a = topology.L * rand(2, 1);
                    b = topology.L * rand(2, 1);   
                end
            end
        end
        a_array1(:, j) = a; b_array1(:, j) = b;
    end
end

