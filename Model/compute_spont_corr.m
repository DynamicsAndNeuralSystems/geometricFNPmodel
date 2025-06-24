function corr_array = compute_spont_corr(topology, dx, dt, hetparam_het, homparam, num_parcel_x, tol)

    % Compute the SCF with reference point p, with given sigma_n

    % Set tolerance if not given
    if nargin < 7 || isempty(tol)
        tol = 1e-7;
    end

    parcel_Nx = topology.Nx / num_parcel_x;
    if mod(topology.Nx, num_parcel_x) ~= 0
        error('Error: N_x must be divisible by number of parcellations per dimension');
    end
    
    % Create stim struct
    sum_xy = zeros(num_parcel_x, num_parcel_x, num_parcel_x, num_parcel_x);
    sum_x2 = zeros(num_parcel_x, num_parcel_x);
    sum_x = zeros(num_parcel_x, num_parcel_x);
    init = 0;
    cov_array = zeros(num_parcel_x, num_parcel_x, num_parcel_x, num_parcel_x);
    std_array = zeros(num_parcel_x, num_parcel_x);
    corr_array = zeros(num_parcel_x, num_parcel_x, num_parcel_x, num_parcel_x);
        
    % Iterate until convergence
    delta = Inf;
    tottime = 0;
    while(delta > tol || isnan(delta))
        tottime = tottime + topology.Nt;
        % Create noisy input
        stiminput = (1/(sqrt(dt) * dx)) * randn(topology.Nx, topology.Nx, topology.Nt);
        % Run noise driven activity
        ts = run_init_cstminput_periodic(topology, homparam, hetparam_het, stiminput, init);
        init = ts(:, :, end-1:end);
        % Parcellate activity
        ts1 = zeros(num_parcel_x, num_parcel_x, topology.Nt);
        for i = 1:num_parcel_x
            for j = 1:num_parcel_x
                ts1(i, j, :) = sum(ts(1+(i-1)*parcel_Nx:i*parcel_Nx, 1+(j-1)*parcel_Nx:j*parcel_Nx, :), [1 2]);
            end
        end
        % Update Correlation structure
        for i = 1:num_parcel_x
            for j = 1:num_parcel_x
                ts_ij = ts1(i, j, :);
                sum_x2(i, j) = sum_x2(i, j) + sum(ts_ij.*ts_ij);
                sum_x(i, j) = sum_x(i, j) + sum(ts_ij);
                % Calculate diagonal sum_xy
                sum_xy(i, j, 1:(i-1), :) = NaN;
                sum_xy(i, j, i, 1:j) = NaN;
                for j1 = (j+1):num_parcel_x
                    ts_ij1 = ts1(i, j1, :);
                    sum_xy(i, j, i, j1) = sum_xy(i, j, i, j1) + sum(ts_ij.*ts_ij1);
                end                
                for i1 = (i+1):num_parcel_x
                    for j1 = 1:num_parcel_x
                        ts_ij1 = ts1(i1, j1, :);
                        sum_xy(i, j, i1, j1) = sum_xy(i, j, i1, j1) + sum(ts_ij.*ts_ij1);
                    end
                end
            end
        end
        cov_array = (sum_xy / tottime) - ...
            (reshape(sum_x, [num_parcel_x, num_parcel_x, 1, 1]) / tottime) .* ...
            (reshape(sum_x, [1, 1, num_parcel_x, num_parcel_x]) / tottime);
        std_array = sqrt((sum_x2 / tottime) - (sum_x / tottime).^2);
        corr_array_old = corr_array;
        corr_array = cov_array ./ ...
            (reshape(std_array, [num_parcel_x, num_parcel_x, 1, 1]) .* ...
            reshape(std_array, [1, 1, num_parcel_x, num_parcel_x])); 
        
        % Compute change in SCF from update
        if all(corr_array(:) == 0) || all(corr_array_old(:) == 0)
            delta = NaN;
        else
            vec1 = corr_array_old(:);
            vec1 = vec1(~isnan(vec1));
            vec2 = corr_array(:);
            vec2 = vec2(~isnan(vec2));
            delta = pdist2(vec1', vec2', 'cosine');
            vec1 = 0; vec2 = 0;
        end
        
        disp(num2str(delta));
    end
end
