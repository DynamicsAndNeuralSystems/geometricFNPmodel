function corr_array = compute_scf(sigma_n, position_centre_grid, topology, dx, dt, hetparam_het, homparam, tol)

    % Compute the SCF with reference point p, with given sigma_n

    % Set tolerance if not given
    if nargin < 8 || isempty(tol)
        tol = 1e-7;
    end
    
    % Create stim struct
    sum_xy = zeros(topology.Nx, topology.Nx);
    sum_x2 = zeros(topology.Nx, topology.Nx);
    sum_x = zeros(topology.Nx, topology.Nx);
    sum_y2 = 0;
    sum_y = 0;
    init = 0;
    cov_array = zeros(topology.Nx, topology.Nx);
    std_array = zeros(topology.Nx, topology.Nx);
    stdy = 0;
    corr_array = zeros(topology.Nx, topology.Nx);
    
    % Reference point coordinates
    i0 = position_centre_grid(1); 
    j0 = position_centre_grid(2);
    
    % Pre-compute stimulation input space
    stimulationinputspace = zeros(topology.Nx, topology.Nx);
    for i = 1:topology.Nx
        for j = 1:topology.Nx
            distx = abs(i - i0);
            if distx > topology.Nx / 2
                distx = topology.Nx - distx;
            end
            disty = abs(j - j0);
            if disty > topology.Nx / 2
                disty = topology.Nx - disty;
            end
            stimulationinputspace(i, j) = exp(-0.5*(distx^2 + disty^2) * (dx^2) / (sigma_n^2));
        end
    end
    
    % Iterate until convergence
    delta = Inf;
    tottime = 0;
    while(delta > tol || isnan(delta))
        tottime = tottime + topology.Nt;
        % Create noisy input
        stiminput = (1/(sqrt(dt) * dx)) * randn(topology.Nx, topology.Nx, topology.Nt);
        if sigma_n < Inf
            stiminput = stiminput .* reshape(stimulationinputspace, [topology.Nx, topology.Nx, 1]);
        end
        % Run noise driven activity
        ts = run_init_cstminput_periodic(topology, homparam, hetparam_het, stiminput, init);
        init = ts(:, :, end-1:end);
        ts_centre = ts(i0, j0, :);

        % Update SCF
        sum_y = sum_y + sum(ts_centre);
        sum_y2 = sum_y2 + sum(ts_centre.^2);
        for i = 1:topology.Nx
            for j = 1:topology.Nx
                ts_point = ts(i, j, :);
                sum_xy(i, j) = sum_xy(i, j) + sum(ts_centre.*ts_point);
                sum_x2(i, j) = sum_x2(i, j) + sum(ts_point.*ts_point);
                sum_x(i, j) = sum_x(i, j) + sum(ts_point);
            end
        end
        cov_array = (sum_xy / tottime) - (sum_x / tottime) * (sum_y / tottime);
        std_array = sqrt((sum_x2 / tottime) - (sum_x / tottime).^2);
        stdy = sqrt((sum_y2 / tottime) - (sum_y / tottime).^2);
        corr_array_old = corr_array;
        corr_array = cov_array ./ (std_array * stdy);
        
        % Compute change in SCF from update
        if all(corr_array(:) == 0) || all(corr_array_old(:) == 0)
            delta = NaN;
        else
            delta = pdist2(corr_array_old(:)', corr_array(:)', 'cosine');
        end
        
        disp(num2str([sigma_n, delta]));
    end
end
