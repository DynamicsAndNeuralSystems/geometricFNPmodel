%% Create schematics of topological constraints

clear; clc;

loadparam;

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

f = figure;
f.Position = [100, 100, 1400, 175];
set(gca,'Color','white')
box on;
t = tiledlayout(1, 28, 'TileSpacing', 'none', 'Padding', 'compact');

num_lrcs = 10;

% Create Ensemble of Connectomes

rng(1, "twister");

a_array = zeros(2, num_lrcs);
b_array = zeros(2, num_lrcs);

for j = 1:num_lrcs
    a = topology.L * rand(2, 1);
    b = topology.L * rand(2, 1);
    a_array(:, j) = a;
    b_array(:, j) = b;
end

% Create schematics of default topology with EDR

lambda_array = [0, 0.5, 1];

dist_array = zeros(1, num_lrcs);

for i = 1:length(lambda_array)

    rng(1, "twister");

    lambda = lambda_array(i);

    a_array1 = a_array; b_array1 = b_array;

    dist_array = zeros(1, num_lrcs);

    % Create connectome

    for j = 1:num_lrcs

        a = a_array1(:, j); b = b_array1(:, j);
        test = 0;
        while (test == 0)
            distx = abs(a(1) - b(1));
            disty = abs(a(2) - b(2));
            dist = norm([distx; disty], 2);
            dist_array(j) = dist;
            test = 1;
            if rand() > (exp(-dist*(100*lambda)))
                test = 0;
                a = topology.L * rand(2, 1);
                b = topology.L * rand(2, 1);                    
            end
        end
        a_array1(:, j) = a; b_array1(:, j) = b;
        dist_array(j) = dist;
    end

    ax = nexttile(3*i - 2, [1 2]);
    hold;
    set(gca, 'Color', 'white');
    ax.Box = "on";
    ax.LineWidth = 1;

    xlim([0 topology.L])
    ylim([0 topology.L])
    xticks([]);
    yticks([]);
    for k = 1:num_lrcs
        quiver(a_array1(1, k), a_array1(2, k), ...
        b_array1(1, k) - a_array1(1, k), ...
        b_array1(2, k) - a_array1(2, k),...
        'Color', 'k', 'LineWidth', 0.5, ...
        'MaxHeadSize', 0.05 / norm(b_array1(:, k) - a_array1(:, k)), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end

    xlabel(append('$\lambda_e = ', num2str(lambda), '$'), 'Interpreter', 'latex', 'FontSize', 16);

    if i == 2
        t1 = title('Exponential Distance Rule', 'Interpreter', 'latex', 'FontSize', 16);
        set(t1, 'Units', 'normalized');
        t1.Position(2) = 1.1;
    end

    if i == 1
        t2 = title('\textbf{i.}', 'Interpreter', 'latex', 'FontSize', 24);
        ax.TitleHorizontalAlignment = 'left';
    end

    hold off;

end


for i = 1:2
    nexttile(3*i);
    hold;
    quiver(0.35, 0.5, 0.3, 0, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 1, 'Marker', '.', 'MarkerSize', 0.0001, 'AutoScale', 'off');
    xlim([0 1])
    ylim([0 1])
    xticks([]);
    yticks([]);
    axis off;
end


% Create schematics of default topology with hubs

% Set hub region position
num_hubs = 4;
hub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
hublength = topology.L / sqrt(34);

lambda_array = [0, 0.5, 1];

for i = 1:length(lambda_array)

    rng(1, "twister");

    lambda = lambda_array(i);

    a_array1 = a_array; b_array1 = b_array;

    % Create connectome

    for j = 1:num_lrcs

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

    ax = nexttile(10 + 3*i - 2, [1 2]);
    hold;
    set(gca, 'Color', 'white');
    ax.Box = "on";
    ax.LineWidth = 1;

    xlim([0 topology.L])
    ylim([0 topology.L])
    xticks([]);
    yticks([]);

    % Draw hub region
    for k = 1:num_hubs
        fill(hub_centres(k, 1) + hublength/2 * [-1 -1 1 1], hub_centres(k, 2) + hublength/2 * [-1 1 1 -1], 'k', 'FaceColor', '#D95319', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    end

    for k = 1:num_lrcs
        quiver(a_array1(1, k), a_array1(2, k), ...
        b_array1(1, k) - a_array1(1, k), ...
        b_array1(2, k) - a_array1(2, k),...
        'Color', 'k', 'LineWidth', 0.5, ...
        'MaxHeadSize', 0.05 / norm(b_array1(:, k) - a_array1(:, k)), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end

    xlabel(append('$\lambda_h = ', num2str(lambda), '$'), 'Interpreter', 'latex', 'FontSize', 16);

    if i == 2
        t1 = title('Hub Specificity', 'Interpreter', 'latex', 'FontSize', 16);
        set(t1, 'Units', 'normalized');
        t1.Position(2) = 1.1;
    end

    if i == 1
        title('\textbf{ii.}', 'Interpreter', 'latex', 'FontSize', 24, 'HorizontalAlignment', 'left');
        ax.TitleHorizontalAlignment = 'left';
    end

    hold off;

end



for i = 1:2
    nexttile(10 + 3*i);
    hold;
    quiver(0.35, 0.5, 0.3, 0, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 1, 'Marker', '.', 'MarkerSize', 0.0001, 'AutoScale', 'off');
    xlim([0 1])
    ylim([0 1])
    xticks([]);
    yticks([]);
    axis off;
end

% Create schematics of default topology with rich club

% Set number of rich club
% Set rich club region position
num_richclubs = 4;
richclub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
richclublength = topology.L / sqrt(34);

lambda_array = [0, 0.5, 1];

for i = 1:length(lambda_array)

    rng(1, "twister");

    lambda = lambda_array(i);

    a_array1 = a_array; b_array1 = b_array;

    % Create connectome

    for j = 1:num_lrcs

        a = a_array1(:, j); b = b_array1(:, j);
        test = 0;
        reject_rv = rand();
        while (test == 0)
            distx = abs(a(1) - b(1));
            disty = abs(a(2) - b(2));
            dist = norm([distx; disty], 2);
            test = 1;
            if reject_rv < lambda
                % Check if LRC is not incident with any hub region
                % Calculate L-inf distance of Ri from hub centre
                dist_a = zeros(1, num_richclubs);
                for iter0 = 1:num_richclubs
                    distx = abs(a(1) - richclub_centres(iter0, 1));
                    disty = abs(a(2) - richclub_centres(iter0, 2));
                    dist_a(iter0) = norm([distx; disty], Inf);
                end
                % Calculate L-inf distance of Rf from hub centre
                dist_b = zeros(1, num_richclubs);
                for iter0 = 1:num_richclubs
                    distx = abs(b(1) - richclub_centres(iter0, 1));
                    disty = abs(b(2) - richclub_centres(iter0, 2));
                    dist_b(iter0) = norm([distx; disty], Inf);
                end
                % Rejection 1: either a or b are outside a hub
                test1 = all(dist_a > richclublength/2) | all(dist_b > richclublength/2);
                % Rejection 2: a and b are both inside the same hub
                test2 = any(dist_a <= richclublength/2 & dist_b <= richclublength/2);
                if test1 || test2
                    test = 0;
                    a = topology.L * rand(2, 1);
                    b = topology.L * rand(2, 1);     
                end
            end
        end
        a_array1(:, j) = a; b_array1(:, j) = b;
    end

    ax = nexttile(20 + 3*i - 2, [1 2]);
    hold;
    set(gca, 'Color', 'white');
    ax.Box = "on";
    ax.LineWidth = 1;

    xlim([0 topology.L])
    ylim([0 topology.L])
    xticks([]);
    yticks([]);

    % Draw hub region
    for k = 1:num_richclubs
        fill(richclub_centres(k, 1) + richclublength/2 * [-1 -1 1 1], richclub_centres(k, 2) + richclublength/2 * [-1 1 1 -1], 'k', 'FaceColor', '#D95319', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    end

    for k = 1:num_lrcs
        quiver(a_array1(1, k), a_array1(2, k), ...
        b_array1(1, k) - a_array1(1, k), ...
        b_array1(2, k) - a_array1(2, k),...
        'Color', 'k', 'LineWidth', 0.5, ...
        'MaxHeadSize', 0.05 / norm(b_array1(:, k) - a_array1(:, k)), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end

    xlabel(append('$\lambda_r = ', num2str(lambda), '$'), 'Interpreter', 'latex', 'FontSize', 16);

    if i == 2
        t1 = title('Rich-club Specificity', 'Interpreter', 'latex', 'FontSize', 16);
        set(t1, 'Units', 'normalized');
        t1.Position(2) = 1.1;
    end

    if i == 1
        title('\textbf{iii.}', 'Interpreter', 'latex', 'FontSize', 24, 'HorizontalAlignment', 'left');
        ax.TitleHorizontalAlignment = 'left';
    end

    hold off;

end



for i = 1:2
    nexttile(20 + 3*i);
    hold;
    quiver(0.35, 0.5, 0.3, 0, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 1, 'Marker', '.', 'MarkerSize', 0.0001, 'AutoScale', 'off');
    xlim([0 1])
    ylim([0 1])
    xticks([]);
    yticks([]);
    axis off;
end

% save as schematictopologyconstraints.svg
print(gcf, 'schematic_nonrandom.svg', '-dsvg');
exportgraphics(gcf, 'schematic_nonrandom.tiff', 'Resolution', 300);

%% BOLD dissimilarity of random models

clear; clc;

loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

tol = 1e-5;

num_lrcs_array = [10, 20, 50, 100];
max_num_lrcs = max(num_lrcs_array);
num_samples = 100;

% Record BOLD divergence
md_array = zeros(length(num_lrcs_array), num_samples);

% Sample Stimulus Position
rng(0, "twister");
stim_positions = topology.L * rand(2, num_samples);

parfor k = 1:num_samples

    stim1 = stim;
    hetparam_het1 = hetparam_het;

    % Set Stimulus Position
    stim1 = stim;
    stim1.stimR = stim_positions(:, k);

    % Simulate homogeneous model
    ts_hom_bold = run_bold(topology, homparam, hetparam_hom, stim1, tol);

    md_subarray = zeros(length(num_lrcs_array), 1);

    % Create Ensemble of Connectomes

    rng(k, "twister");

    a_array = zeros(2, max_num_lrcs);
    b_array = zeros(2, max_num_lrcs);

    for j = 1:max_num_lrcs
        a = topology.L * rand(2, 1);
        b = topology.L * rand(2, 1);
        a_array(:, j) = a;
        b_array(:, j) = b;
    end

    for i = 1:length(num_lrcs_array)

        num_lrcs = num_lrcs_array(i);

        hetparam_het1.m = num_lrcs;
        hetparam_het1.c = (homparam.r)^2 * ones(1, num_lrcs);
        hetparam_het1.tau = zeros(1, num_lrcs);
        hetparam_het1.a = a_array(:, 1:num_lrcs);
        hetparam_het1.b = b_array(:, 1:num_lrcs);

        ts_het_bold = run_bold(topology, homparam, hetparam_het1, stim1, tol);

        dissim_bold = pdist2(...
            ts_hom_bold(:)', ts_het_bold(:)', 'cosine');

        md_subarray(i) = dissim_bold;
        disp([i, k, dissim_bold]);

    end

    md_array(:, k) = md_subarray;

end

save('dissimboldensemble_multilrcs.mat', "md_array", "num_samples", "num_lrcs_array", "topology");


%% BOLD dissimilarity versus EDR decay
% Fix number of shortcuts to 50

clear; clc;
 
loadparam;
 
dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

tol = 1e-5;

num_lrcs_array = [10, 20, 50, 100];
max_num_lrcs = max(num_lrcs_array);
num_samples = 100;

lambda_array = 0.2:0.2:1;
 
md_array = zeros(length(lambda_array), length(num_lrcs_array), num_samples);

% Sample Stimulus Position
rng(0, "twister");
stim_positions = topology.L * rand(2, num_samples);
 
parfor k = 1:num_samples

    stim1 = stim;
    hetparam_het1 = hetparam_het;

    % Set Stimulus Position
    stim1 = stim;
    stim1.stimR = stim_positions(:, k);

    % Simulate homogeneous model
    ts_hom_bold = run_bold(topology, homparam, hetparam_hom, stim1, tol);

    md_subarray = zeros(length(lambda_array), length(num_lrcs_array));

    % Create Ensemble of Connectomes
    
    rng(k, "twister");
    
    a_array = zeros(2, max_num_lrcs);
    b_array = zeros(2, max_num_lrcs);
    
    for j = 1:max_num_lrcs
        a = topology.L * rand(2, 1);
        b = topology.L * rand(2, 1);
        a_array(:, j) = a;
        b_array(:, j) = b;
    end

    for i = 1:length(lambda_array)

        lambda = lambda_array(i);

        a_array1 = a_array;
        b_array1 = b_array;
    
        % Create topological constrained connectome for given lambda

        for j = 1:max_num_lrcs
            
            a = a_array1(:, j); b = b_array1(:, j);
            test = 0;
            while(test == 0)
                test = 1;
                distx = abs(a(1) - b(1));
                distx = min(distx, topology.L - distx);
                disty = abs(a(2) - b(2));
                disty = min(disty, topology.L - disty);
                dist = norm([distx; disty], 2);
                % Rejection: exp(-lambda*l) is less than unif[0, 1]
                test1 = (rand() > (exp(-dist*(100*lambda))));
                if test1
                    test = 0;
                    a = topology.L * rand(2, 1);
                    b = topology.L * rand(2, 1);      
                end
            end
            a_array1(:, j) = a;
            b_array1(:, j) = b;
        end 
        for j = 1:length(num_lrcs_array)
    
            num_lrcs = num_lrcs_array(j);
    
            hetparam_het1.m = num_lrcs;
            hetparam_het1.c = (homparam.r)^2 * ones(1, num_lrcs);
            hetparam_het1.tau = zeros(1, num_lrcs);
            hetparam_het1.a = a_array1(:, 1:num_lrcs);
            hetparam_het1.b = b_array1(:, 1:num_lrcs);
    
            ts_het_bold = run_bold(topology, homparam, hetparam_het1, stim1, tol);
    
            dissim_bold = pdist2(...
                ts_hom_bold(:)', ts_het_bold(:)', 'cosine');
        
            md_subarray(i, j) = dissim_bold;
            disp(num2str([i, j, k, dissim_bold]));
    
        end
    
        md_array(:, :, k) = md_subarray;
    
    end

end
%
save('dissimboldensemble_edr.mat', "md_array", "num_samples", "lambda_array", "topology");
% 


%% BOLD Divergence versus hub specificity

clear; clc;
 
loadparam;
 
dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

tol = 1e-5;

num_lrcs_array = [10, 20, 50, 100];
max_num_lrcs = max(num_lrcs_array);
num_samples = 100;

% Set hub region position
num_hubs = 4;
hub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
hublength = topology.L / sqrt(34);

lambda_array = 0.2:0.2:1;

md_array = zeros(length(lambda_array), length(num_lrcs_array), num_samples);

% Sample Stimulus Position
rng(0, "twister");
stim_positions = topology.L * rand(2, num_samples);

parfor k = 1:num_samples

    stim1 = stim;
    hetparam_het1 = hetparam_het;

    % Set Stimulus Position
    stim1 = stim;
    stim1.stimR = stim_positions(:, k);

    % Simulate homogeneous model
    ts_hom_bold = run_bold(topology, homparam, hetparam_hom, stim1, tol);

    md_subarray = zeros(length(lambda_array), length(num_lrcs_array));

    % Create Ensemble of Connectomes
    
    rng(k, "twister");
    
    a_array = zeros(2, max_num_lrcs);
    b_array = zeros(2, max_num_lrcs);
    
    for j = 1:max_num_lrcs
        a = topology.L * rand(2, 1);
        b = topology.L * rand(2, 1);
        a_array(:, j) = a;
        b_array(:, j) = b;
    end
    
    for i = 1:length(lambda_array)       
        
        lambda = lambda_array(i);
        
        a_array1 = a_array;
        b_array1 = b_array;

        % Create topological constrained connectome for given lambda
    
        for j = 1:max_num_lrcs
    
            a = a_array1(:, j); b = b_array1(:, j);
            test = 0;
            reject_rv = rand();
            while (test == 0)
                test = 1;
                if reject_rv < lambda
                    % Check if LRC is not incident with any hub region
                    % Calculate L-inf distance of a from all hub centres
                    dist_a = zeros(1, num_hubs);
                    for iter0 = 1:num_hubs
                        distx = abs(a(1) - hub_centres(iter0, 1));
                        disty = abs(a(2) - hub_centres(iter0, 2));
                        dist_a(iter0) = norm([distx; disty], Inf);
                    end
                    % Calculate L-inf distance of b from all hub centres
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

        for j = 1:length(num_lrcs_array)
    
            num_lrcs = num_lrcs_array(j);
    
            hetparam_het1.m = num_lrcs;
            hetparam_het1.c = (homparam.r)^2 * ones(1, num_lrcs);
            hetparam_het1.tau = zeros(1, num_lrcs);
            hetparam_het1.a = a_array1(:, 1:num_lrcs);
            hetparam_het1.b = b_array1(:, 1:num_lrcs);
    
            ts_het_bold = run_bold(topology, homparam, hetparam_het1, stim1, tol);
    
            dissim_bold = pdist2(...
                ts_hom_bold(:)', ts_het_bold(:)', 'cosine');
        
            md_subarray(i, j) = dissim_bold;
            disp(num2str([i, j, k, dissim_bold]));
    
        end
    
        md_array(:, :, k) = md_subarray;
    
    end

end
%
save('dissimboldensemble_hub.mat', "md_array", "num_samples", "lambda_array", "num_hubs", "hub_centres", "hublength", "topology");
% 

%% BOLD Divergence versus rich-club specificity

clear; clc;
 
loadparam;
 
dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

tol = 1e-5;

num_lrcs_array = [10, 20, 50, 100];
max_num_lrcs = max(num_lrcs_array);
num_samples = 100;

% Set number of hubs: superior parietal, preceuneus, superior frontal 
% (van den heuvel jneurosci)
num_richclubs = 4;
richclub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
richclublength = topology.L / sqrt(34);

lambda_array = 0.2:0.2:1;

md_array = zeros(length(lambda_array), length(num_lrcs_array), num_samples);

% Sample Stimulus Position
rng(0, "twister");
stim_positions = topology.L * rand(2, num_samples);

parfor k = 1:num_samples

    stim1 = stim;
    hetparam_het1 = hetparam_het;

    % Set Stimulus Position
    stim1 = stim;
    stim1.stimR = stim_positions(:, k);
    
    % Simulate homogeneous model
    ts_hom_bold = run_bold(topology, homparam, hetparam_hom, stim1, tol);
    
    md_subarray = zeros(length(lambda_array), length(num_lrcs_array));
    
    % Create Ensemble of Connectomes
    
    rng(k, "twister");
    
    a_array = zeros(2, max_num_lrcs);
    b_array = zeros(2, max_num_lrcs);
    
    for j = 1:max_num_lrcs
        a = topology.L * rand(2, 1);
        b = topology.L * rand(2, 1);
        a_array(:, j) = a;
        b_array(:, j) = b;
    end
    
    for i = 1:length(lambda_array)

        lambda = lambda_array(i);

        a_array1 = a_array;
        b_array1 = b_array;

        % Create topological constrained connectome for given lambda

        for j = 1:max_num_lrcs
    
            a = a_array1(:, j); b = b_array1(:, j);
            test = 0;
            reject_rv = rand();
            while (test == 0)
                test = 1;
                if reject_rv < lambda
                    % Check if LRC is not incident with any hub region
                    % Calculate L-inf distance of Ri from hub centre
                    dist_a = zeros(1, num_richclubs);
                    for iter0 = 1:num_richclubs
                        distx = abs(a(1) - richclub_centres(iter0, 1));
                        disty = abs(a(2) - richclub_centres(iter0, 2));
                        dist_a(iter0) = norm([distx; disty], Inf);
                    end
                    % Calculate L-inf distance of Rf from hub centre
                    dist_b = zeros(1, num_richclubs);
                    for iter0 = 1:num_richclubs
                        distx = abs(b(1) - richclub_centres(iter0, 1));
                        disty = abs(b(2) - richclub_centres(iter0, 2));
                        dist_b(iter0) = norm([distx; disty], Inf);
                    end
                    % Rejection 1: either a or b are outside a hub
                    test1 = all(dist_a > richclublength/2) | all(dist_b > richclublength/2);
                    % Rejection 2: a and b are both inside the same hub
                    test2 = any(dist_a <= richclublength/2 & dist_b <= richclublength/2);
                    if test1 || test2
                        test = 0;
                        a = topology.L * rand(2, 1);
                        b = topology.L * rand(2, 1);     
                    end
                end
            end
            a_array1(:, j) = a; b_array1(:, j) = b;
        end
        
        for j = 1:length(num_lrcs_array)
    
            num_lrcs = num_lrcs_array(j);
    
            hetparam_het1.m = num_lrcs;
            hetparam_het1.c = (homparam.r)^2 * ones(1, num_lrcs);
            hetparam_het1.tau = zeros(1, num_lrcs);
            hetparam_het1.a = a_array1(:, 1:num_lrcs);
            hetparam_het1.b = b_array1(:, 1:num_lrcs);
    
            ts_het_bold = run_bold(topology, homparam, hetparam_het1, stim1, tol);
    
            dissim_bold = pdist2(...
                ts_hom_bold(:)', ts_het_bold(:)', 'cosine');
        
            md_subarray(i, j) = dissim_bold;
            disp(num2str([i, j, k, dissim_bold]));
    
        end
    
        md_array(:, :, k) = md_subarray;
    
    end

end
%
save('dissimboldensemble_core.mat', "md_array", "num_samples", "lambda_array", "num_richclubs", "richclub_centres", "richclublength", "topology");
% 

%% Plot mean BOLD dissimilarity versus three topological control parameters

clear; clc;
loadparam;

f = figure;
f.Position = [100, 100, 1400, 500];
hold;
title(['Cosine Distance vs. Probability of Hub Incidence'])
set(gcf, 'Color', 'white')
color_array = get(gca,'colororder');


% Load random connectome statistics
load('dissimboldensemble_multilrcs.mat', 'md_array', 'num_lrcs_array')
% Keep only 10, 20, 50 and 100 LRCs
md_array0 = md_array;
clear md_array;
md_array0 = reshape(md_array0, [1 size(md_array0)]);

lambda = 188;
reciprobability = 0.31;
hubprobability = 1 - sqrt(12 / 68); % 12 hub regions identified from 68 across both hemispheres in VDH 2011

t = tiledlayout(1, 28, 'TileSpacing', 'None', 'Padding', 'Compact');

ax = nexttile(10*1 - 9, [1 8]); hold;

load('dissimboldensemble_edr.mat');

% Append random connectome statistics
md_array = cat(1, md_array0, md_array);
lambda_array = [0 lambda_array];

mean_dissimilarity = squeeze(mean(md_array, 3));
std_dissimilarity = squeeze(std(md_array, 1, 3));

set(ax,'ColorOrderIndex',1);
for i = 1:length(num_lrcs_array)
    h = errorbar(lambda_array, mean_dissimilarity(:, i), std_dissimilarity(:, i), 'o', 'MarkerSize', 0.000001, 'LineWidth', 1);
    alpha = 0.2;   
    % Set transparency (undocumented)
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
end
for i = 1:length(num_lrcs_array)
    plot(lambda_array, mean_dissimilarity(:, i), 'Marker', 'o', 'Color', color_array(i, :), 'MarkerFaceColor', color_array(i, :));
end

xlabel('$\lambda_e$', 'Interpreter', 'latex')
ylabel('$C_z$', 'Interpreter', 'latex', 'Rotation', 0)
ax.TickLabelInterpreter = 'latex';
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.XAxis.LabelFontSizeMultiplier  = 1.5;
ax.YAxis.LabelFontSizeMultiplier  = 2;

xlim([-0.05*max(xlim) 1.05*max(xlim)])
yticks([0:0.02:0.1])
xticks(lambda_array)

title('\textbf{i.}', 'Interpreter', 'latex', 'FontSize', 24);
ax.TitleHorizontalAlignment = 'left';

hold off;

ax = nexttile(10*2 - 9, [1 8]); hold;

load('dissimboldensemble_hub.mat');


% Append random connectome statistics
md_array = cat(1, md_array0, md_array);
lambda_array = [0 lambda_array];

mean_dissimilarity = squeeze(mean(md_array, 3));
std_dissimilarity = squeeze(std(md_array, 1, 3));

set(ax,'ColorOrderIndex',1);
for i = 1:length(num_lrcs_array)
    h = errorbar(lambda_array, mean_dissimilarity(:, i), std_dissimilarity(:, i), 'o', 'MarkerSize', 0.000001, 'LineWidth', 1);
    alpha = 0.2;   
    % Set transparency (undocumented)
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
end
for i = 1:length(num_lrcs_array)
    plot(lambda_array, mean_dissimilarity(:, i), 'Marker', 'o', 'Color', color_array(i, :), 'MarkerFaceColor', color_array(i, :));
end

xlabel('$\lambda_h$', 'Interpreter', 'latex');
ax.TickLabelInterpreter = 'latex';
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.XAxis.LabelFontSizeMultiplier  = 1.5;

xlim([-0.05*max(xlim) 1.05*max(xlim)])
yticks([0:0.02:0.1])
xticks(lambda_array);

title('\textbf{ii.}', 'Interpreter', 'latex', 'FontSize', 24);
ax.TitleHorizontalAlignment = 'left';

hold off;

ax = nexttile(10*3 - 9, [1 8]); hold;

load('dissimboldensemble_core.mat');

% Append random connectome statistics
md_array = cat(1, md_array0, md_array);
lambda_array = [0 lambda_array];

mean_dissimilarity = squeeze(mean(md_array, 3));
std_dissimilarity = squeeze(std(md_array, 1, 3));

set(ax,'ColorOrderIndex',1);
for i = 1:length(num_lrcs_array)
    h = errorbar(lambda_array, mean_dissimilarity(:, i), std_dissimilarity(:, i), 'o', 'MarkerSize', 0.000001, 'LineWidth', 1);
    alpha = 0.2;   
    % Set transparency (undocumented)
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
end
for i = 1:length(num_lrcs_array)
    plot(lambda_array, mean_dissimilarity(:, i), 'Marker', 'o', 'Color', color_array(i, :), 'MarkerFaceColor', color_array(i, :));
end

xlabel('$\lambda_r$', 'Interpreter', 'latex');
ax.TickLabelInterpreter = 'latex';
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.XAxis.LabelFontSizeMultiplier  = 1.5;

xlim([-0.05*max(xlim) 1.05*max(xlim)])
yticks([0:0.02:0.1])
xticks(lambda_array);

title('\textbf{iii.}', 'Interpreter', 'latex', 'FontSize', 24);
ax.TitleHorizontalAlignment = 'left';

hold off;

Legend = cell(1, 2*length(num_lrcs_array));
for i = 1:length(num_lrcs_array)
    Legend{i} = "";
    Legend{length(num_lrcs_array) + i} = append('$N = ', num2str(num_lrcs_array(i)), '$');
end

l = legend(Legend, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 16);
l.Layout.Tile = 'south';
l.Orientation = 'horizontal';

% save as dissimboldensemble_topologyconstraints.svg
print(gcf, 'dissimboldensemble_nonrandom.svg', '-dsvg');
exportgraphics(gcf, 'dissimboldensemble_nonrandom.tiff', 'Resolution', 300);

%% BOLD dissimilarity of sampled random, hub and rich-club connectomes

clear; clc;

loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

num_lrcs = 100;
num_models = 3;

% Set hubs: superior parietal, preceuneus, superior frontal 
% (van den heuvel jneurosci)
num_hubs = 4;
hub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
hublength = topology.L / sqrt(34);
% Set rich clubs: superior parietal, preceuneus, superior frontal 
% (van den heuvel jneurosci)
num_richclubs = 4;
richclub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
richclublength = topology.L / sqrt(34);

% Create Ensemble of Connectomes

rng(1, "twister");

a_array = zeros(2, num_lrcs);
b_array = zeros(2, num_lrcs);

for j = 1:num_lrcs
    a = topology.L * rand(2, 1);
    b = topology.L * rand(2, 1);
    a_array(:, j) = a;
    b_array(:, j) = b;
end

hetparam_het_array = cell(1, num_models);

% Create random, hub and core connectome

lambda = 1;

for iter = 1:num_models
    a_array1 = a_array; b_array1 = b_array;
    if iter == 2
        rng(1, "twister")
        for j = 1:num_lrcs
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
    elseif iter == 3
        rng(1, "twister")
        for j = 1:num_lrcs
            a = a_array1(:, j); b = b_array1(:, j);
            test = 0;
            reject_rv = rand();
            while (test == 0)
                distx = abs(a(1) - b(1));
                disty = abs(a(2) - b(2));
                dist = norm([distx; disty], 2);
                test = 1;
                if reject_rv < lambda
                    % Check if LRC is not incident with any hub region
                    % Calculate L-inf distance of Ri from hub centre
                    dist_a = zeros(1, num_richclubs);
                    for iter0 = 1:num_richclubs
                        distx = abs(a(1) - richclub_centres(iter0, 1));
                        disty = abs(a(2) - richclub_centres(iter0, 2));
                        dist_a(iter0) = norm([distx; disty], Inf);
                    end
                    % Calculate L-inf distance of Rf from hub centre
                    dist_b = zeros(1, num_richclubs);
                    for iter0 = 1:num_richclubs
                        distx = abs(b(1) - richclub_centres(iter0, 1));
                        disty = abs(b(2) - richclub_centres(iter0, 2));
                        dist_b(iter0) = norm([distx; disty], Inf);
                    end
                    % Rejection 1: either a or b are outside a hub
                    test1 = all(dist_a > richclublength/2) | all(dist_b > richclublength/2);
                    % Rejection 2: a and b are both inside the same hub
                    test2 = any(dist_a <= richclublength/2 & dist_b <= richclublength/2);
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
    hetparam_het1 = hetparam_het;
    hetparam_het1.m = num_lrcs;
    hetparam_het1.c = (homparam.r)^2 * ones(1, num_lrcs);
    hetparam_het1.tau = zeros(1, num_lrcs);
    hetparam_het1.a = a_array1(:, 1:num_lrcs);
    hetparam_het1.b = b_array1(:, 1:num_lrcs);
    hetparam_het_array{iter} = hetparam_het1;
end


% Set Stimulus Positions
num_samples_x = 40;
num_samples = num_samples_x.^2;
[stim_x, stim_y] = meshgrid(1:num_samples_x, 1:num_samples_x);
stim_positions = (topology.L / num_samples_x)*([stim_x(:), stim_y(:)]');

tol = 1e-5;

md_array = zeros(2 + num_models, num_samples);
parfor k = 1:num_samples

    stim1 = stim;
    stim1.stimR = stim_positions(:, k);

    % Simulate homogeneous model
    ts_hom_bold = run_bold(topology, homparam, hetparam_hom, stim1, tol);

    md_subarray = zeros(2 + num_models, 1);
    md_subarray(1:2) = stim1.stimR;
    
    for iter = 1:num_models
        hetparam_het = hetparam_het_array{iter};
        ts_het_bold = run_bold(topology, homparam, hetparam_het, stim1, tol);
        md_subarray(2 + iter) = pdist2(ts_hom_bold(:)', ts_het_bold(:)', 'cosine');
        disp([k, iter]);
    end
    

    md_array(:, k) = md_subarray;
    
end
%
save('bolddissimrandomvshubclub.mat', "md_array", "num_models", "hetparam_het_array", "num_samples", "num_samples_x", "num_hubs", "hub_centres", "hublength", "num_richclubs", "richclub_centres", "richclublength", "topology");

%% Plot spatial distribution of BOLD dissimilarity of sampled random, hub and rich-club connectomes

clear; clc;
loadparam;

Colormap = [linspace(1, 0, 256)' linspace(1, 0, 256)', linspace(1, 1, 256)'];

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

load('bolddissimrandomvshubclub.mat');

[xq,yq] = ndgrid(dx*(1:topology.Nx));
[stim_x, stim_y] = ndgrid(linspace(0, topology.L, num_samples_x + 1), linspace(0, topology.L, num_samples_x + 1));
dissim0 = reshape(md_array(3:end, :, :), [num_models, num_samples_x, num_samples_x]);

f = figure;
f.Position = [100 100 760 450];
num_subcols = 8;
num_subrows = 3;
tot_subcols = 2 + num_models*num_subcols;
tot_subrows = 1 + 2*num_subrows;
t = tiledlayout(tot_subrows, tot_subcols, 'TileSpacing', 'tight', 'Padding', 'compact');

modeltitles = ["Random", "Hub ($\lambda_h = 1$)", "Rich-club ($\lambda_r = 1$)"];
connectivitytitles = "\textbf{" + ["i", "ii", "iii"] + ".}";
dissimtitles = "\textbf{" + ["iv", "v", "vi"] + ".}";

% title(t, '$C_z$ of sample models by stimulus position', 'Interpreter', 'latex', 'FontSize', 24);


for iter = 1:num_models
    ax = nexttile(3 + (iter - 1)*(num_subcols), [1, num_subcols]);
    hold on;
    axis off;
    text(0.5, 0, modeltitles(iter), 'Interpreter', 'latex', 'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Units', 'normalized');
    hold off;
end

ax = nexttile(tot_subcols + 1, [num_subrows, 2]);
axis off;
text(0, 0.5, 'Connectivity', 'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90, 'Units', 'normalized')

for iter = 1:num_models
    ax = nexttile(tot_subcols + 3 + (iter - 1)*(num_subcols), [num_subrows, num_subcols]);
    hold on;
    box on;
    ax.LineWidth = 1;
    ax.Color = 'k';
    set(gca, 'Color', 'white');
    hetparam_het = hetparam_het_array{iter};
    if iter == 2
        % Draw hub region
        for k = 1:num_hubs
            fill(hub_centres(k, 1) + hublength/2 * [-1 -1 1 1], hub_centres(k, 2) + hublength/2 * [-1 1 1 -1], 'k', 'FaceColor', '#D95319', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        end
    elseif iter == 3
        % Draw hub region
        for k = 1:num_richclubs
            fill(richclub_centres(k, 1) + richclublength/2 * [-1 -1 1 1], richclub_centres(k, 2) + richclublength/2 * [-1 1 1 -1], 'k', 'FaceColor', '#D95319', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        end
    end
    for m = 1:hetparam_het.m
        q = quiver(hetparam_het.a(1, m), hetparam_het.a(2, m), ...
        hetparam_het.b(1, m) - hetparam_het.a(1, m), ...
        hetparam_het.b(2, m) - hetparam_het.a(2, m),...
        'Color', 'k', 'LineWidth', 0.1, ...
        'MaxHeadSize', 0.05 / norm(hetparam_het.a(:, m) - hetparam_het.b(:, m)), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end
    ylabel([connectivitytitles(iter), '\\', '\\', '\\'], 'Interpreter', 'latex', 'FontSize', 24, 'Rotation', 0);
    xlim([0, topology.L + dx]);
    ylim([0, topology.L + dx]);
    xticks([]);
    yticks([]);
    hold off;
end

ax = nexttile((num_subrows + 1)*tot_subcols + 1, [num_subrows, 2]);
axis off;
text(0, 0.5, ['$C_z$ by stimulus', newline, 'position'], 'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90, 'Units', 'normalized')

for iter = 1:num_models
    ax = nexttile((num_subrows + 1)*tot_subcols + 3 + (iter - 1)*(num_subcols), [num_subrows, num_subcols]);
    hold on;
    box on;
    % ax.LineWidth = 1;
    % ax.Color = 'k';
    hetparam_het = hetparam_het_array{iter};
    dissim = squeeze(dissim0(iter, :, :));
    % Create boundary conditions, introduce data points with x = 0 OR y = 0
    dissim = [dissim(end, end) dissim(end, :);
        dissim(:, end) dissim];
    % F = griddedInterpolant(stim_x,stim_y,dissim, 'spline');
    % vq = F(xq,yq);
    imagesc((topology.L / num_samples_x)*(0:num_samples_x),(topology.L / num_samples_x)*(0:num_samples_x),dissim); 
    for m = 1:hetparam_het.m
        scatter(hetparam_het.a(1, m), hetparam_het.a(2, m), 5, 'k', 'filled');
    end
    clim([0  max(md_array(3:end, :, :), [], 'all')]); 
    shading flat; view(0, 90); colormap(Colormap)
    set(gca, 'YDir', 'normal');
    ylabel([dissimtitles(iter), '\\', '\\', '\\'], 'Interpreter', 'latex', 'FontSize', 24, 'Rotation', 0);
    xlim([0, topology.L + dx]);
    ylim([0, topology.L + dx]);
    xticks([]);
    yticks([]);
    hold off;
end

cb = colorbar; 
cb.FontSize = 15;
cb.TickLabelInterpreter = 'latex';
cb.Ticks = [0, 0.1];
cb.TickLabels = {'0', '0.1'};

for iter = 1:num_models
    ax = nexttile((num_subrows + 1)*tot_subcols + 3 + (iter - 1)*(num_subcols), [num_subrows, num_subcols]);
    hold on;
    pos = get(ax, 'Position');
    annotation('rectangle', [pos(1), pos(2), pos(1) + pos(3) - pos(1), pos(2) + pos(4) - pos(2)], ...
        'LineWidth', 1, 'Color', 'k');
end

% save as dissimboldsamplerandomhubcore_spatialplot.svg
print(gcf, 'bolddissimrandomvshubclub_spatialplot.svg', '-dsvg');
exportgraphics(gcf, 'bolddissimrandomvshubclub_spatialplot.tiff', 'Resolution', 300);

%% Distribution of BOLD dissimilarity of sampled random, hub and rich-club connectomes

clear; clc;
loadparam;

Colormap = [linspace(1, 0, 256)' linspace(1, 0, 256)', linspace(1, 1, 256)'];

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

load('bolddissimrandomvshubclub.mat');

md_array = md_array(3:end, :);

f = figure;
f.Position = [100 100 600 430];
t = tiledlayout(1, 6, 'TileSpacing', 'tight', 'Padding', 'compact');

ax = nexttile(1, [1 1]);
axis off;
text(0, 1, '\textbf{vii.}', 'Interpreter', 'latex', 'FontSize', 24, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Units', 'normalized')

ax = nexttile(2, [1 5]);
hold;

colors = get(gca, 'ColorOrder');  % Retrieve the default color order

for i = 1:num_models

    xline(mean(md_array(i, :)), 'Color', colors(i, :), 'LineWidth', 2);

end

for i = 1:num_models

    data = md_array(i, :);

    % Calculate the probability density (using the 'ksdensity' function)
    [f, xi] = ksdensity(data, 'Support', 'positive');

    % Fill the area under the density curve with the same color as the line
    fill(xi, f, colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

end

xticks([0 : 0.05 : 0.5]); xlim([0, 0.15]);

ylabel("Density", 'Interpreter', 'latex', 'Rotation', 90);
xlabel('$C_z$', 'Interpreter', 'latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.XAxis.LabelFontSizeMultiplier  = 2;
ax.YAxis.LabelFontSizeMultiplier  = 1.5;


Legend = {"Random", "Hub ($\lambda_h = 1$)", "Rich-club ($\lambda_r = 1$)"};
l = legend(Legend, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 12);

% save as dissimboldsamplerandomhubclub_distribution.svg
print(gcf, 'dissimboldsamplerandomhubclub_distribution.svg', '-dsvg');
exportgraphics(gcf, 'dissimboldsamplerandomhubclub_distribution.tiff', 'Resolution', 300);