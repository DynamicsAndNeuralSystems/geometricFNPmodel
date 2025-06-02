%% Create schematics of connectomes with nonrandom properties

clear; clc;
% Load model parameters
loadparam;

% Set number of FNPs in schematics
N = 10;

% Create Random ensemble of Connectomes with specified number of FNPs
rng(1, "twister");
a_array = zeros(2, N);
b_array = zeros(2, N);
for j = 1:N
    a = topology.L * rand(2, 1);
    b = topology.L * rand(2, 1);
    a_array(:, j) = a;
    b_array(:, j) = b;
end

% Set lambda values
lambda_array = [0, 0.5, 1];

% Set hub region position
num_hubs = 4;
hub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
hublength = topology.L / sqrt(34);

% Create ensemble of connectomes with nonrandom properties
a_array_concat = zeros(2, N, 3, 3);
b_array_concat = zeros(2, N, 3, 3);

% Create ensemble of connectomes with EDR, hub, and richclubs

for j = 1:3
    for i = 1:length(lambda_array)
        rng(1, "twister");
        lambda = lambda_array(i);
        if j == 1
            % Do not account for periodic conditions for aesthetic purposes
            isperiodic = 0;
            [a_array1, b_array1] = generate_connectome_edr(a_array, b_array, lambda, topology, isperiodic);
        elseif j == 2    
            [a_array1, b_array1] = generate_connectome_hub(a_array, b_array, lambda, topology, hub_centres, hublength);
        else
            [a_array1, b_array1] = generate_connectome_richclub(a_array, b_array, lambda, topology, hub_centres, hublength);
        end
        a_array_concat(:, :, i, j) = a_array1; b_array_concat(:, :, i, j) = b_array1;
    end
end


% Create tiledlayout figure
fig = figure;
fig.Position = [100, 100, 1400, 175];
set(gca,'Color','white')
box on;
t = tiledlayout(1, 28, 'TileSpacing', 'none', 'Padding', 'compact');
lambda_names = ["lambda_e", "lambda_h", "lambda_r"];
property_names = ["Exponential Distance Rule", "Hub Specificity", "Rich-club Specificity"];
label_names = ["i.", "ii.", "iii."];

% Create schematics of default topology with EDR

for j = 1:3
    for i = 1:length(lambda_array)
        a_array1 = a_array_concat(:, :, i, j); b_array1 = b_array_concat(:, :, i, j);
    
        ax = nexttile(3*i + (j-1)*10 - 2, [1 2]);
        hold;
        set(gca, 'Color', 'white');
        ax.Box = "on";
        ax.LineWidth = 1;
    
        xlim([0 topology.L])
        ylim([0 topology.L])
        xticks([]);
        yticks([]);

        % Draw hub regions
        if j == 2 || j == 3
            for k = 1:num_hubs
                fill(hub_centres(k, 1) + hublength/2 * [-1 -1 1 1], hub_centres(k, 2) + hublength/2 * [-1 1 1 -1], 'k', 'FaceColor', '#D95319', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
            end
        end

        for k = 1:N
            quiver(a_array1(1, k), a_array1(2, k), ...
            b_array1(1, k) - a_array1(1, k), ...
            b_array1(2, k) - a_array1(2, k),...
            'Color', 'k', 'LineWidth', 0.5, ...
            'MaxHeadSize', 0.05 / norm(b_array1(:, k) - a_array1(:, k)), ...
            'Marker', '.', 'MarkerSize', 0.0001, ...
            'AutoScale','off');
        end
        
        %Labels and titles
        xlabel(append('$\', lambda_names(j), ' = ', num2str(lambda), '$'), 'Interpreter', 'latex', 'FontSize', 16);

        if i == 2
            t1 = title(property_names(j), 'Interpreter', 'latex', 'FontSize', 16);
            set(t1, 'Units', 'normalized');
            t1.Position(2) = 1.1;
        end
        if i == 1
            t2 = title(append('\textbf{', label_names(j), '}'), 'Interpreter', 'latex', 'FontSize', 24);
            ax.TitleHorizontalAlignment = 'left';
        end
        hold off;
    
    end
end

% Draw arrows
for i = 1:2
    for j = [3*i, 10 + 3*i, 20 + 3*i]
        nexttile(j);
        hold;
        quiver(0.35, 0.5, 0.3, 0, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 1, 'Marker', '.', 'MarkerSize', 0.0001, 'AutoScale', 'off');
        xlim([0 1])
        ylim([0 1])
        xticks([]);
        yticks([]);
        axis off;
    end
end

% save as schematic_nonrandom.tiff
exportgraphics(gcf, '7_schematic_nonrandom.tiff', 'Resolution', 300);
close(fig);

%% BOLD dissimilarity of random connectomes

clear; clc;

loadparam;

% Set number of FNPs
N_array = [10, 20, 50, 100];
max_N = max(N_array);

% Set number samples - Change for 100 for figures in paper
num_samples = 10;

% Set tolerance for computing bold response
tol = 1e-5;

% Sample Stimulus Position
rng(0, "twister");
stim_positions = topology.L * rand(2, num_samples);

% Record BOLD dissimilarity
bold_dissim_array = zeros(length(N_array), num_samples);

parfor k = 1:num_samples

    % Create duplicates for parallelization
    stim1 = stim;
    hetparam_het1 = hetparam_het;

    % Set Stimulus Position
    stim1.stimR = stim_positions(:, k);

    % Simulate homogeneous model
    ts_hom_bold = run_bold(topology, homparam, hetparam_hom, stim1, tol);

    bold_dissim_subarray = zeros(length(N_array), 1);

    % Create Random Connectome
    rng(k, "twister");
    a_array = zeros(2, max_N);
    b_array = zeros(2, max_N);
    for j = 1:max_N
        a = topology.L * rand(2, 1);
        b = topology.L * rand(2, 1);
        a_array(:, j) = a;
        b_array(:, j) = b;
    end

    % Compute BOLD dissimilarity for each N
    for i = 1:length(N_array)

        N = N_array(i);

        hetparam_het1.m = N;
        hetparam_het1.c = (homparam.r)^2 * ones(1, N);
        hetparam_het1.tau = zeros(1, N);
        hetparam_het1.a = a_array(:, 1:N);
        hetparam_het1.b = b_array(:, 1:N);

        ts_het_bold = run_bold(topology, homparam, hetparam_het1, stim1, tol);

        dissim_bold = pdist2(...
            ts_hom_bold(:)', ts_het_bold(:)', 'cosine');

        bold_dissim_subarray(i) = dissim_bold;
        disp([i, k, dissim_bold]);

    end

    bold_dissim_array(:, k) = bold_dissim_subarray;

end

save('dissim_bold_multilrcs.mat', "topology", "N_array", "num_samples", "bold_dissim_array");


%% BOLD dissimilarity versus EDR decay rate
% Fix number of shortcuts to 50

clear; clc;
 
loadparam;

% Set number of FNPs
N_array = [10, 20, 50, 100];
max_N = max(N_array);

% Set number samples - Change to 100 for figures in paper
num_samples = 10;

% Set tolerance for computing bold responses
tol = 1e-5;

% Sample Stimulus Position
rng(0, "twister");
stim_positions = topology.L * rand(2, num_samples);

% Set lambda (EDR decay rate)
lambda_array = 0.2:0.2:1;
 

% Compute BOLD dissimilarity
bold_dissim_array = zeros(length(lambda_array), length(N_array), num_samples);

parfor k = 1:num_samples

    % Set seed for worker
    rng(k, "twister");

    % Create duplicates for parallelization
    stim1 = stim;
    hetparam_het1 = hetparam_het;

    % Set Stimulus Position
    stim1.stimR = stim_positions(:, k);

    % Simulate homogeneous model
    ts_hom_bold = run_bold(topology, homparam, hetparam_hom, stim1, tol);

    bold_dissim_subarray = zeros(length(lambda_array), length(N_array));

    % Create Random Connectome
    a_array = zeros(2, max_N);
    b_array = zeros(2, max_N);
    for j = 1:max_N
        a = topology.L * rand(2, 1);
        b = topology.L * rand(2, 1);
        a_array(:, j) = a;
        b_array(:, j) = b;
    end

    % Create Nonrandom connectome for each lambda
    for i = 1:length(lambda_array)

        lambda = lambda_array(i);
        [a_array1, b_array1] = generate_connectome_edr(a_array, b_array, lambda, topology);
        
        % Compute BOLD dissimilarity for each N
        for j = 1:length(N_array)
    
            N = N_array(j);
    
            hetparam_het1.m = N;
            hetparam_het1.c = (homparam.r)^2 * ones(1, N);
            hetparam_het1.tau = zeros(1, N);
            hetparam_het1.a = a_array1(:, 1:N);
            hetparam_het1.b = b_array1(:, 1:N);
    
            ts_het_bold = run_bold(topology, homparam, hetparam_het1, stim1, tol);
    
            dissim_bold = pdist2(...
                ts_hom_bold(:)', ts_het_bold(:)', 'cosine');
        
            bold_dissim_subarray(i, j) = dissim_bold;
            disp(num2str([i, j, k, dissim_bold]));
    
        end
    
        bold_dissim_array(:, :, k) = bold_dissim_subarray;
    
    end

end
%
save('dissim_bold_edr.mat', "topology", "N_array", "lambda_array", "num_samples", "bold_dissim_array");
% 


%% BOLD dissimilarity versus level of hub specificity

clear; clc;
 
loadparam;

tol = 1e-5;

% Set number of FNPs
N_array = [10, 20, 50, 100];
max_N = max(N_array);

% Set number samples - Change to 100 for figures in paper
num_samples = 10;

% Set number of hubs, positions, and dimensions
num_hubs = 4;
hub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
hublength = topology.L / sqrt(34);

% Set lambda (level of hub specificity)
lambda_array = 0.2:0.2:1;

bold_dissim_array = zeros(length(lambda_array), length(N_array), num_samples);

% Sample Stimulus Position
rng(0, "twister");
stim_positions = topology.L * rand(2, num_samples);

parfor k = 1:num_samples

    % Set seed for worker
    rng(k, "twister");

    % Create duplicates for parallelization
    stim1 = stim;
    hetparam_het1 = hetparam_het;

    % Set Stimulus Position
    stim1.stimR = stim_positions(:, k);

    % Simulate homogeneous model
    ts_hom_bold = run_bold(topology, homparam, hetparam_hom, stim1, tol);

    bold_dissim_subarray = zeros(length(lambda_array), length(N_array));

    % Create Random Connectome
      
    a_array = zeros(2, max_N);
    b_array = zeros(2, max_N);
    
    for j = 1:max_N
        a = topology.L * rand(2, 1);
        b = topology.L * rand(2, 1);
        a_array(:, j) = a;
        b_array(:, j) = b;
    end
    
    % Create Nonrandom connectome for each lambda
    for i = 1:length(lambda_array)       
        
        lambda = lambda_array(i);
        [a_array1, b_array1] = generate_connectome_hub(a_array, b_array, lambda, topology, hub_centres, hublength);

        % Compute BOLD dissimilarity for each N
        for j = 1:length(N_array)
    
            N = N_array(j);
    
            hetparam_het1.m = N;
            hetparam_het1.c = (homparam.r)^2 * ones(1, N);
            hetparam_het1.tau = zeros(1, N);
            hetparam_het1.a = a_array1(:, 1:N);
            hetparam_het1.b = b_array1(:, 1:N);
    
            ts_het_bold = run_bold(topology, homparam, hetparam_het1, stim1, tol);
    
            dissim_bold = pdist2(...
                ts_hom_bold(:)', ts_het_bold(:)', 'cosine');
        
            bold_dissim_subarray(i, j) = dissim_bold;
            disp(num2str([i, j, k, dissim_bold]));
    
        end
    
        bold_dissim_array(:, :, k) = bold_dissim_subarray;
    
    end

end
%
save('dissim_bold_hub.mat', "topology", "N_array", "num_hubs", "hub_centres", "hublength", "lambda_array", "num_samples", "bold_dissim_array");
% 

%% BOLD dissimilarity versus level of rich-club specificity

clear; clc;
 
loadparam;

tol = 1e-5;

% Set number of FNPs
N_array = [10, 20, 50, 100];
max_N = max(N_array);

% Set number samples - Change to 100 for figures in paper
num_samples = 10;

% Set number of hubs, positions, and dimensions
num_hubs = 4;
hub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
hublength = topology.L / sqrt(34);

% Set lambda (level of rich-club specificity)
lambda_array = 0.2:0.2:1;

bold_dissim_array = zeros(length(lambda_array), length(N_array), num_samples);

% Sample Stimulus Position
rng(0, "twister");
stim_positions = topology.L * rand(2, num_samples);

parfor k = 1:num_samples

    % Set seed for worker
    rng(k, "twister");

    % Create duplicates for parallelization
    stim1 = stim;
    hetparam_het1 = hetparam_het;

    % Set Stimulus Position
    stim1.stimR = stim_positions(:, k);
    
    % Simulate homogeneous model
    ts_hom_bold = run_bold(topology, homparam, hetparam_hom, stim1, tol);
    
    bold_dissim_subarray = zeros(length(lambda_array), length(N_array));
    
    % Create Random Connectome
        
    a_array = zeros(2, max_N);
    b_array = zeros(2, max_N);
    
    for j = 1:max_N
        a = topology.L * rand(2, 1);
        b = topology.L * rand(2, 1);
        a_array(:, j) = a;
        b_array(:, j) = b;
    end
    
    % Create Nonrandom connectome for each lambda
    for i = 1:length(lambda_array)

        lambda = lambda_array(i);
        [a_array1, b_array1] = generate_connectome_richclub(a_array, b_array, lambda, topology, hub_centres, hublength);

        % Compute BOLD dissimilarity for each N  
        for j = 1:length(N_array)
    
            N = N_array(j);
    
            hetparam_het1.m = N;
            hetparam_het1.c = (homparam.r)^2 * ones(1, N);
            hetparam_het1.tau = zeros(1, N);
            hetparam_het1.a = a_array1(:, 1:N);
            hetparam_het1.b = b_array1(:, 1:N);
    
            ts_het_bold = run_bold(topology, homparam, hetparam_het1, stim1, tol);
    
            dissim_bold = pdist2(...
                ts_hom_bold(:)', ts_het_bold(:)', 'cosine');
        
            bold_dissim_subarray(i, j) = dissim_bold;
            disp(num2str([i, j, k, dissim_bold]));
    
        end
    
        bold_dissim_array(:, :, k) = bold_dissim_subarray;
    
    end

end
%
save('dissim_bold_core.mat', "topology", "N_array", "num_hubs", "hub_centres", "hublength", "lambda_array", "num_samples", "bold_dissim_array");
% 

%% Plot mean BOLD dissimilarity versus three topological control parameters

clear; clc;
loadparam;

fig = figure;
fig.Position = [100, 100, 1400, 500];
hold;
title(['Cosine Distance vs. Probability of Hub Incidence'])
set(gcf, 'Color', 'white')
color_array = get(gca,'colororder');


% Load random connectome statistics
load('dissim_bold_multilrcs.mat', 'bold_dissim_array')
bold_dissim_array0 = bold_dissim_array;
clear bold_dissim_array;
bold_dissim_array0 = reshape(bold_dissim_array0, [1 size(bold_dissim_array0)]);

t = tiledlayout(1, 28, 'TileSpacing', 'None', 'Padding', 'Compact');

ax = nexttile(10*1 - 9, [1 8]); hold;

load('dissim_bold_edr.mat', 'N_array', 'lambda_array', 'bold_dissim_array');

% Append random to nonrandom connectome statistics [random, edr]
bold_dissim_array = cat(1, bold_dissim_array0, bold_dissim_array);
lambda_array = [0 lambda_array];

mean_dissimilarity = squeeze(mean(bold_dissim_array, 3));
std_dissimilarity = squeeze(std(bold_dissim_array, 1, 3));

set(ax,'ColorOrderIndex',1);
for i = 1:length(N_array)
    h = errorbar(lambda_array, mean_dissimilarity(:, i), std_dissimilarity(:, i), 'o', 'MarkerSize', 0.000001, 'LineWidth', 1);
    alpha = 0.2;   
    % Set transparency (undocumented)
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
end
for i = 1:length(N_array)
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

load('dissim_bold_hub.mat', "N_array", 'lambda_array', 'bold_dissim_array');


% Append random to nonrandom connectome statistics [random, hub]
bold_dissim_array = cat(1, bold_dissim_array0, bold_dissim_array);
lambda_array = [0 lambda_array];

mean_dissimilarity = squeeze(mean(bold_dissim_array, 3));
std_dissimilarity = squeeze(std(bold_dissim_array, 1, 3));

set(ax,'ColorOrderIndex',1);
for i = 1:length(N_array)
    h = errorbar(lambda_array, mean_dissimilarity(:, i), std_dissimilarity(:, i), 'o', 'MarkerSize', 0.000001, 'LineWidth', 1);
    alpha = 0.2;   
    % Set transparency (undocumented)
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
end
for i = 1:length(N_array)
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

load('dissim_bold_core.mat', 'N_array', 'lambda_array', 'bold_dissim_array');

% Append random to nonrandom connectome statistics [random, rich-club]
bold_dissim_array = cat(1, bold_dissim_array0, bold_dissim_array);
lambda_array = [0 lambda_array];

mean_dissimilarity = squeeze(mean(bold_dissim_array, 3));
std_dissimilarity = squeeze(std(bold_dissim_array, 1, 3));

set(ax,'ColorOrderIndex',1);
for i = 1:length(N_array)
    h = errorbar(lambda_array, mean_dissimilarity(:, i), std_dissimilarity(:, i), 'o', 'MarkerSize', 0.000001, 'LineWidth', 1);
    alpha = 0.2;   
    % Set transparency (undocumented)
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
end
for i = 1:length(N_array)
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

Legend = cell(1, 2*length(N_array));
for i = 1:length(N_array)
    Legend{i} = "";
    Legend{length(N_array) + i} = append('$N = ', num2str(N_array(i)), '$');
end

l = legend(Legend, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 16);
l.Layout.Tile = 'south';
l.Orientation = 'horizontal';

% save as dissim_bold_nonrandom.tiff
exportgraphics(gcf, '7_dissim_bold_nonrandom.tiff', 'Resolution', 300);
close(fig);

%% BOLD dissimilarity of sampled random, hub and rich-club connectomes

clear; clc;

loadparam;

% Set number of fnps
N = 100;
num_models = 3;

% Set tolerance for computing bold dissimilarity
tol = 1e-5;

% Set hubs, positions and dimensions
num_hubs = 4;
hub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
hublength = topology.L / sqrt(34);

% Set Stimulus Positions - change num_samples_x to 40 for resolution in
% paper
num_samples_x = 8;
num_samples = num_samples_x.^2;
[stim_x, stim_y] = meshgrid(1:num_samples_x, 1:num_samples_x);
stim_positions = (topology.L / num_samples_x)*([stim_x(:), stim_y(:)]');

% Set lambda for nonrandom connectomes (hub + richclub)
lambda = 1;

% Create Random, Hub and Rich club Connectomes

rng(1, "twister");
a_array = zeros(2, N);
b_array = zeros(2, N);
for j = 1:N
    a = topology.L * rand(2, 1);
    b = topology.L * rand(2, 1);
    a_array(:, j) = a;
    b_array(:, j) = b;
end
hetparam_het_array = cell(1, num_models);
for iter = 1:num_models
    if iter == 1
        a_array1 = a_array; b_array1 = b_array;
    elseif iter == 2
        rng(1, "twister")
        [a_array1, b_array1] = generate_connectome_hub(a_array, b_array, lambda, topology, hub_centres, hublength);
    elseif iter == 3
        rng(1, "twister")
        [a_array1, b_array1] = generate_connectome_richclub(a_array, b_array, lambda, topology, hub_centres, hublength);       
    end
    hetparam_het1 = hetparam_het;
    hetparam_het1.m = N;
    hetparam_het1.c = (homparam.r)^2 * ones(1, N);
    hetparam_het1.tau = zeros(1, N);
    hetparam_het1.a = a_array1(:, 1:N);
    hetparam_het1.b = b_array1(:, 1:N);
    hetparam_het_array{iter} = hetparam_het1;
end

% Compute bold dissimilarity for all stimulus positions on grid

bold_dissim_array = zeros(2 + num_models, num_samples);
parfor k = 1:num_samples

    stim1 = stim;
    stim1.stimR = stim_positions(:, k);

    % Simulate homogeneous model
    ts_hom_bold = run_bold(topology, homparam, hetparam_hom, stim1, tol);
    bold_dissim_subarray = zeros(2 + num_models, 1);
    bold_dissim_subarray(1:2) = stim1.stimR;
    
    for iter = 1:num_models
        hetparam_het = hetparam_het_array{iter};
        ts_het_bold = run_bold(topology, homparam, hetparam_het, stim1, tol);
        bold_dissim_subarray(2 + iter) = pdist2(ts_hom_bold(:)', ts_het_bold(:)', 'cosine');
        disp([k, iter]);
    end

    bold_dissim_array(:, k) = bold_dissim_subarray;
    
end
%
save('bold_dissim_randomvshubclub.mat', "topology", "num_models", "num_hubs", "hub_centres", "hublength", "hetparam_het_array", "num_samples", "num_samples_x", "bold_dissim_array");

%% Plot spatial distribution of BOLD dissimilarity of sampled random, hub and rich-club connectomes

clear; clc;
loadparam;

Colormap = [linspace(1, 0, 256)' linspace(1, 0, 256)', linspace(1, 1, 256)'];

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

load('bold_dissim_randomvshubclub.mat');

[xq,yq] = ndgrid(dx*(1:topology.Nx));
[stim_x, stim_y] = ndgrid(linspace(0, topology.L, num_samples_x + 1), linspace(0, topology.L, num_samples_x + 1));
dissim0 = reshape(bold_dissim_array(3:end, :, :), [num_models, num_samples_x, num_samples_x]);

% Create figure and properties
fig = figure;
fig.Position = [100 100 760 450];
num_subcols = 8;
num_subrows = 3;
tot_subcols = 2 + num_models*num_subcols;
tot_subrows = 1 + 2*num_subrows;
t = tiledlayout(tot_subrows, tot_subcols, 'TileSpacing', 'tight', 'Padding', 'compact');
modeltitles = ["Random", "Hub ($\lambda_h = 1$)", "Rich-club ($\lambda_r = 1$)"];
connectivitytitles = "\textbf{" + ["i", "ii", "iii"] + ".}";
dissimtitles = "\textbf{" + ["iv", "v", "vi"] + ".}";

% Title connectivites
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

% Draw connectivities
for iter = 1:num_models
    ax = nexttile(tot_subcols + 3 + (iter - 1)*(num_subcols), [num_subrows, num_subcols]);
    hold on;
    box on;
    ax.LineWidth = 1;
    ax.Color = 'k';
    set(gca, 'Color', 'white');
    hetparam_het = hetparam_het_array{iter};
    % Draw hub region
    if iter == 2 || iter == 3
        for k = 1:num_hubs
            fill(hub_centres(k, 1) + hublength/2 * [-1 -1 1 1], hub_centres(k, 2) + hublength/2 * [-1 1 1 -1], 'k', 'FaceColor', '#D95319', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
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

% Label dissimilarity for each connectivity
ax = nexttile((num_subrows + 1)*tot_subcols + 1, [num_subrows, 2]);
axis off;
text(0, 0.5, ['$C_z$ by stimulus', newline, 'position'], 'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90, 'Units', 'normalized')

% Plot heatmaps of BOLD dissimilarities
for iter = 1:num_models
    ax = nexttile((num_subrows + 1)*tot_subcols + 3 + (iter - 1)*(num_subcols), [num_subrows, num_subcols]);
    hold on;
    box on;
    hetparam_het = hetparam_het_array{iter};
    dissim = squeeze(dissim0(iter, :, :));
    % Create boundary conditions, introduce data points with x = 0 OR y = 0
    dissim = [dissim(end, end) dissim(end, :);
        dissim(:, end) dissim];
    imagesc((topology.L / num_samples_x)*(0:num_samples_x),(topology.L / num_samples_x)*(0:num_samples_x),dissim); 
    for m = 1:hetparam_het.m
        scatter(hetparam_het.a(1, m), hetparam_het.a(2, m), 5, 'k', 'filled');
    end
    clim([0  max(bold_dissim_array(3:end, :, :), [], 'all')]); 
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

% save as bold_dissim_randomvshubclub_spatialplot.svg
exportgraphics(gcf, '7_bold_dissim_randomvshubclub_spatialplot.tiff', 'Resolution', 300);
close(fig);

%% Distribution of BOLD dissimilarity of sampled random, hub and rich-club connectomes

clear; clc;
loadparam;

% Load dissimilarity data - throw out stimulus positions
load('bold_dissim_randomvshubclub.mat');
bold_dissim_array = bold_dissim_array(3:end, :);

% Create figure
fig = figure;
fig.Position = [100 100 600 430];
t = tiledlayout(1, 6, 'TileSpacing', 'tight', 'Padding', 'compact');

% Label plot
ax = nexttile(1, [1 1]);
axis off;
text(0, 1, '\textbf{vii.}', 'Interpreter', 'latex', 'FontSize', 24, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Units', 'normalized')

% Plot distributions
ax = nexttile(2, [1 5]);
hold;
colors = get(gca, 'ColorOrder');  % Retrieve the default color order
% Draw mean of distributions
for i = 1:num_models
    xline(mean(bold_dissim_array(i, :)), 'Color', colors(i, :), 'LineWidth', 2);
end
% Plot densities
for i = 1:num_models
    data = bold_dissim_array(i, :);
    [f, xi] = ksdensity(data, 'Support', 'positive');
    fill(xi, f, colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

xticks([0:0.05:0.5]); xlim([0, 0.15]);
ylabel("Density", 'Interpreter', 'latex', 'Rotation', 90);
xlabel('$C_z$', 'Interpreter', 'latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.XAxis.LabelFontSizeMultiplier  = 2;
ax.YAxis.LabelFontSizeMultiplier  = 1.5;

% Legend
Legend = {"Random", "Hub ($\lambda_h = 1$)", "Rich-club ($\lambda_r = 1$)"};
l = legend(Legend, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 12);

% save as dissim_bold_randomhubclub_distribution.svg
exportgraphics(gcf, '7_dissim_bold_randomhubclub_distribution.tiff', 'Resolution', 300);
close(fig);