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
fig.Position = [100, 100, 1390, 175];
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
        xlabel(append('$\', lambda_names(j), ' = ', num2str(lambda_array(i)), '$'), 'Interpreter', 'latex', 'FontSize', 16);

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

%% Max dissimilarity of random connectomes

clear; clc;

loadparam;

% Set number of FNPs
N_array = [10, 20, 50, 100];
max_N = max(N_array);

% Set number samples - Change for 100 for figures in paper
num_samples = 200;

% Sample Stimulus Position
rng(0, "twister");
stim_positions = topology.L * rand(2, num_samples);

% Record BOLD dissimilarity
max_dissim_array = zeros(length(N_array), num_samples);

parfor k = 1:num_samples

    % Create duplicates for parallelization
    stim1 = stim;
    hetparam_het1 = hetparam_het;

    % Set Stimulus Position
    stim1.stimR = stim_positions(:, k);

    % Simulate homogeneous model
    ts_hom = run_periodic(topology, homparam, hetparam_hom, stim1);

    % Initiate dissimcurve array
    dissimcurve = zeros(1, topology.Nt);

    max_dissim_subarray = zeros(length(N_array), 1);

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

        ts_het = run_periodic(topology, homparam, hetparam_het1, stim1);

        for n = 1:topology.Nt
            dissimcurve(n) = pdist2(...
                reshape(ts_hom(:, :, n), [], 1)', reshape(ts_het(:, :, n), [], 1)', 'cosine');
        end

        dissim_max = max(dissimcurve);
        max_dissim_subarray(i) = dissim_max;
        disp([i, k, dissim_max]);

    end

    max_dissim_array(:, k) = max_dissim_subarray;

end

save('dissim_max_multilrcs.mat', "topology", "N_array", "num_samples", "stim_positions", "max_dissim_array");


%% Max dissimilarity versus EDR decay rate

clear; clc;
 
loadparam;

% Set number of FNPs
N_array = [10, 20, 50, 100];
max_N = max(N_array);

% Set number samples - Change to 100 for figures in paper
num_samples = 200;

% Sample Stimulus Position
rng(0, "twister");
stim_positions = topology.L * rand(2, num_samples);

% Set lambda (EDR decay rate)
lambda_array = 0.2:0.2:1;
 

% Compute BOLD dissimilarity
max_dissim_array = zeros(length(lambda_array), length(N_array), num_samples);

parfor k = 1:num_samples

    % Set seed for worker
    rng(k, "twister");

    % Create duplicates for parallelization
    stim1 = stim;
    hetparam_het1 = hetparam_het;

    % Set Stimulus Position
    stim1.stimR = stim_positions(:, k);

    % Simulate homogeneous model
    ts_hom = run_periodic(topology, homparam, hetparam_hom, stim1);

    % Initiate dissimcurve array
    dissimcurve = zeros(1, topology.Nt);

    max_dissim_subarray = zeros(length(lambda_array), length(N_array));

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
    
            ts_het = run_periodic(topology, homparam, hetparam_het1, stim1);
    
            for n = 1:topology.Nt
                dissimcurve(n) = pdist2(...
                    reshape(ts_hom(:, :, n), [], 1)', reshape(ts_het(:, :, n), [], 1)', 'cosine');
            end
    
            dissim_max = max(dissimcurve);
            max_dissim_subarray(i, j) = dissim_max;
            disp(num2str([i, j, k, dissim_max]));
    
        end
    
        max_dissim_array(:, :, k) = max_dissim_subarray;
    
    end

end
%
save('dissim_max_edr.mat', "topology", "N_array", "lambda_array", "num_samples", "stim_positions", "max_dissim_array");
%

%% Max dissimilarity versus level of hub specificity

clear; clc;
 
loadparam;

% Set number of FNPs
N_array = [10, 20, 50, 100];
max_N = max(N_array);

% Set number samples - Change to 100 for figures in paper
num_samples = 200;

% Set number of hubs, positions, and dimensions
num_hubs = 4;
hub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
hublength = topology.L / sqrt(34);

% Set lambda (level of hub specificity)
lambda_array = 0.2:0.2:1;

max_dissim_array = zeros(length(lambda_array), length(N_array), num_samples);

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
    ts_hom = run_periodic(topology, homparam, hetparam_hom, stim1);

    % Initiate dissimcurve array
    dissimcurve = zeros(1, topology.Nt);
    
    max_dissim_subarray = zeros(length(lambda_array), length(N_array));

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
    
            ts_het = run_periodic(topology, homparam, hetparam_het1, stim1);
    
            for n = 1:topology.Nt
                dissimcurve(n) = pdist2(...
                    reshape(ts_hom(:, :, n), [], 1)', reshape(ts_het(:, :, n), [], 1)', 'cosine');
            end
    
            dissim_max = max(dissimcurve);
            max_dissim_subarray(i, j) = dissim_max;
            disp(num2str([i, j, k, dissim_max]));
    
        end
    
        max_dissim_array(:, :, k) = max_dissim_subarray;
    
    end

end
%
save('dissim_max_hub.mat', "topology", "N_array", "num_hubs", "hub_centres", "hublength", "lambda_array", "num_samples", "stim_positions", "max_dissim_array");
% 

%% Max dissimilarity versus level of rich-club specificity

clear; clc;
 
loadparam;

% Set number of FNPs
N_array = [10, 20, 50, 100];
max_N = max(N_array);

% Set number samples - Change to 100 for figures in paper
num_samples = 200;

% Set number of hubs, positions, and dimensions
num_hubs = 4;
hub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
hublength = topology.L / sqrt(34);

% Set lambda (level of rich-club specificity)
lambda_array = 0.2:0.2:1;

max_dissim_array = zeros(length(lambda_array), length(N_array), num_samples);

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
    ts_hom = run_periodic(topology, homparam, hetparam_hom, stim1);
    
    % Initiate dissimcurve array
    dissimcurve = zeros(1, topology.Nt);
    
    max_dissim_subarray = zeros(length(lambda_array), length(N_array));
    
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
    
            ts_het = run_periodic(topology, homparam, hetparam_het1, stim1);
    
            for n = 1:topology.Nt
                dissimcurve(n) = pdist2(...
                    reshape(ts_hom(:, :, n), [], 1)', reshape(ts_het(:, :, n), [], 1)', 'cosine');
            end
    
            dissim_max = max(dissimcurve);
            max_dissim_subarray(i, j) = dissim_max;
            disp(num2str([i, j, k, dissim_max]));
    
        end
    
        max_dissim_array(:, :, k) = max_dissim_subarray;
    
    end

end
%
save('dissim_max_core.mat', "topology", "N_array", "num_hubs", "hub_centres", "hublength", "lambda_array", "num_samples", "stim_positions", "max_dissim_array");
% 

%% Plot max dissimilarity versus three topological control parameters

clear; clc;
loadparam;

fig = figure;
fig.Position = [100, 100, 1400, 500];
hold;
title(['Cosine Distance vs. Probability of Hub Incidence'])
set(gcf, 'Color', 'white')
color_array = get(gca,'colororder');


% Load random connectome statistics
load('dissim_max_multilrcs.mat', 'max_dissim_array')
max_dissim_array0 = max_dissim_array;
clear max_dissim_array;
max_dissim_array0 = reshape(max_dissim_array0, [1 size(max_dissim_array0)]);

t = tiledlayout(1, 28, 'TileSpacing', 'None', 'Padding', 'Compact');

ax = nexttile(10*1 - 9, [1 8]); hold;

load('dissim_max_edr.mat', 'N_array', 'lambda_array', 'max_dissim_array');

% Append random to nonrandom connectome statistics [random, edr]
max_dissim_array = cat(1, max_dissim_array0, max_dissim_array);
lambda_array = [0 lambda_array];

mean_dissimilarity = squeeze(mean(max_dissim_array, 3));
std_dissimilarity = squeeze(std(max_dissim_array, 1, 3));

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
ylabel('$C_{\max}$', 'Interpreter', 'latex', 'Rotation', 0)
ax.TickLabelInterpreter = 'latex';
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.XAxis.LabelFontSizeMultiplier  = 1.5;
ax.YAxis.LabelFontSizeMultiplier  = 1.5;

xlim([-0.05*max(xlim) 1.05*max(xlim)])
yticks([0:0.1:0.4]); ylim([0 Inf]);
xticks(lambda_array)

title('\textbf{i.}', 'Interpreter', 'latex', 'FontSize', 24);
ax.TitleHorizontalAlignment = 'left';

hold off;

ax = nexttile(10*2 - 9, [1 8]); hold;

load('dissim_max_hub.mat', "N_array", 'lambda_array', 'max_dissim_array');


% Append random to nonrandom connectome statistics [random, hub]
max_dissim_array = cat(1, max_dissim_array0, max_dissim_array);
lambda_array = [0 lambda_array];

mean_dissimilarity = squeeze(mean(max_dissim_array, 3));
std_dissimilarity = squeeze(std(max_dissim_array, 1, 3));

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
yticks([0:0.1:0.4]); ylim([0 Inf]);
xticks(lambda_array);

title('\textbf{ii.}', 'Interpreter', 'latex', 'FontSize', 24);
ax.TitleHorizontalAlignment = 'left';

hold off;

ax = nexttile(10*3 - 9, [1 8]); hold;

load('dissim_max_core.mat', 'N_array', 'lambda_array', 'max_dissim_array');

% Append random to nonrandom connectome statistics [random, rich-club]
max_dissim_array = cat(1, max_dissim_array0, max_dissim_array);
lambda_array = [0 lambda_array];

mean_dissimilarity = squeeze(mean(max_dissim_array, 3));
std_dissimilarity = squeeze(std(max_dissim_array, 1, 3));

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
yticks([0:0.1:0.4]); ylim([0 Inf]);
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

% save as dissim_max_nonrandom.tiff
exportgraphics(gcf, '7_dissim_max_nonrandom.tiff', 'Resolution', 300);
close(fig);

%% Max dissimilarity of sampled random, hub and rich-club connectomes

clear; clc;

loadparam;

% Set number of fnps
N = 100;
num_models = 3;

% Set hubs, positions and dimensions
num_hubs = 4;
hub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
hublength = topology.L / sqrt(34);

% Set Stimulus Positions - change num_samples_x to 40 for resolution in
% paper
num_samples_x = 40;
num_samples = num_samples_x.^2;
[stim_x, stim_y] = meshgrid(1:num_samples_x, 1:num_samples_x);
stim_positions = (topology.L / num_samples_x)*([stim_x(:), stim_y(:)]');

% Set lambda for nonrandom connectomes (hub + richclub)
lambda = 1;

% Create Random, Hub and Rich club Connectomes

hetparam_het_array = cell(1, num_models);
for iter = 1:num_models
    rng(1, "twister");
    a_array = zeros(2, N);
    b_array = zeros(2, N);
    for j = 1:N
        a = topology.L * rand(2, 1);
        b = topology.L * rand(2, 1);
        a_array(:, j) = a;
        b_array(:, j) = b;
    end
    if iter == 1
        a_array1 = a_array; b_array1 = b_array;
    elseif iter == 2
        [a_array1, b_array1] = generate_connectome_hub(a_array, b_array, lambda, topology, hub_centres, hublength);
    elseif iter == 3
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

% Compute max dissimilarity for all stimulus positions on grid

max_dissim_array = zeros(2 + num_models, num_samples);
parfor k = 1:num_samples

    stim1 = stim;
    stim1.stimR = stim_positions(:, k);

    % Simulate homogeneous model
    ts_hom = run_periodic(topology, homparam, hetparam_hom, stim1);

    % Initiate dissimcurve array
    dissimcurve = zeros(1, topology.Nt);

    max_dissim_subarray = zeros(2 + num_models, 1);
    max_dissim_subarray(1:2) = stim1.stimR;
    
    for iter = 1:num_models
        hetparam_het = hetparam_het_array{iter};
        ts_het = run_periodic(topology, homparam, hetparam_het, stim1);
        for n = 1:topology.Nt
            dissimcurve(n) = pdist2(...
                reshape(ts_hom(:, :, n), [], 1)', reshape(ts_het(:, :, n), [], 1)', 'cosine');
        end
        max_dissim_subarray(2 + iter) = max(dissimcurve);
        disp([k, iter]);
    end

    max_dissim_array(:, k) = max_dissim_subarray;
    
end
%
save('max_dissim_randomvshubclub.mat', "topology", "num_models", "num_hubs", "hub_centres", "hublength", "hetparam_het_array", "num_samples", "num_samples_x", "max_dissim_array");

%% Plot spatial distribution of max dissimilarity of sampled random, hub and rich-club connectomes

clear; clc;
loadparam;

colors = [0 0 0; 1 0 0; 0 0.5 0];

Colormap = [linspace(1, 0, 256)' linspace(1, 0, 256)', linspace(1, 1, 256)'];

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

load('max_dissim_randomvshubclub.mat');

dissim0 = reshape(max_dissim_array(3:end, :, :), [num_models, num_samples_x, num_samples_x]);

% Create figure and properties
fig = figure;
fig.Position = [100 100 770 450];
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
text(0, 0.5, ['$C_{\max}$ by stimulus', newline, 'position'], 'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90, 'Units', 'normalized')

% Plot heatmaps of max dissimilarities
for iter = 1:num_models
    ax = nexttile((num_subrows + 1)*tot_subcols + 3 + (iter - 1)*(num_subcols), [num_subrows, num_subcols]);
    hold on;
    hetparam_het = hetparam_het_array{iter};
    dissim = squeeze(dissim0(iter, :, :));
    % Create boundary conditions, introduce data points with x = 0 OR y = 0
    dissim = [dissim(end, end) dissim(end, :);
        dissim(:, end) dissim];
    imagesc((topology.L / num_samples_x)*(0:num_samples_x),(topology.L / num_samples_x)*(0:num_samples_x),dissim); 
    for m = 1:hetparam_het.m
        scatter(hetparam_het.a(1, m), hetparam_het.a(2, m), 5, 'k', 'filled');
    end
    clim([0  0.9]); 
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
cb.Ticks = [0, 0.9];
cb.TickLabels = {'0', '0.9'};
ylabel(cb, '$C_{\max}$', 'Interpreter', 'latex', 'Rotation', 0, 'Position',[2.5 0.495 0]);

% Draw borders
for iter = 1:num_models
    ax = nexttile(tot_subcols + 3 + (iter - 1)*(num_subcols), [num_subrows, num_subcols]);
    hold on;
    pos = get(ax, 'Position');
    annotation('rectangle', [pos(1), pos(2), pos(1) + pos(3) - pos(1), pos(2) + pos(4) - pos(2)], ...
        'LineWidth', 1.5, 'Color', colors(iter, :));
    ax = nexttile((num_subrows + 1)*tot_subcols + 3 + (iter - 1)*(num_subcols), [num_subrows, num_subcols]);
    hold on;
    pos = get(ax, 'Position');
    annotation('rectangle', [pos(1), pos(2), pos(1) + pos(3) - pos(1), pos(2) + pos(4) - pos(2)], ...
        'LineWidth', 1.5, 'Color', colors(iter, :));
end

% save as max_dissim_randomvshubclub_spatialplot.svg
exportgraphics(gcf, '7_max_dissim_singlerandomvshubclub_spatialplot.tiff', 'Resolution', 300);
close(fig);

%% Distribution of max dissimilarity of ensemble of random, hub (lambda_h = 1) and rich-club (lambda_r = 1) connectomes

clear; clc;
loadparam;

num_models = 3;

% Load random dissimilarity data
load('dissim_max_multilrcs.mat');
max_dissim_array0 = max_dissim_array(end, :);
% Load hub and rich club dissimilarity data (lambda = 1)
load('dissim_max_hub.mat');
max_dissim_array = reshape(max_dissim_array(end, end, :), [1, num_samples]);
max_dissim_array0 = cat(1, max_dissim_array0, max_dissim_array);
load('dissim_max_core.mat');
max_dissim_array = reshape(max_dissim_array(end, end, :), [1, num_samples]);
max_dissim_array0 = cat(1, max_dissim_array0, max_dissim_array);

max_dissim_array = max_dissim_array0; clear max_dissim_array0;

% Sample Stimulus Position
rng(0, "twister");
stim_positions = topology.L * rand(2, num_samples);

% Compute dist between stim and nearest hub
dist_array = zeros(1, num_samples);
dist_subarray = zeros(1, num_hubs);
for i = 1:num_samples
    for j = 1:num_hubs
        distx = abs(stim_positions(1, i) - hub_centres(j, 1));
        disty = abs(stim_positions(2, i) - hub_centres(j, 2));
        distx = min(distx, topology.L - distx);
        disty = min(disty, topology.L - disty);
        dist_subarray(j) = norm([distx; disty], 2);  
    end
    dist_array(i) = min(dist_subarray);
end

% Create figure
fig = figure; hold;
fig.Position = [100 100 550 430];
colors = [0 0 0; 1 0 0; 0 0.5 0];
ax = gca;

% Plot dummy scatter for legend
for i = 2:num_models
    scatter(NaN, NaN, 5, 'filled', 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', 'none');
end

% Plot scatter of distance and max dissimilarity
for i = 2:num_models
    scatter(dist_array, max_dissim_array(i, :), 5, 'filled', 'MarkerFaceColor', 0.5*([1 1 1] + colors(i, :)), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7);
end

xlim([0 0.16]); ylim([0 0.01*ceil(max(max_dissim_array, [], 'all')/0.01)]);
xlabel(['Distance between stimulus', newline, 'and nearest hub'], 'Interpreter', 'latex');
ylabel('$C_{\max}$', 'Interpreter', 'latex', 'Rotation', 0);
ax.TickLabelInterpreter = 'latex';
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.XAxis.LabelFontSizeMultiplier  = 1;
ax.YAxis.LabelFontSizeMultiplier  = 1.5;

% Draw mean of random connectome distribution
yline(mean(max_dissim_array(1, :)), 'Color', colors(1, :), 'LineStyle', '--', 'LineWidth', 2, 'Label', 'Random', 'Interpreter', 'latex', 'FontSize', 12);

% Compute linear fit and annotate
for i = 2:num_models
    coefficients = polyfit(dist_array, log(max_dissim_array(i, :)), 1);
    plot(linspace(min(xlim), max(xlim), 100), exp(polyval(coefficients, (linspace(min(xlim), max(xlim), 100)))), 'Color', colors(i, :), 'LineStyle', '-', 'LineWidth', 2);
end

% Draw mean of distributions
for i = 2:num_models
    yline(mean(max_dissim_array(i, :)), 'Color', colors(i, :), 'LineStyle', '--', 'LineWidth', 1);
end

% Legend
Legend = {"Hub ($\lambda_h = 1$)", "Rich-club ($\lambda_r = 1$)"};
l = legend(Legend, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 12);

% save as max_dissim_randomvshubclub_spatialplot.tiff
exportgraphics(gcf, '7_max_dissim_randomhubclub.tiff', 'Resolution', 300);
close(fig);