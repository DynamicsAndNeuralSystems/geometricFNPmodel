%% Create schematics of multiple FNPs

clear; clc;

loadparam;

% Number of FNPs
N_array = [10, 20, 50, 100];
max_N = max(N_array);

% Create tiledlayout figure
fig = figure;
fig.Position = [100, 100, 1400, 160];
set(gca,'Color','white')
% box on;
num_subcols = 9;
num_subcols_gap = 7;
t = tiledlayout(1, length(N_array)*(num_subcols+2*num_subcols_gap), 'TileSpacing', 'none', 'Padding', 'compact');

% Create Ensemble of Connectomes containing max specified number of FNPs

% Set seed
rng(1, "twister");

a_array = zeros(2, max_N);
b_array = zeros(2, max_N);

for j = 1:max_N
    a = topology.L * rand(2, 1);
    b = topology.L * rand(2, 1);
    a_array(:, j) = a;
    b_array(:, j) = b;
end

% Plot connectome containing specified number of FNPs. 
% For N FNPs, plot the first N sampled FNPs above
for j = 1:length(N_array)

    N = N_array(j);

    ax = nexttile(1 + num_subcols_gap + (num_subcols+2*num_subcols_gap)*(j - 1), [1 num_subcols]);
    hold;
    set(gca, 'Color', 'white');
    ax.Box = "on";
    ax.LineWidth = 1;

    xlim([0 topology.L])
    ylim([0 topology.L])
    xticks([]);
    yticks([]);
    for k = 1:N
        quiver(a_array(1, k), a_array(2, k), ...
        b_array(1, k) - a_array(1, k), ...
        b_array(2, k) - a_array(2, k),...
        'Color', 'k', 'LineWidth', 0.5, ...
        'MaxHeadSize', 0.05 / norm(a_array(:, k) - b_array(:, k)), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end

    xlabel(append('$N = ', num2str(N), '$'), 'Interpreter', 'latex', 'FontSize', 16);

    hold off;

end

% Plot arrow between plots
for j = 1:length(N_array)-1
    nexttile(1 + num_subcols + num_subcols_gap + (num_subcols+2*num_subcols_gap)*(j - 1), [1, 2*num_subcols_gap]);
    hold;
    quiver(0.4, 0.5, 0.2, 0, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 1, 'Marker', '.', 'MarkerSize', 0.0001, 'AutoScale', 'off');
    xlim([0 1])
    ylim([0 1])
    xticks([]);
    yticks([]);
    axis off;
end

% save as schematicmultifnps.tiff
exportgraphics(gcf, '6_schematicmultifnps.tiff', 'Resolution', 300);
close(fig);

%% Compute dissimiliarity curves of RNG models over time

clear; clc;

loadparam;

N_array = [10, 20, 50, 100];
max_N = max(N_array);

num_samples = 10; % change to 100 for used in figure

% Record total length used and times to thresholds
dissimcurve_array = zeros(topology.Nt, length(N_array), num_samples);

% Sample Stimulus Position
rng(0, "twister");
stim_positions = topology.L * rand(2, num_samples);

parfor k = 1:num_samples

    stim1 = stim;
    hetparam_het1 = hetparam_het;

    dissimcurve_subarray = zeros(topology.Nt, length(N_array));

    % Create Ensemble of Connectomes

    rng(k, "twister");

    a_array = zeros(2, max_N);
    b_array = zeros(2, max_N);

    for j = 1:max_N
        a = topology.L * rand(2, 1);
        b = topology.L * rand(2, 1);
        a_array(:, j) = a;
        b_array(:, j) = b;
    end

    % Set Stimulus Position as source of the first FNP
    stim1.stimR = a_array(:, 1);

    % Simulate homogeneous model
    ts_hom = run_periodic(topology, homparam, hetparam_hom, stim1);

    for i = 1:length(N_array)

        N = N_array(i);

        hetparam_het1.m = N;
        hetparam_het1.c = (homparam.r)^2 * ones(1, N);
        hetparam_het1.tau = zeros(1, N);
        hetparam_het1.a = a_array(:, 1:N);
        hetparam_het1.b = b_array(:, 1:N);

        ts_het = run_periodic(topology, homparam, hetparam_het1, stim1);

        for n = 1:topology.Nt
            dissimcurve_subarray(n, i) = pdist2(...
                reshape(ts_hom(:, :, n), [], 1)', reshape(ts_het(:, :, n), [], 1)', 'cosine');
        end

        disp(num2str([i, k]))

    end

    dissimcurve_array(:, :, k) = dissimcurve_subarray;
    
end

save('dissim_curve_random.mat', "dissimcurve_array", "N_array", "num_samples", "topology");

% save data as dissimovertime.mat

%% Plot dissimilarity curve across ensemble with multiple randomly positioned FNPs

clear; clc;
% Load model parameters
loadparam;

% Compute timestep
dt = topology.T / topology.Nt;

% Load simulated data
load('dissim_curve_random.mat');
meandissimcurve = squeeze(mean(dissimcurve_array, 3));

% Create figure
fig = figure;
fig.Position = [100, 100, 1400, 500];

ax = gca;
hold on;

% Plot dissimilarity curves
time_array = dt * (1:topology.Nt);
for i = 1:length(N_array)
    for j = 1:num_samples
        currentColor = get(gca, 'ColorOrder');
        currentColor = currentColor(i, :);
        mixedColor = 0.3*currentColor + 0.7*[1 1 1];
        plot(1000*(time_array - stim.stimt), dissimcurve_array(:, i, j), 'Color', mixedColor, 'LineWidth', 0.1);
    end
end
% Plot mean dissimilarity curve
set(ax,'ColorOrderIndex',1);
for i = 1:length(N_array)
    plot(1000*(time_array - stim.stimt), meandissimcurve(:, i), 'LineWidth', 2);
end

% Reorder layering of curves, put mean curves on top in with least N value curve on top
h = get(gca,'Children');
set(gca, 'Children', [h(length(N_array):-1:1); h(length(N_array)+1:length(h))]);

xlim([-5, Inf]);
xticks([0:10:100]);
yticks([0:0.1:1]);
xlabel("$t \ (\mathrm{ms})$", 'Interpreter', 'latex');
ylabel('$C_{\phi}$', 'Interpreter', 'latex', 'Rotation', 0);
ax.TickLabelInterpreter = 'latex';
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.XAxis.LabelFontSizeMultiplier  = 1.5;
ax.YAxis.LabelFontSizeMultiplier  = 2;

hold off;

% Legend
Legend = cell(1, (num_samples + 1)*length(N_array));
for i = 1:num_samples*length(N_array)
    Legend{i} = "";
end
for i = 1:length(N_array)
    Legend{num_samples*length(N_array) + i} = append('$N = ', num2str(N_array(i)), '$');
end
l = legend(Legend, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 16);
l.Orientation = 'vertical';


% save as dissim_curve_multifnps.svg
exportgraphics(gcf, '6_dissim_curve_multifnps.tiff', 'Resolution', 300);
close(fig);