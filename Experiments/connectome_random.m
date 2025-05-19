%% Create schematics of multiple FNPs

clear; clc;

loadparam;

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

num_lrcs_array = [10, 20, 50, 100];
max_num_lrcs = max(num_lrcs_array);

num_subcols = 9;
num_subcols_gap = 7;

f = figure;
f.Position = [100, 100, 1400, 160];
set(gca,'Color','white')
box on;
t = tiledlayout(1, length(num_lrcs_array)*(num_subcols+2*num_subcols_gap), 'TileSpacing', 'none', 'Padding', 'compact');

% Create Ensemble of Connectomes

rng(1, "twister");

a_array = zeros(2, max_num_lrcs);
b_array = zeros(2, max_num_lrcs);

for j = 1:max_num_lrcs
    a = topology.L * rand(2, 1);
    b = topology.L * rand(2, 1);
    a_array(:, j) = a;
    b_array(:, j) = b;
end

for j = 1:length(num_lrcs_array)

    num_lrcs = num_lrcs_array(j);

    ax = nexttile(1 + num_subcols_gap + (num_subcols+2*num_subcols_gap)*(j - 1), [1 num_subcols]);
    hold;
    set(gca, 'Color', 'white');
    ax.Box = "on";
    ax.LineWidth = 1;

    xlim([0 topology.L])
    ylim([0 topology.L])
    xticks([]);
    yticks([]);
    for k = 1:num_lrcs
        quiver(a_array(1, k), a_array(2, k), ...
        b_array(1, k) - a_array(1, k), ...
        b_array(2, k) - a_array(2, k),...
        'Color', 'k', 'LineWidth', 0.5, ...
        'MaxHeadSize', 0.05 / norm(a_array(:, k) - b_array(:, k)), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end

    xlabel(append('$N = ', num2str(num_lrcs), '$'), 'Interpreter', 'latex', 'FontSize', 16);

    hold off;

end

for j = 1:length(num_lrcs_array)-1
    nexttile(1 + num_subcols + num_subcols_gap + (num_subcols+2*num_subcols_gap)*(j - 1), [1, 2*num_subcols_gap]);
    hold;
    quiver(0.4, 0.5, 0.2, 0, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 1, 'Marker', '.', 'MarkerSize', 0.0001, 'AutoScale', 'off');
    xlim([0 1])
    ylim([0 1])
    xticks([]);
    yticks([]);
    axis off;
end

% save as schematicmultifnps.svg
print(gcf, 'schematicmultifnps.svg', '-dsvg');
exportgraphics(gcf, 'schematicmultifnps.tiff', 'Resolution', 300);

%% Compute dissimiliarity curves of RNG models over time

clear; clc;

loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

num_fnps_array = [10, 20, 50, 100];
max_num_lrcs = max(num_fnps_array);
num_samples = 100;

% Record total length used and times to thresholds
dissimcurve_array = zeros(topology.Nt, length(num_fnps_array), num_samples);

% Sample Stimulus Position
rng(0, "twister");
stim_positions = topology.L * rand(2, num_samples);

parfor k = 1:num_samples

    stim1 = stim;
    hetparam_het1 = hetparam_het;

    % % Set Stimulus Position
    % stim1.stimR = stim_positions(:, k);

    % % Simulate homogeneous model
    % ts_hom = run_periodic(topology, homparam, hetparam_hom, stim1);

    dissimcurve_subarray = zeros(topology.Nt, length(num_fnps_array));

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

    % Set Stimulus Position as source of the first FNP
    stim1.stimR = a_array(:, 1);

    % Simulate homogeneous model
    ts_hom = run_periodic(topology, homparam, hetparam_hom, stim1);

    for i = 1:length(num_fnps_array)

        num_lrcs = num_fnps_array(i);

        hetparam_het1.m = num_lrcs;
        hetparam_het1.c = (homparam.r)^2 * ones(1, num_lrcs);
        hetparam_het1.tau = zeros(1, num_lrcs);
        hetparam_het1.a = a_array(:, 1:num_lrcs);
        hetparam_het1.b = b_array(:, 1:num_lrcs);

        ts_het = run_periodic(topology, homparam, hetparam_het1, stim1);

        for n = 1:topology.Nt
            dissimcurve_subarray(n, i) = pdist2(...
                reshape(ts_hom(:, :, n), [], 1)', reshape(ts_het(:, :, n), [], 1)', 'cosine');
        end

        disp(num2str([i, k]))

    end

    dissimcurve_array(:, :, k) = dissimcurve_subarray;
    
end

save('dissimcurveensemble_random.mat', "dissimcurve_array", "num_fnps_array", "num_samples", "topology");

% save data as dissimovertime.mat

%% Plot dissimilarity curve across ensemble with multiple randomly positioned FNPs

clear; clc;
loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

time_array = dt * (1:topology.Nt);

% Plot divergence of default model and random models on same plot

f = figure;
f.Position = [100, 100, 1400, 500];

% num_subcols = 9;
% t = tiledlayout(1, 3*num_subcols + 1);

% ax = nexttile(t, [1 2*num_subcols]);
ax = gca;
hold on;

load('dissimcurveensemble_random.mat');

meandissimcurve = squeeze(mean(dissimcurve_array, 3));

for i = 1:length(num_fnps_array)
    for j = 1:num_samples
        currentColor = get(gca, 'ColorOrder');
        currentColor = currentColor(i, :);
        mixedColor = 0.3*currentColor + 0.7*[1 1 1];
        plot(1000*(time_array - stim.stimt), dissimcurve_array(:, i, j), 'Color', mixedColor, 'LineWidth', 0.1);
    end
end

set(ax,'ColorOrderIndex',1);
for i = 1:length(num_fnps_array)
    plot(1000*(time_array - stim.stimt), meandissimcurve(:, i), 'LineWidth', 2);
end

h = get(gca,'Children');
set(gca,'Children',[h(4); h(3); h(2); h(1); h(5:404)]);

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

Legend = cell(1, (num_samples + 1)*length(num_fnps_array));
for i = 1:num_samples*length(num_fnps_array)
    Legend{i} = "";
end
for i = 1:length(num_fnps_array)
    Legend{num_samples*length(num_fnps_array) + i} = append('$N = ', num2str(num_fnps_array(i)), '$');
end
l = legend(Legend, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 16);
% l.Layout.Tile = 'eastoutside';
l.Orientation = 'vertical';

hold off;

% save as dissimcurveboldensemble_multifnps.svg
print(gcf, 'dissimcurveensemble_multifnps.svg', '-dsvg');
exportgraphics(gcf, 'dissimcurveensemble_multifnps.tiff', 'Resolution', 300);