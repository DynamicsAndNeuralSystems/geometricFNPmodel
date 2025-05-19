%% Stimulus evoked response of Model I, and Model II with 50 FNPs

clear; clc;

loadparam;

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

% Create discrete connectome of 50 FNPs for Model II
num_fnps = 50;
rng(1, "twister");
a_array = zeros(2, num_fnps);
b_array = zeros(2, num_fnps);
for j = 1:num_fnps   
    a = topology.L*rand(2, 1);
    b = topology.L*rand(2, 1);
    a_array(:, j) = a;
    b_array(:, j) = b;    
end

% Create FNP struct for Model II
hetparam_het.m = num_fnps;
hetparam_het.c = (homparam.r)^2 * ones(1, num_fnps);
hetparam_het.tau = zeros(1, num_fnps);
hetparam_het.a = a_array;
hetparam_het.b = b_array;

hetparam_array = {hetparam_hom, hetparam_het};

% Simulate evoked response
ts_full = zeros(topology.Nx, topology.Nx, topology.Nt, 2);
ts_full(:, :, :, 1) = run_periodic(topology,homparam,hetparam_hom,stim);
ts_full(:, :, :, 2) = run_periodic(topology,homparam,hetparam_het,stim);
%
maxval = 650;
timepoints = 0:0.002:0.010;
num_snapshots = length(timepoints);

Colormap = [linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];

f = figure;
f.Position = [100, 100, 1400, 410];

num_subcols = 5;
num_subrows = 10;
tot_subcols = 4+num_subcols*(1+num_snapshots);
tot_subrows = 2+num_subrows*2;

t = tiledlayout(tot_subrows,tot_subcols,"TileSpacing", "tight", "Padding", "compact");

modeltitles = "\textbf{" + ["Geometric", "Hybrid"] + "}";
stimulustitles = "\textbf{" + ["i", "ii"] + ".}";
responsetitles = "\textbf{" + ["iii", "iv"] + ".}";
coltitles = strings(1, 1+num_snapshots);
coltitles(1) = "Connectivity";
coltitles(2) = append("$\phi(\mathbf{r}, t = ", num2str(1000*timepoints(1)), "\ \mathrm{ms})$");
for i = 2:num_snapshots
    coltitles(i+1) = append("$t = ", num2str(1000*timepoints(i)), "\ \mathrm{ms}$");
end

for iter = 1:2
    hetparam_het = hetparam_array{iter};
    ts = squeeze(ts_full(:, :, :, iter));
    % Label Models
    nexttile(1 + (iter - 1)*(2+num_subrows)*tot_subcols, [num_subrows, 1]);
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(0.5, 0.5, modeltitles(iter), 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 16, 'Units', 'normalized', 'Rotation', 90);
    % % Label Connectivity Subplot
    % nexttile(2 + (iter - 1)*(2+num_subrows)*tot_subcols, [num_subrows, 1])
    % axis off;
    % set(gca, 'Color', 'none');
    % xticks([]); yticks([]);
    % text(1, 1, stimulustitles(iter), 'Interpreter', 'latex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 24, 'Units', 'normalized');
    % Plot Connectivity Subplot
    ax = nexttile(3 + (iter - 1)*(2+num_subrows)*tot_subcols, [num_subrows, num_subcols]);
    hold on;
    ax.Box = "on";
    ax.LineWidth = 2;
    set(gca, 'color', 'white');
    % Draw FNPs in Model II
    if iter == 2
        for k = 1:hetparam_het.m
            quiver(hetparam_het.a(1, k), hetparam_het.a(2, k), ...
            hetparam_het.b(1, k) - hetparam_het.a(1, k), ...
            hetparam_het.b(2, k) - hetparam_het.a(2, k),...
            'Color', 'k', 'LineWidth', 1, ...
            'MaxHeadSize', 0.05 / norm(hetparam_het.a(:, k) - hetparam_het.b(:, k)), ...
            'Marker', '.', 'MarkerSize', 0.0001, ...
            'AutoScale','off');
        end
    end
    view(0,90);
    xlim([0, topology.L + dx]);
    ylim([0, topology.L + dx]);
    xticks([]);
    yticks([]);
    % Label stimulus
    scatter(stim.stimR(1), stim.stimR(2), 15, 'r', 'filled');
    if iter == 1
        text(stim.stimR(1) + 0.03, stim.stimR(2) - 0.03, "$\mathbf{r}_{os}$", 'Interpreter', 'latex', 'Color', 'r', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 15);
        text(topology.L - 0.01, 0, "$\Omega$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontSize', 15);
    end
    % Title Connectivity Subplot
    if iter == 1
        title(coltitles(1), 'Interpreter', 'latex', 'FontSize', 16);
    end
    hold off;
    % % Label Response Subplot
    % nexttile(3 + num_subcols + (iter - 1)*(2+num_subrows)*tot_subcols, [num_subrows, 2])
    % axis off;
    % set(gca, 'Color', 'none');
    % xticks([]); yticks([]);
    % text(1, 1, responsetitles(iter), 'Interpreter', 'latex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 24, 'Units', 'normalized');
    % Plot Response Subplot
    for i = 1:num_snapshots
        ax = nexttile(5 + num_subcols + (iter - 1)*(2+num_subrows)*tot_subcols + (i - 1)*num_subcols, [num_subrows, num_subcols]);
        hold on;
        ax.Box = "on";
        ax.LineWidth = 1;
        imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), ts(:,:,ceil((stim.stimt + timepoints(i)) / dt))');
        set(ax,'YDir','normal');
        if iter == 1 && i == 1
            text(topology.L - 0.01, 0, "$\Omega$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontSize', 15);
        end
        if i > 1
            for k = 1:hetparam_het.m
                if ts(ceil(hetparam_het.a(1, k)/dx),ceil(hetparam_het.a(2, k)/dx),ceil((stim.stimt + timepoints(i))/dt)) >= 0.2*maxval && ...
                        ts(ceil(hetparam_het.b(1, k)/dx),ceil(hetparam_het.b(2, k)/dx),ceil((stim.stimt + timepoints(i))/dt)) >= 0.2*maxval && ...
                        ~(ts(ceil(hetparam_het.a(1, k)/dx),ceil(hetparam_het.a(2, k)/dx),ceil((stim.stimt + timepoints(i-1))/dt)) >= 0.2*maxval && ...
                        ts(ceil(hetparam_het.b(1, k)/dx),ceil(hetparam_het.b(2, k)/dx),ceil((stim.stimt + timepoints(i-1))/dt)) >= 0.2*maxval)
                    quiver(hetparam_het.a(1, k), hetparam_het.a(2, k), ...
                    hetparam_het.b(1, k) - hetparam_het.a(1, k), ...
                    hetparam_het.b(2, k) - hetparam_het.a(2, k),...
                    'Color', 'k', 'LineWidth', 1, ...
                    'MaxHeadSize', 0.05 / norm(hetparam_het.a(:, k) - hetparam_het.b(:, k)), ...
                    'Marker', '.', 'MarkerSize', 0.0001, ...
                    'AutoScale','off');
                end
            end
        end
        colormap(ax, Colormap);
        shading flat;
        clim([0 maxval]);
        view(0,90);
        xlim([0, topology.L + dx]);
        ylim([0, topology.L + dx]);
        xticks([]);
        yticks([]);
        % Title Heatmaps
        if iter == 1
            title(coltitles(i+1), 'Interpreter', 'latex', 'FontSize', 16);
        end
        hold off;
    end
end

% Colorbar
cb = colorbar; 
cb.Layout.Tile = 'east'; 
cb.FontSize = 15;
cb.TickLabelInterpreter = 'latex';
cbt = title(cb, "$\phi(\mathbf{r}, t)$", 'Interpreter', 'latex');

% Draw borders around connectivity and response subplots
for iter = 1:2
    ax = nexttile(5 + num_subcols + (iter - 1)*(2+num_subrows)*tot_subcols + (1 - 1)*num_subcols, [num_subrows, num_subcols]);
    pos1 = get(ax, 'Position');
    ax = nexttile(5 + num_subcols + (iter - 1)*(2+num_subrows)*tot_subcols + (num_snapshots - 1)*num_subcols, [num_subrows, num_subcols]);
    pos2 = get(ax, 'Position');
    annotation('rectangle', [pos1(1), pos1(2), pos2(1) + pos2(3) - pos1(1), pos2(2) + pos2(4) - pos1(2)], ...
    'LineWidth', 2, 'Color', 'k');
end

% save as proofofprinciple.svg
print(gcf, 'proofofprinciple.svg', '-dsvg');
exportgraphics(gcf, 'proofofprinciple.tiff', 'Resolution', 300);