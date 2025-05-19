%% 4A. Visualize default model with nonoptimal stimulation with snapshots over time

clear; clc;
loadparam;

% Nonoptimal stimulation
stim.stimR = [0.25; 0.15];

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

ts_full = zeros(topology.Nx, topology.Nx, topology.Nt, 2);
ts_full(:, :, :, 1) = run_periodic(topology,homparam,hetparam_hom,stim);
ts_full(:, :, :, 2) = run_periodic(topology,homparam,hetparam_het,stim);

maxval = 650;
timepoints = [0, 0.005, 0.01, 0.02, 0.03, 0.05];
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
coltitles(2) = append("$\phi(\mathbf{r}, t = 0 \ \mathrm{ms})$");
for i = 2:num_snapshots
    coltitles(i+1) = append("$t = ", num2str(1000*timepoints(i)), "\ \mathrm{ms}$");
end

for iter = 1:2
    ts = squeeze(ts_full(:, :, :, iter));
    % Label Connectivities
    nexttile(1 + (iter - 1)*(2+num_subrows)*tot_subcols, [num_subrows, 1]);
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(0.5, 0.5, modeltitles(iter), 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 16, 'Units', 'normalized', 'Rotation', 90);
    % Label Connectivity Subplot
    nexttile(2 + (iter - 1)*(2+num_subrows)*tot_subcols, [num_subrows, 1])
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(1, 1, stimulustitles(iter), 'Interpreter', 'latex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 24, 'Units', 'normalized');
    % Plot Connectivity Subplot
    ax = nexttile(3 + (iter - 1)*(2+num_subrows)*tot_subcols, [num_subrows, num_subcols]);
    hold on;
    ax.Box = "on";
    ax.LineWidth = 2;
    set(gca, 'color', 'white');
    % Draw FNP in Model II
    if iter == 2
        quiver(hetparam_het.a(1), hetparam_het.a(2), ...
        hetparam_het.b(1) - hetparam_het.a(1), ...
        hetparam_het.b(2) - hetparam_het.a(2),...
        'Color', 'k', 'LineWidth', 1, ...
        'MaxHeadSize', 0.05 / norm(hetparam_het.a - hetparam_het.b), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end
    view(0,90);
    xlim([0, topology.L + dx]);
    ylim([0, topology.L + dx]);
    xticks([]);
    yticks([]);
    % Label Stimulus
    scatter(stim.stimR(1), stim.stimR(2), 15, 'r', 'filled');
    if iter == 1
        text(topology.L - 0.01, 0, "$\Omega$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontSize', 15);
        text(stim.stimR(1) + 0.01, stim.stimR(2) - 0.02, "$\mathbf{r}_{os}$", 'Interpreter', 'latex', 'Color', 'r', 'HorizontalAlignment','left','VerticalAlignment','middle', 'FontSize', 15);
    elseif iter == 2
        text(hetparam_het.a(1) - 0.02, hetparam_het.a(2) - 0.02, "$\mathbf{p}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 15);
        text(hetparam_het.b(1) + 0.02, hetparam_het.b(2) + 0.02, "$\mathbf{q}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 15);            
    end
    % Title Connectivity Subplot
    if iter == 1
        title(coltitles(1), 'Interpreter', 'latex', 'FontSize', 16);
    end
    hold off;
    % Label Response Subplot
    nexttile(3 + num_subcols + (iter - 1)*(2+num_subrows)*tot_subcols, [num_subrows, 2])
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(1, 1, responsetitles(iter), 'Interpreter', 'latex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 24, 'Units', 'normalized');
    % Plot Response Subplot
    for i = 1:num_snapshots
        ax = nexttile(5 + num_subcols + (iter - 1)*(2+num_subrows)*tot_subcols + (i - 1)*num_subcols, [num_subrows, num_subcols]);
        hold on;
        ax.Box = "on";
        ax.LineWidth = 1;
        imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), ts(:,:,ceil((stim.stimt + timepoints(i)) / dt))');
        set(ax, 'YDir', 'normal')
        if iter == 1 && i == 1
            text(topology.L - 0.01, 0, "$\Omega$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontSize', 15);
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

% save as snapshotsdefaultnonoptimal.svg
print(gcf, 'snapshotsdefaultnonoptimal.svg', '-dsvg');
exportgraphics(gcf, 'snapshotsdefaultnonoptimal.tiff', 'Resolution', 300);

%% Plot positions of impulse stimulus

clear; clc;
loadparam;

f = figure;
f.Position = [0 100 310 295];

num_positions = 5;
positions = [linspace(hetparam_het.a(1), 0.25, num_positions); hetparam_het.a(2) * ones(1, num_positions)];
positiontexts = "\mathbf{r}_" + string(1:num_positions);

color_array = get(gca,'colororder');

ax = gca;
hold on;
ax.Box = "on";
ax.LineWidth = 2;
ax.TickLabelInterpreter = 'latex';
set(gca, 'color', 'white');
xticks([]);
yticks([]);
xlim([0 topology.L]);
ylim([0 topology.L]);
text(hetparam_het.a(1) - 0.02, hetparam_het.a(2) - 0.02, "$\mathbf{p}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 20);
text(hetparam_het.b(1) + 0.02, hetparam_het.b(2) + 0.02, "$\mathbf{q}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 20);            
scatter(hetparam_het.a(1), hetparam_het.a(2), 20, 'k', 'filled');
scatter(hetparam_het.b(1), hetparam_het.b(2), 20, 'k', 'filled');
for iter = 1:num_positions
    scatter(positions(1, iter), positions(2, iter), 20, color_array(iter, :), 'filled');
end
hold off;

% save as dissimcurvepositionlegend.svg
print(gcf, 'dissimcurvepositionlegend.svg', '-dsvg');
exportgraphics(gcf, 'dissimcurvepositionlegend.tiff', 'Resolution', 300);

%% Plot dissimilarity of default model over time, varying stimulus positions

clear; clc;
loadparam;

num_positions = 5;
stim_pos_array = [linspace(hetparam_het.a(1), 0.25, num_positions); hetparam_het.a(2) * ones(1, num_positions)];


dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

time_array = dt * (1:topology.Nt);

dissim_time = zeros(size(stim_pos_array, 5), topology.Nt);

for i = 1:size(stim_pos_array, 2)

    stim.stimR = stim_pos_array(:, i);

    ts_hom = run_periodic(topology,homparam,hetparam_hom,stim);
    ts_het = run_periodic(topology,homparam,hetparam_het,stim);
    
    for n = 1:topology.Nt
        dissim_time(i, n) = pdist2(...
            reshape(ts_hom(:, :, n), [], 1)', reshape(ts_het(:, :, n), [], 1)', 'cosine');
    end

end

f = figure;
f.Position = [100, 100, 1400, 500];
hold;

color_array = get(gca,'colororder');

for i = num_positions:-1:1
    plot(1000*(time_array - stim.stimt), dissim_time(i, :)', 'Color', color_array(i, :), 'LineWidth', 2);
end

ylim([0, 0.09]);
xlim([-5, Inf]);
xticks([0:10:40]);
yticks([0:0.02:1]);
xlabel("$t \ (\mathrm{ms})$", 'Interpreter', 'latex');
ylabel('$C_{\phi}$', 'Interpreter', 'latex', 'Rotation', 0)

ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.XAxis.LabelFontSizeMultiplier  = 1.5;
ax.YAxis.LabelFontSizeMultiplier  = 2;

% save as dissimcurveposition.svg
print(gcf, 'dissimcurveposition.svg', '-dsvg');
exportgraphics(gcf, 'dissimcurveposition.tiff', 'Resolution', 300);