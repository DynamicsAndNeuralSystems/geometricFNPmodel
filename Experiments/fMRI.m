%% Visualization of BOLD response as a sum of responses

clear; clc;
loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

Colormap = [linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];

timepoints = 0:0.005:0.01;
num_snapshots = length(timepoints);

num_subtiles = 3;

f = figure;
f.Position = [100, 100, 1550, 215];
t = tiledlayout(1,1+num_subtiles+1+(2*num_subtiles)*num_snapshots+num_subtiles,"TileSpacing","none","Padding","compact");
%title(t, 'Spatiotemporal Evoked Response')

figtitles = "\textbf{" + ["i", "ii"] + ".}";
coltitles = strings(1, num_snapshots);
coltitles(1) = append("$\phi(\mathbf{r}, t = ", num2str(1000*timepoints(1)), "\ \mathrm{ms})$");
for i = 2:num_snapshots
    coltitles(i) = append("$t = ", num2str(1000*timepoints(i)), "\ \mathrm{ms}$");
end


ts = run_periodic(topology,homparam,hetparam_hom,stim);
maxval = 600;

tol = 1e-5;
ts_bold = run_bold(topology,homparam,hetparam_hom,stim,tol);
maxvalbold = max(ts_bold, [], 'all');

ax = nexttile;
axis off;
set(gca, 'Color', 'none');
text(1, 1, figtitles(1), 'Interpreter', 'latex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 24, 'Units', 'normalized');

ax = nexttile([1 num_subtiles]);
hold on;
ax.Box = "on";
ax.LineWidth = 1;
set(gca, 'color', 'white');
scatter(stim.stimR(1) + dx/2 , stim.stimR(2) + dx/2, 15, 'r', 'filled');
view(0,90);
xlim([0, topology.L + dx]);
ylim([0, topology.L + dx]);
xticks([]);
yticks([]);
title("Connectivity", 'Interpreter', 'latex', 'FontSize', 16);
hold off;

nexttile;
axis off;
text(1, 1, figtitles(2), 'Interpreter', 'latex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 24, 'Units', 'normalized');

for i = 1:num_snapshots
    ax = nexttile([1 num_subtiles]);
    hold on;
    ax.Box = "on";
    ax.LineWidth = 1;
    imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), ts(:,:,ceil((stim.stimt + timepoints(i)) / dt))');
    view(0,90);
    colormap(ax, Colormap);
    shading flat;
    clim([0 maxval]);
    xlim([0 topology.L + dx]);
    ylim([0 topology.L + dx]);
    xticks([]);
    yticks([]);
    title(coltitles(i), 'Interpreter', 'latex', 'FontSize', 16);
    hold off;
    if i < num_snapshots
        nexttile([1 num_subtiles]);
        hold on;
        axis off;
        text(0.5, 0.5, "$\times \Delta t \ + \ \ldots \ + $", 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 16, 'Units', 'normalized', 'FontWeight','bold');
        hold off;
    else
        nexttile([1 num_subtiles]);
        hold on;
        axis off;
        text(0.5, 0.5, "$\times \Delta t \ + \ \ldots \ = $", 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 16, 'Units', 'normalized', 'FontWeight','bold');
        hold off;
    end
end

ax = nexttile([1 num_subtiles]);
hold on;
ax.Box = "on";
ax.LineWidth = 1;
imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), ts_bold');
set(ax, 'YDir', 'normal');
clim([0 maxvalbold]);
view(0,90);
colormap(ax, Colormap);
shading flat;
xticks([]);
yticks([]);
xlim([0 topology.L + dx]);
ylim([0 topology.L + dx]);
title('$z(\mathbf{r})$', 'Interpreter', 'latex', 'FontSize', 16);
hold off;



% save as schematicbold.svg
print(gcf, 'schematicbold.svg', '-dsvg');
exportgraphics(gcf, 'schematicbold.tiff', 'Resolution', 300);

%% BOLD evoked response of heterogeneous vs homogeneous model

clear; clc;
loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

% Optimal Stimulus
stim.stimR = hetparam_het.a;

% Calculate BOLD response
tol = 1e-5;
bold_full = zeros(topology.Nx, topology.Nx, 2);
bold_full(:, :, 1) = run_bold(topology,homparam,hetparam_hom,stim,tol);
bold_full(:, :, 2) = run_bold(topology,homparam,hetparam_het,stim,tol);

maxvalbold = max(bold_full, [], 'all');

Colormap = [linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];
%
f = figure;
f.Position = [100, 100, 530, 410];

num_subcols = 6;
num_subrows = 10;
tot_subcols = 4 + num_subcols*2;
tot_subrows = 2 + num_subrows*2;

t = tiledlayout(tot_subrows, tot_subcols, "TileSpacing","tight", "Padding", "compact");

modeltitles = "\textbf{" + ["Geometric", "Hybrid"] + "}";
stimulustitles = "\textbf{" + ["i", "ii"] + ".}";
responsetitles = "\textbf{" + ["iii", "iv"] + ".}";
coltitles = strings(1, 2);
coltitles(1) = "Connectivity";
coltitles(2) = append("$z(\mathbf{r})$");

for iter = 1:2
    bold = squeeze(bold_full(:, :, iter));
    % Label Heterogeneous or Homogeneous Model
    nexttile(1 + (iter - 1)*(2+num_subrows)*tot_subcols, [num_subrows, 1]);
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(0.5, 0.5, modeltitles(iter), 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 16, 'Units', 'normalized', 'Rotation', 90);
    % Label Stimulus Subplot
    nexttile(2 + (iter - 1)*(2+num_subrows)*tot_subcols, [num_subrows, 1])
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(1, 1, stimulustitles(iter), 'Interpreter', 'latex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 24, 'Units', 'normalized');
    % Plot Stimulus Subplot
    ax = nexttile(3 + (iter - 1)*(2+num_subrows)*tot_subcols, [num_subrows, num_subcols]);
    hold on;
    ax.Box = "on";
    ax.LineWidth = 1;
    set(gca, 'color', 'white');
    if iter == 2
        quiver(hetparam_het.a(1), hetparam_het.a(2), ...
        hetparam_het.b(1) - hetparam_het.a(1), ...
        hetparam_het.b(2) - hetparam_het.a(2),...
        'Color', 'k', 'LineWidth', 1, ...
        'MaxHeadSize', 0.05 / norm(hetparam_het.a - hetparam_het.b), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end
    scatter(stim.stimR(1), stim.stimR(2), 15, 'r', 'filled');
    text(stim.stimR(1) - 0.02, stim.stimR(2) - 0.02, "$\mathbf{p}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 15);
    if iter == 2
        text(hetparam_het.b(1, 1) + 0.02, hetparam_het.b(2, 1) + 0.02, "$\mathbf{q}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 15);            
    end
    view(0,90);
    xlim([0, topology.L + dx]);
    ylim([0, topology.L + dx]);
    xticks([]);
    yticks([]);
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
    ax = nexttile(5 + num_subcols + (iter - 1)*(2+num_subrows)*tot_subcols, [num_subrows, num_subcols]);
    hold on;
    ax.Box = "on";
    ax.LineWidth = 1;
    imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), bold');
    set(ax, 'YDir', 'normal')
    colormap(ax, Colormap);
    shading flat;
    clim([0 maxvalbold]);
    view(0,90);
    xlim([0, topology.L + dx]);
    ylim([0, topology.L + dx]);
    xticks([]);
    yticks([]);
    if iter == 1
        title(coltitles(2), 'Interpreter', 'latex', 'FontSize', 16);
    end
    hold off;
end

cb = colorbar; 
cb.Layout.Tile = 'east';
cb.FontSize = 15;
cb.TickLabelInterpreter = 'latex';

% save as bolddefaultoptimal.svg
print(gcf, 'bolddefaultoptimal.svg', '-dsvg');
exportgraphics(gcf, 'bolddefaultoptimal.tiff', 'Resolution', 300);


%% Plot default heterogeneous model divergence horizontal line on top of model divergence vs. time curve

clear; clc;
loadparam;

% Optimal stimulus
stim.stimR = hetparam_het.a;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

time_array = dt * (1:topology.Nt);

ts_hom = run_periodic(topology,homparam,hetparam_hom,stim);
ts_het = run_periodic(topology,homparam,hetparam_het,stim);

% Calculate BOLD divergence
tol = 1e-5;
ts_hom_bold = run_bold(topology, homparam, hetparam_hom, stim, tol);
ts_het_bold = run_bold(topology, homparam, hetparam_het, stim, tol);
dissim_bold = pdist2(ts_hom_bold(:)', ts_het_bold(:)', 'cosine');

dissim_time = zeros(1, topology.Nt);
for n = 1:topology.Nt
    dissim_time(n) = pdist2(...
        reshape(ts_hom(:, :, n), [], 1)', reshape(ts_het(:, :, n), [], 1)', 'cosine');
end

% Plot divergence of default model

f = figure;
f.Position = [100, 100, 800, 470];
hold;

plot(1000*(time_array - stim.stimt), dissim_time', 'Color', 'k', 'LineWidth', 2);

ylim([0, 0.09]);
xlim([-2.5, 40]);
xticks([0:10:40]);
yticks([0:0.02:1]);
xlabel("$t \ (\mathrm{ms})$", 'Interpreter', 'latex');
ylabel('$C_\phi$', 'Interpreter', 'latex', 'Rotation', 0)

dissim_max = max(dissim_time);
yline(dissim_max, 'Label', ...
    append('$\max C_{\phi} = ', num2str(dissim_max, "%.3f"),'$'), 'Interpreter', 'latex', 'FontSize', 21, 'Color', 'k', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1, 'LineStyle', '--');
yline(dissim_bold, 'Label', ...
    append('$C_z = ', num2str(dissim_bold, "%.3f"), '$'), 'Interpreter', 'latex', 'FontSize', 21, 'Color', [0.5, 0.5, 0.5], 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1);
% text(1, 0.05, '$C_{\phi}$', 'Interpreter', 'latex', 'FontSize', 28, 'Color', 'b', 'HorizontalAlignment','center','VerticalAlignment','middle');

ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.XAxis.LabelFontSizeMultiplier  = 1.5;
ax.YAxis.LabelFontSizeMultiplier  = 2;

% save as dissimcurveboldoptimal.svg
print(gcf, 'dissimcurveboldoptimal.svg', '-dsvg');
exportgraphics(gcf, 'dissimcurveboldoptimal.tiff', 'Resolution', 300);
