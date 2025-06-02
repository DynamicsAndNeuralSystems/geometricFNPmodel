%% Cosine dissimilarity + BOLD dissimilarity vs stimulus position

clear; clc;

loadparam;

% Create ensemble of 200 stimulus positions
num_samples = 200;
% Set seed
rng(1, "twister");
stim_array = topology.L * rand(2, num_samples);

% tolerance for bold calculation
tol = 1e-5;

% Record total length used and times to thresholds
dissim_time_array = zeros(topology.Nt, num_samples);
dissim_bold_array = zeros(1, num_samples);

parfor k = 1:num_samples
    
    stim1 = stim;
    stim1.stimR = stim_array(:, k);

    % Compute dissimilarity
    ts_hom = run_periodic(topology, homparam, hetparam_hom, stim1);
    ts_het = run_periodic(topology, homparam, hetparam_het, stim1);
    dissim_time = zeros(topology.Nt, 1);
    for n = 1:topology.Nt
        dissim_time(n) = pdist2(...
            reshape(ts_hom(:, :, n), [], 1)', reshape(ts_het(:, :, n), [], 1)', 'cosine');
    end

    % Compute bold dissimilarity
    ts_hom_bold = run_bold(topology, homparam, hetparam_hom, stim1, tol);
    ts_het_bold = run_bold(topology, homparam, hetparam_het, stim1, tol);

    disp(k);

    dissim_time_array(:, k) = dissim_time;
    dissim_bold_array(k) = pdist2(...
            ts_hom_bold(:)', ts_het_bold(:)', 'cosine');
    
end
%
save('dissim_stimposition.mat', "topology", "num_samples", "stim_array", "dissim_time_array", "dissim_bold_array");

%% F. Plot dissimilarity of default model over time, optimal and nonoptimal stimulus

clear; clc;
loadparam;

% Plot dissimilarity with optimal Position
stim.stimR = hetparam_het.a(:, 1);
ts_hom = run_periodic(topology,homparam,hetparam_hom,stim);
ts_het = run_periodic(topology,homparam,hetparam_het,stim);
dissim_time_optimal = zeros(1, topology.Nt);
for n = 1:topology.Nt
    dissim_time_optimal(n) = pdist2(...
        reshape(ts_hom(:, :, n), [], 1)', reshape(ts_het(:, :, n), [], 1)', 'cosine');
end


dt = topology.T / topology.Nt;

% Create figure
fig = figure;
fig.Position = [100, 100, 1400, 500];
hold;

% Non-optimal Stimulation
load('dissim_stimposition.mat');
for i = 1:num_samples
    plot(1000*(dt * (1:topology.Nt) - stim.stimt), dissim_time_array(:, i), 'Color', 0.3*[0 0 1] + 0.7*[1 1 1]);
end

% Optimal Stimulation
plot(1000*(dt * (1:topology.Nt) - stim.stimt), dissim_time_optimal', 'Color', 'k', 'LineWidth', 1);

% Aggregate of all non-optimal stimulation
mean_dissim_time = squeeze(mean(dissim_time_array, 2));
plot(1000*(dt * (1:topology.Nt) - stim.stimt), mean_dissim_time, 'Color', 'b', 'LineWidth', 1);

Legend = cell(1, 2 + num_samples);
for i = 1:num_samples
    Legend{i} = "";
end
Legend{num_samples + 1} = "Stimulus positioned at FNP source ($\mathbf{r}_0 = \mathbf{p}$)";
Legend{num_samples + 2} = append("Average over $", num2str(num_samples), "$ uniformly sampled", newline, "stimulus positions ($\mathbf{r}_0 \sim \mathcal{U}(\Omega)$)");

xlim([-5, Inf]);
ylim([-0.0005, 0.09]);
xticks([0:10:40]);
yticks([0:0.02:0.1]);
xlabel("$t \ (\mathrm{ms})$", 'Interpreter', 'latex');
ylabel('$C_{\phi}$', 'Interpreter', 'latex', 'Rotation', 0);
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.XAxis.LabelFontSizeMultiplier  = 1.5;
ax.YAxis.LabelFontSizeMultiplier  = 2;

l = legend(Legend, 'Interpreter', 'latex', 'Location', 'northeast', 'Box', 'off', 'FontSize', 16);

% save as SI_dissimcurve_stimposition.tiff
exportgraphics(gcf, 'S1_dissimcurve_stimposition.tiff', 'Resolution', 300);
close(fig);

%% Plot (max dissim + bold dissim) vs (distance of stim position from p)

clear; clc;
loadparam;

% Compute max dissim + dissim bold with optimal Position
stim.stimR = hetparam_het.a(:, 1);
ts_hom = run_periodic(topology,homparam,hetparam_hom,stim);
ts_het = run_periodic(topology,homparam,hetparam_het,stim);
dissim_time_optimal = zeros(1, topology.Nt);
for n = 1:topology.Nt
    dissim_time_optimal(n) = pdist2(...
        reshape(ts_hom(:, :, n), [], 1)', reshape(ts_het(:, :, n), [], 1)', 'cosine');
end
max_dissim_opt = max(dissim_time_optimal);

tol = 1e-5;
ts_hom_bold = run_bold(topology, homparam, hetparam_hom, stim, tol);
ts_het_bold = run_bold(topology, homparam, hetparam_het, stim, tol);
dissim_bold_opt = pdist2(...
            ts_hom_bold(:)', ts_het_bold(:)', 'cosine');


load('dissim_stimposition.mat');
% Calculate dist of sampled stim positions from p
dist_array = zeros(1, num_samples);
for i = 1:num_samples
    distx = abs(stim_array(1, i) - hetparam_het.a(1));
    distx = min(distx, topology.L - distx);
    disty = abs(stim_array(2, i) - hetparam_het.a(2));
    disty = min(disty, topology.L - disty);
    dist_array(i) = norm([distx; disty], 2);
end

% Compute max + time to max
max_dissim_array = max(dissim_time_array, [], 1);


dt = topology.T / topology.Nt;

% Create figure and properties
fig = figure;
fig.Position = [100 100 1400 350];
Titles = "\textbf{" + ["i", "ii"] + ".}";

for i = 1:2
    subplot(1, 2, i); hold;
    if i == 1
        scatter([0 dist_array], [max_dissim_opt max_dissim_array], 'k', 'filled');
    else
        scatter([0 dist_array], [dissim_bold_opt dissim_bold_array], 'k', 'filled');
    end
    xlabel("$\|\mathbf{r}_0 - \mathbf{p}\| \ (\mathrm{m})$", 'Interpreter', 'latex');
    if i == 1
        ylabel('$\max_t \{ C_{\phi}(t) \}$', 'Interpreter', 'latex', 'Rotation', 0);
    else
        ylabel('$ C_z $', 'Interpreter', 'latex', 'Rotation', 0);
    end
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.XAxis.FontSize = 14;
    ax.YAxis.FontSize = 14;
    ax.XAxis.LabelFontSizeMultiplier  = 1.5;
    ax.YAxis.LabelFontSizeMultiplier  = 1.5;
    title(Titles(i), 'Interpreter', 'latex', 'FontSize', 24);
    hold off;
end

% save as 5_dissim_stimposition.tiff
exportgraphics(gcf, 'S1_dissim_stimposition.tiff', 'Resolution', 300);
close(fig);

%% Plot SCF with geometric connectivity vs input width
clear; clc;
loadparam;

dx = topology.L / topology.Nx;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

load("scfvssigma_n_geometric.mat")

% Specify coordinates of p
position_centre = hetparam_het.a;
position_centre_grid = position_centre / dx;

% Remove corr(p, p) = 1
for count = 1:length(sigma_n_array)
    corr_array(position_centre_grid(1), position_centre_grid(2), count) = Inf;
end
    
% Create tiledlayout figure and properties
fig = figure;
fig.Position = [100, 50, 1170, 465];
num_subrows = 4;
tot_subrows = 2 + num_subrows;
t = tiledlayout(tot_subrows, num_experiments, 'TileSpacing', 'compact', 'Padding', 'tight');
titles = "\textbf{" + ["i", "ii", "iii", "iv"] + ".}";
sigma_n_titles = "$\sigma_n=" + ["3L/80", "L/8", "L/4", "\infty"] + "$";

% Title SCF heatmaps
nexttile(1, [1, num_experiments])
axis off;
text(0, 1, "$\gamma_\mathbf{p}^\mathrm{geo}(\mathbf{r})$ over $\Omega$ - varying $\sigma_n$", 'Interpreter', 'latex', 'FontSize', 24, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')

% Plot SCF heatmaps
for count = 1:num_experiments
    sigma_n = sigma_n_array(count);
    ax = nexttile(num_experiments + count, [num_subrows, 1]);
    box on;
    ax.LineWidth = 1;
    ax.Color = 'k';
    hold on;
    imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), corr_array(:, :, count)');
    set(ax, 'YDir', 'normal');    
    colormap(flipud(pink))
    view(0, 90);
    shading flat;
    % Draw circle/square
    if sigma_n < Inf
        theta = linspace(0,2*pi,1000);
        x = sigma_n * cos(theta) + position_centre(1);
        y = sigma_n * sin(theta) + position_centre(2);
        fill(x, y, 'r', 'FaceAlpha', 0., 'EdgeColor', 'r', 'LineWidth', 1);
    else
        fill([dx dx topology.L topology.L], [dx topology.L topology.L dx], 'r', 'FaceAlpha', 0., 'EdgeColor', 'r', 'LineWidth', 1)
    end
    xlim([0 topology.L + dx])
    ylim([0 topology.L + dx])
    clim([0 Inf]);
    xticks([]);
    yticks([]);
    ax.TickLength = [0 0];
    title(titles(count), 'Interpreter', 'latex', 'FontSize', 24, 'Units', 'normalized');
    ax.TitleHorizontalAlignment = 'left';
    if count == 1
        text(hetparam_het.a(1) - 0.005, hetparam_het.a(2) - 0.005, '$\mathbf{p}$', 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        text(topology.L - 0.01, 0, "$\Omega$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontSize', 16);
    end
    text(0.02, 0.99, sigma_n_titles(count), 'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 20);
    cb = colorbar;
    cb.FontSize = 10;
    cb.TickLabelInterpreter = 'latex';
    cb.Location = 'southoutside';
    title(cb, "$\gamma_{\mathbf{p}}^{\mathrm{geo}}(\mathbf{r})$", 'Interpreter', 'latex', 'Position', [100 -40 0], 'FontSize', 15);
    hold off;   
end

% save as scfvssigma_n_geometric.tiff
exportgraphics(gcf, 'S2_scfvssigma_n_geometric.tiff', 'Resolution', 300);
close(fig);

%% Generate 1d schematic for derivation

clear; clc;
num_nodes = 121;
theta_array = linspace(0, 2*pi, num_nodes);
fig = figure;
hold on;
fig.Position = [200 200 1350 225];
max_line_length = 10;
for j = max_line_length:-1:1
    color = ((max_line_length - j) * [0 0 1] + j * [1 1 1]) / max_line_length;
    for i = 1:num_nodes
        k = i + j;
        if k > num_nodes
            k = k - num_nodes;
        end
        line([cos(theta_array(i)), cos(theta_array(k))], [sin(theta_array(i)), sin(theta_array(k))], 'Color', color);
    end
end
for i = 1:num_nodes
    line([cos(theta_array(i)), 1.05*cos(theta_array(i))], [sin(theta_array(i)), 1.05*sin(theta_array(i))], 'Color', 'k');
end
plot(cos(theta_array), sin(theta_array), 'Marker', 'o', 'MarkerSize', 4, 'Color', 'blue', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue', 'LineWidth', 2);
plot(1.05*cos(theta_array), 1.05*sin(theta_array), 'Marker', 'o', 'MarkerSize', 2, 'Color', 'red', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red', 'LineWidth', 1);

% plot(0.95*cos(2*pi/3), 0.95*sin(2*pi/3), 'Marker', 'o', 'Color', [0 0.5 0], 'MarkerSize', 3, 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0]);
% plot(0.9*cos(pi/3), 0.9*sin(pi/3), 'Marker', 'o', 'Color', [0 0.5 0], 'MarkerSize', 3, 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0]);

%scatter(1.05*cos(theta_array), 1.05*sin(theta_array), 10, 'red', 'o', 'filled');
xlim([-1 1]); ylim([0.7 Inf]);
xticks([]); yticks([]);

% export as physiologicalschematic.tiff
exportgraphics(gcf, 'S3_physiologicalschematic.tiff', 'Resolution', 300);
close(fig);