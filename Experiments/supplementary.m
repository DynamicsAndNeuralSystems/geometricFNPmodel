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

%% Plot dissimilarity of default model over time, optimal and nonoptimal stimulus

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
ylabel('$C(t)$', 'Interpreter', 'latex', 'Rotation', 0);
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
        ylabel('$C_{\max}$', 'Interpreter', 'latex', 'Rotation', 0);
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

%% Compute correlation structure of geometric model

clear; clc;

loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

tol = 1e-7;

num_parcel_x = 10;

% Set seed for generating noisy input
rng(0, "twister");
corr_array = compute_spont_corr(topology, dx, dt, hetparam_hom, homparam, num_parcel_x, tol);

save('corr_parcel_geometric.mat', 'num_parcel_x', 'corr_array', 'topology')

%% Compute correlation structure of random connectome
% N = 100

clear; clc;

loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

tol = 1e-7;

load('corr_parcel_geometric.mat', 'num_parcel_x');

% Set number of FNPs
N_array = [100];
max_N = max(N_array);

% Create Random Connectome
rng(0, "twister");
a_array = zeros(2, max_N);
b_array = zeros(2, max_N);
for j = 1:max_N
    a = topology.L * rand(2, 1);
    b = topology.L * rand(2, 1);
    a_array(:, j) = a;
    b_array(:, j) = b;
end

% Compute BOLD dissimilarity

N = N_array(1);

hetparam_het.m = N;
hetparam_het.c = (homparam.r)^2 * ones(1, N);
hetparam_het.tau = zeros(1, N);
hetparam_het.a = a_array(:, 1:N);
hetparam_het.b = b_array(:, 1:N);

corr_array = compute_spont_corr(topology, dx, dt, hetparam_het, homparam, num_parcel_x, tol);


save('corr_parcel_hybrid.mat', 'num_parcel_x', 'hetparam_het', 'corr_array', 'topology')

%% Plot SCF of geometric and sample hybrid model with random connectome
% At parcel (5, 5)

clear; clc;
loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

num_models = 2;

% Set reference point
ref_pt = [5, 5];

load('corr_parcel_geometric.mat');
corr_array_hom = corr_array;
load('corr_parcel_hybrid.mat', 'corr_array', 'hetparam_het');
corr_array = cat(5, corr_array_hom, corr_array);
for k = 1:num_models
    for i = 1:num_parcel_x
        for j = 1:num_parcel_x
            for i1 = 1:num_parcel_x
                for j1 = 1:num_parcel_x
                    if isnan(corr_array(i, j, i1, j1, k)) && ~isnan(corr_array(i1, j1, i, j, k))
                        corr_array(i, j, i1, j1, k) = corr_array(i1, j1, i, j, k);
                    elseif ~isnan(corr_array(i, j, i1, j1, k)) && isnan(corr_array(i1, j1, i, j, k))
                        corr_array(i1, j1, i, j, k) = corr_array(i, j, i1, j1, k);
                    elseif isnan(corr_array(i, j, i1, j1, k)) && isnan(corr_array(i1, j1, i, j, k))
                        corr_array(i, j, i1, j1, k) = Inf;
                    end
                end
            end
        end
    end
end
corr_array = squeeze(corr_array(ref_pt(1), ref_pt(2), :, :, :));

hetparam_het_array = {hetparam_hom, hetparam_het};

% Create figure and properties
fig = figure;
fig.Position = [100 100 520 450];
num_subcols = 8;
num_subrows = 3;
tot_subcols = 2 + num_models*num_subcols;
tot_subrows = 1 + 2*num_subrows;
t = tiledlayout(tot_subrows, tot_subcols, 'TileSpacing', 'tight', 'Padding', 'compact');
modeltitles = ["Geometric", "Hybrid"];
connectivitytitles = "\textbf{" + ["i", "ii"] + ".}";
dissimtitles = "\textbf{" + ["iii", "iv"] + ".}";

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
text(0.5, 0.5, 'Connectivity', 'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline', 'Rotation', 90, 'Units', 'normalized')

% Draw connectivities
for iter = 1:num_models
    ax = nexttile(tot_subcols + 3 + (iter - 1)*(num_subcols), [num_subrows, num_subcols]);
    hold on;
    box on;
    ax.LineWidth = 1;
    ax.Color = 'k';
    set(gca, 'Color', 'white');
    hetparam_het = hetparam_het_array{iter};
    for m = 1:hetparam_het.m
        q = quiver(hetparam_het.a(1, m), hetparam_het.a(2, m), ...
        hetparam_het.b(1, m) - hetparam_het.a(1, m), ...
        hetparam_het.b(2, m) - hetparam_het.a(2, m),...
        'Color', 'k', 'LineWidth', 0.1, ...
        'MaxHeadSize', 0.05 / norm(hetparam_het.a(:, m) - hetparam_het.b(:, m)), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end
    if iter == 1
        text(topology.L - 0.01, 0, "$\Omega$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontSize', 15);
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
text(0.5, 0.5, ['$\gamma$ with reference', newline, 'parcel $\Omega_{', num2str(ref_pt(1)),',', num2str(ref_pt(2)),'}$'], 'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline', 'Rotation', 90, 'Units', 'normalized')

% Plot heatmaps of max dissimilarities
for iter = 1:num_models
    ax = nexttile((num_subrows + 1)*tot_subcols + 3 + (iter - 1)*(num_subcols), [num_subrows, num_subcols]);
    hold on;
    box on;
    ax.LineWidth = 1;
    ax.Color = 'k';
    set(gca, 'Color', 'white');
    hetparam_het = hetparam_het_array{iter};
    imagesc((topology.L / num_parcel_x)*(-0.5 + 1:num_parcel_x),(topology.L / num_parcel_x)*(-0.5 + 1:num_parcel_x),squeeze(corr_array(:, :, iter))); 
    clim([0  1]); 
    shading flat; view(0, 90); 
    colormap(flipud(pink));
    set(gca, 'YDir', 'normal');
    ylabel([dissimtitles(iter), '\\', '\\', '\\'], 'Interpreter', 'latex', 'FontSize', 24, 'Rotation', 0);
    xlim([-dx, topology.L + dx]);
    ylim([-dx, topology.L + dx]);
    xticks([]);
    yticks([]);
    hold off;
end

cb = colorbar; 
cb.FontSize = 15;
cb.TickLabelInterpreter = 'latex';
cb.Ticks = [0, 1];
cb.TickLabels = {'0', '1'};

% save as scf_random.svg
exportgraphics(gcf, 'S3_scf_random.tiff', 'Resolution', 300);
close(fig);

%% Plot SCF of geometric and sample hybrid model with random connectome
% At parcel (5, 5)

clear; clc;
loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

num_models = 2;

% Set reference point
ref_pt = [5, 5];

load('corr_parcel_geometric.mat');
corr_array_hom = corr_array;
load('corr_parcel_hybrid.mat', 'corr_array', 'hetparam_het');
corr_array = cat(5, corr_array_hom, corr_array);


% Create figure and properties
fig = figure;
hold;
fig.Position = [100 100 400 450];
ax = gca;

vec1 = reshape(corr_array(:, :, :, :, 1), [], 1);
vec2 = reshape(corr_array(:, :, :, :, 2), [], 1);
vec1 = vec1(~isnan(vec1)); vec2 = vec2(~isnan(vec2));
xlim([floor(min(vec1, [], 'all')/0.1) ceil(max(vec1, [], 'all')/0.1)] * 0.1);
ylim([floor(min(vec2, [], 'all')/0.1) ceil(max(vec2, [], 'all')/0.1)] * 0.1);
% Compute correlation between SCF's and annotate
corr = 1 - pdist2(vec1', vec2', 'correlation');
text(max(xlim), min(ylim), append("r = $",num2str(corr, 2), "$"), 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
% Compute linear fit and annotate
coefficients = polyfit(vec1, vec2, 1);
plot(xlim, polyval(coefficients, xlim), 'Color', [0.6, 0.6, 0.6], 'LineStyle', '--', 'LineWidth', 1);
% Plot SCF points on top of linear fit
scatter(vec1, vec2, 2, 'black', 'filled');
ax.TickLabelInterpreter = 'latex';
xlabel('$\gamma^\mathrm{geo}$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$\gamma^\mathrm{hyb}$', 'Interpreter', 'latex', 'FontSize', 15);
ax.XAxis.FontSize = 10;
ax.YAxis.FontSize = 10;
ax.LabelFontSizeMultiplier  = 1.5;
ax.TitleHorizontalAlignment = 'left';
hold off;

% save as corr_scatter_random.svg
exportgraphics(gcf, 'S3_corr_scatter_random.tiff', 'Resolution', 300);
close(fig);


%% Compute corr dissimilarity of random connectome
% N = 100

clear; clc;

loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

tol = 1e-5;

load('corr_parcel_geometric.mat', 'num_parcel_x', 'corr_array');
corr_array_hom = corr_array; clear corr_array;

% Set number of FNPs
N_array = [100];
max_N = max(N_array);

% Set number samples - Change for 100 for figures in paper
num_samples = 20;

% Record BOLD dissimilarity
corr_dissim_array = zeros(length(N_array), num_samples);

parfor k = 1:num_samples

    % Create duplicates for parallelization
    hetparam_het1 = hetparam_het;

    corr_dissim_subarray = zeros(length(N_array), 1);

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

        corr_array = compute_spont_corr(topology, dx, dt, hetparam_het1, homparam, num_parcel_x, tol);

        vec1 = corr_array_hom(:);
        vec1 = vec1(~isnan(vec1));
        vec2 = corr_array(:);
        vec2 = vec2(~isnan(vec2));
        dissim_corr = pdist2(vec1', vec2', 'correlation');
        vec1 = 0; vec2 = 0;

        corr_dissim_subarray(i) = dissim_corr;
        disp([i, k, dissim_corr]);

    end

    corr_dissim_array(:, k) = corr_dissim_subarray;

end

save('dissim_corr_multilrcs.mat', "topology", "N_array", "num_samples", "corr_dissim_array");


%% Compute corr dissimilarity versus level of hub specificity
% N = 100

clear; clc;
 
loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

tol = 1e-5;

load('corr_parcel_geometric.mat', 'num_parcel_x', 'corr_array');
corr_array_hom = corr_array; clear corr_array;

% Set number of FNPs
N_array = [100];
max_N = max(N_array);

% Set number samples - Change to 100 for figures in paper
num_samples = 20;

% Set number of hubs, positions, and dimensions
num_hubs = 4;
hub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
hublength = topology.L / sqrt(34);

% Set lambda (level of hub specificity)
lambda_array = 0.2:0.2:1;

corr_dissim_array = zeros(length(lambda_array), length(N_array), num_samples);

parfor k = 1:num_samples

    % Set seed for worker
    rng(k, "twister");

    % Create duplicates for parallelization
    hetparam_het1 = hetparam_het;

    corr_dissim_subarray = zeros(length(lambda_array), length(N_array));

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
    
            corr_array = compute_spont_corr(topology, dx, dt, hetparam_het1, homparam, num_parcel_x, tol);
    
            vec1 = corr_array_hom(:);
            vec1 = vec1(~isnan(vec1));
            vec2 = corr_array(:);
            vec2 = vec2(~isnan(vec2));
            dissim_corr = pdist2(vec1', vec2', 'correlation');
            vec1 = 0; vec2 = 0;
        
            corr_dissim_subarray(i, j) = dissim_corr;
            disp(num2str([i, j, k, dissim_corr]));
    
        end
    
        corr_dissim_array(:, :, k) = corr_dissim_subarray;
    
    end

end
%
save('dissim_corr_hub.mat', "topology", "N_array", "num_hubs", "hub_centres", "hublength", "lambda_array", "num_samples", "corr_dissim_array");
% 

%% Compute corr dissimilarity versus level of rich-club specificity
% N = 100

clear; clc;
 
loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

tol = 1e-5;

load('corr_parcel_geometric.mat', 'num_parcel_x', 'corr_array');
corr_array_hom = corr_array; clear corr_array;

% Set number of FNPs
N_array = [100];
max_N = max(N_array);

% Set number samples - Change to 100 for figures in paper
num_samples = 20;

% Set number of hubs, positions, and dimensions
num_hubs = 4;
hub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
hublength = topology.L / sqrt(34);

% Set lambda (level of hub specificity)
lambda_array = 0.2:0.2:1;

corr_dissim_array = zeros(length(lambda_array), length(N_array), num_samples);

parfor k = 1:num_samples

    % Set seed for worker
    rng(k, "twister");

    % Create duplicates for parallelization
    hetparam_het1 = hetparam_het;

    corr_dissim_subarray = zeros(length(lambda_array), length(N_array));

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
    
            corr_array = compute_spont_corr(topology, dx, dt, hetparam_het1, homparam, num_parcel_x, tol);
    
            vec1 = corr_array_hom(:);
            vec1 = vec1(~isnan(vec1));
            vec2 = corr_array(:);
            vec2 = vec2(~isnan(vec2));
            dissim_corr = pdist2(vec1', vec2', 'correlation');
            vec1 = 0; vec2 = 0;
        
            corr_dissim_subarray(i, j) = dissim_corr;
            disp(num2str([i, j, k, dissim_corr]));
    
        end
    
        corr_dissim_array(:, :, k) = corr_dissim_subarray;
    
    end

end
%
save('dissim_corr_core.mat', "topology", "N_array", "num_hubs", "hub_centres", "hublength", "lambda_array", "num_samples", "corr_dissim_array");
% 


%% Plot corr dissimilarity versus hub and rich-club specificity

clear; clc;
loadparam;

fig = figure;
fig.Position = [100, 100, 500, 450];
hold;
colors = [1 0 0; 0 0.5 0];
ax = gca;


% Load random connectome statistics
load('dissim_corr_multilrcs.mat', 'corr_dissim_array')
corr_dissim_array0 = squeeze(corr_dissim_array);
clear corr_dissim_array;

load('dissim_corr_hub.mat', "N_array", 'lambda_array', 'corr_dissim_array');


% Append random to nonrandom connectome statistics [random, hub]
corr_dissim_array = cat(1, corr_dissim_array0, squeeze(corr_dissim_array));
lambda_array = [0 lambda_array];

mean_dissimilarity = squeeze(mean(corr_dissim_array, 2));
std_dissimilarity = squeeze(std(corr_dissim_array, 1, 2));

for i = 1:length(N_array)
    h = errorbar(lambda_array, mean_dissimilarity(:, i), std_dissimilarity(:, i), 'o', 'MarkerSize', 0.000001, 'LineWidth', 1, 'Color', colors(1, :));
    alpha = 0.2;   
    % Set transparency (undocumented)
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
end
for i = 1:length(N_array)
    plot(lambda_array, mean_dissimilarity(:, i), 'Marker', 'o', 'Color', colors(1, :), 'MarkerFaceColor', colors(1, :));
end

load('dissim_corr_core.mat', 'N_array', 'lambda_array', 'corr_dissim_array');

% Append random to nonrandom connectome statistics [random, rich-club]
corr_dissim_array = cat(1, corr_dissim_array0, squeeze(corr_dissim_array));
lambda_array = [0 lambda_array];

mean_dissimilarity = squeeze(mean(corr_dissim_array, 2));
std_dissimilarity = squeeze(std(corr_dissim_array, 1, 2));

for i = 1:length(N_array)
    h = errorbar(lambda_array, mean_dissimilarity(:, i), std_dissimilarity(:, i), 'o', 'MarkerSize', 0.000001, 'LineWidth', 1, 'Color', colors(2, :));
    alpha = 0.2;   
    % Set transparency (undocumented)
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
end
for i = 1:length(N_array)
    plot(lambda_array, mean_dissimilarity(:, i), 'Marker', 'o', 'Color', colors(2, :), 'MarkerFaceColor', colors(2, :));
end

xlabel('$\lambda$', 'Interpreter', 'latex');
ylabel('$1 - \mathrm{r}$', 'Interpreter', 'latex', 'Rotation', 0);
ax.TickLabelInterpreter = 'latex';
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.XAxis.LabelFontSizeMultiplier  = 1.5;

xlim([-0.05*max(xlim) 1.05*max(xlim)])
yticks([0:0.1:0.6]); ylim([0.15 0.55]);
xticks(lambda_array);

hold off;

Legend = {"", "Hub Specificity ($\lambda_h$)", "", "Rich-club Specificity ($\lambda_r$)"};

l = legend(Legend, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 16, 'Location', 'southeast');
l.Orientation = 'vertical';

% save as dissim_max_nonrandom.tiff
exportgraphics(gcf, 'S3_dissim_corr_nonrandom.tiff', 'Resolution', 300);
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

xlim([-1 1]); ylim([0.7 Inf]);
xticks([]); yticks([]);

% export as physiologicalschematic.tiff
exportgraphics(gcf, 'S4_physiologicalschematic.tiff', 'Resolution', 300);
close(fig);