%% Create schematic of correlation plot

clear; clc; 
% Load model parameters
loadparam; 

% Compute grid spacing
dx = topology.L / topology.Nx;

fig = figure;
fig.Position = [100, 100, 520, 500];
ax = gca;
set(ax, 'color', 'white');
axis off;
hold on;
xticks([]);
yticks([]);
xlim([0 topology.L]);
ylim([0 topology.L]);

% Plot targeted noisy input with specified characteristic width (sigma_n)
xposition = hetparam_het.a(1);
yposition = hetparam_het.a(2);
sigma_n = 0.04;
stimulationinputspace = zeros(topology.Nx, topology.Nx);
i0 = xposition / dx; j0 = yposition / dx;
for i = 1:topology.Nx
    for j = 1:topology.Nx
        distx = abs(i - i0);
        if distx > topology.Nx / 2
            distx = topology.Nx - distx;
        end
        disty = abs(j - j0);
        if disty > topology.Nx / 2
            disty = topology.Nx - disty;
        end
        stimulationinputspace(i, j) = exp(-0.5*(distx^2 + disty^2) * (dx^2) / (sigma_n^2));
    end
end
% Set seed for generating noisy input
rng(0, "twister");
% Compute noisy input
stimulationinputspace = stimulationinputspace.* randn(topology.Nx);
imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), stimulationinputspace');
set(gca, 'YDir', 'normal');
colormap([linspace(0, 1, 256)' linspace(0, 1, 256)', linspace(1, 1, 256)';
    linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)']);
clim([-max(stimulationinputspace, [], 'all') max(stimulationinputspace, [], 'all')])

% Draw circle with sigma radius
scatter(xposition , yposition, 20, 'r', 'filled');
theta = linspace(0,2*pi,1000);
x = sigma_n * cos(theta) + xposition;
y = sigma_n * sin(theta) + yposition;
fill(x, y, 'r', 'FaceAlpha', 0., 'EdgeColor', 'r', 'LineWidth', 1);
quiver(xposition - sigma_n, yposition - 2*sigma_n, 2*sigma_n, 0, ...
    'Color', 'k', 'LineWidth', 1, ...
    'MaxHeadSize', 0.3 / norm(hetparam_het.a - hetparam_het.b), ...
    'Marker', '.', 'MarkerSize', 0.0001, ...
    'AutoScale','off');
quiver(xposition + sigma_n, yposition - 2*sigma_n, -2*sigma_n, 0, ...
'Color', 'k', 'LineWidth', 1, ...
'MaxHeadSize', 0.3 / norm(hetparam_het.a - hetparam_het.b), ...
'Marker', '.', 'MarkerSize', 0.0001, ...
'AutoScale','off');
text(xposition, yposition - 2*sigma_n - 0.01, "$2\sigma_n$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','top', 'FontSize', 14);

% Draw circle on p, the source of the FNP
scatter(hetparam_het.a(1) , hetparam_het.a(2), 40, 'k', 'filled');
text(hetparam_het.a(1), hetparam_het.a(2) - 0.015, "$\mathbf{p}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 18);
text(xposition - 0.05, yposition + 0.02, "$\mathbf{i}\textbf{.}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 24);


% Draw circle on r, a random point on Omega
randpt = [0.25, 0.3];
scatter(randpt(1), randpt(2), 20, 'k', 'filled');
text((hetparam_het.a(1) + randpt(1))/2 + 0.01, (hetparam_het.a(2) + randpt(2))/2 - 0.01, ["$\gamma_{\mathbf{p}}(\mathbf{r}) = \mathrm{corr}($","$\hspace{1em} \phi(\mathbf{r}, t), \phi(\mathbf{p}, t))$"], 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','left','VerticalAlignment','middle', 'FontSize', 15);
text(randpt(1), randpt(2) + 0.015, "$\mathbf{r}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 18);

% Annotate correlation between FNP source and random point
plot([hetparam_het.a(1), randpt(1)], [hetparam_het.a(2), randpt(2)], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
text((hetparam_het.a(1) + randpt(1))/2, (hetparam_het.a(2) + randpt(2))/2 - 0.01, "$\mathbf{ii}\textbf{.}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontSize', 24);
text(topology.L - 0.01, 0, "$\Omega$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontSize', 24);

% Draw border
pos = get(ax, 'Position');
annotation('rectangle', [pos(1), pos(2), pos(1) + pos(3) - pos(1), pos(2) + pos(4) - pos(2)], ...
    'LineWidth', 1, 'Color', 'k');

% save as schematiccorrelation.tiff
exportgraphics(gcf, '5_schematiccorrelation.tiff', 'Resolution', 300);
close(fig);

%% Calculate SCF of p with hybrid connectivity vs spatial width of noisy input

clear; clc;
% Load model parameters
loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

% Specify p, the source of the FNP
position_centre = hetparam_het.a;
position_centre_grid = position_centre / dx;

% Specify selected values of sigma to investigate
sigma_n_array = [0.015, 0.05, 0.10, Inf];
num_experiments = length(sigma_n_array);
corr_array = zeros(topology.Nx, topology.Nx, length(sigma_n_array));

% Compute SCF with reference point p
parfor count = 1:length(sigma_n_array)
    % Set seed for generating noisy input
    rng(0, "twister");
    sigma_n = sigma_n_array(count);
    tol = 1e-4; % change to 1e-7 for tolerance used in figures
    corr_array(:, :, count) = compute_scf(sigma_n, position_centre_grid, topology, dx, dt, hetparam_het, homparam, tol);
end

save('scfvssigma_n.mat', 'num_experiments', 'sigma_n_array', 'corr_array', 'topology')

%% Calculate SCF of p with geometric connectivity vs spatial width of noisy input

clear; clc;
% Load model parameters
loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

% Specify p, the source of the FNP
position_centre = hetparam_het.a;
position_centre_grid = position_centre / dx;

% Specify selected values of sigma to investigate
sigma_n_array = [0.015, 0.05, 0.10, Inf];
num_experiments = length(sigma_n_array);
corr_array = zeros(topology.Nx, topology.Nx, length(sigma_n_array));

% Compute SCF with reference point p
parfor count = 1:length(sigma_n_array)
    % Set seed for generating noisy input
    rng(0, "twister");
    sigma_n = sigma_n_array(count);
    tol = 1e-4; % change to 1e-7 for tolerance used in figures
    corr_array(:, :, count) = compute_scf(sigma_n, position_centre_grid, topology, dx, dt, hetparam_hom, homparam, tol);
end

save('scfvssigma_n_geometric.mat', 'num_experiments', 'sigma_n_array', 'corr_array', 'topology')

%% Plot SCF of p with hybrid connectivity vs input width

clear; clc;
% Load model parameters
loadparam;

% Compute grid spacing
dx = topology.L / topology.Nx;

% Load SCF under geometric and hybrid connectivity
load('scfvssigma_n.mat')
corr_array_het = corr_array;
load("scfvssigma_n_geometric.mat")
corr_array_hom = corr_array;
clear corr_array;

% Specify coordinates of p, the source of the FNP
position_centre = hetparam_het.a;
position_centre_grid = position_centre / dx;
% Specify coordinates of q, the target of the FNP
position_b_grid = hetparam_het.b / dx;

% Remove corr(p, p) = 1
for count = 1:length(sigma_n_array)
    corr_array_het(position_centre_grid(1), position_centre_grid(2), count) = Inf;
    corr_array_hom(position_centre_grid(1), position_centre_grid(2), count) = Inf;
end

% Create tiledlayout figure and properties
fig = figure;
fig.Position = [100, 50, 1170, 910];
num_subcols = 8;
num_subrows = 4;
tot_subcols = num_subcols*num_experiments;
tot_subrows = 2 + 2*num_subrows;
t = tiledlayout(tot_subrows, tot_subcols, 'TileSpacing', 'loose', 'Padding', 'tight');
titles = "\textbf{" + ["i", "ii", "iii", "iv"] + ".}";
sigma_n_titles = "$\sigma_n=" + ["3L/80", "L/8", "L/4", "\infty"] + "$";

% Label SCF heatmaps
nexttile(1, [1, tot_subcols])
axis off;
text(0, 1, "$\gamma_\mathbf{p}^\mathrm{hyb}(\mathbf{r})$ over $\Omega$ - varying $\sigma_n$", 'Interpreter', 'latex', 'FontSize', 24, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')

% Plot SCF heatmaps
for count = 1:num_experiments
    sigma_n = sigma_n_array(count);
    ax = nexttile(1 + tot_subcols + (count - 1)*num_subcols, [num_subrows num_subcols]);
    box on;
    ax.LineWidth = 1;
    ax.Color = 'k';
    hold on;
    imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), corr_array_het(:, :, count)');
    set(ax, 'YDir', 'normal');
    colormap(flipud(pink))
    view(0, 90);
    shading flat;
    % Draw circle with center p radius sigma_n, or Omega if sigma_n = Inf
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
    % Annotate p and q (target of FNP) in first plot
    if count == 1
        text(hetparam_het.a(1) - 0.005, hetparam_het.a(2) - 0.005, '$\mathbf{p}$', 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        text(hetparam_het.b(1) + 0.005, hetparam_het.b(2) + 0.005, '$\mathbf{q}$', 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        text(topology.L - 0.01, 0, "$\Omega$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontSize', 16);
    end
    text(0.02, 0.99, sigma_n_titles(count), 'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 20);
    cb = colorbar;
    cb.FontSize = 10;
    cb.TickLabelInterpreter = 'latex';
    cb.Location = 'southoutside';
    title(cb, "$\gamma_{\mathbf{p}}^{\mathrm{hyb}}(\mathbf{r})$", 'Interpreter', 'latex', 'Position', [105 -40 0], 'FontSize', 15);
    hold off;   
end

% Label SCF/SCF scatter plot
nexttile(1 + tot_subcols * (1 + num_subrows), [1, tot_subcols])
axis off;
text(0, 1, "$(\gamma_\mathbf{p}^\mathrm{geo}(\mathbf{r}), \gamma_\mathbf{p}^\mathrm{hyb}(\mathbf{r}))$ - varying $\sigma_n$", 'Interpreter', 'latex', 'FontSize', 24, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')



for count = 1:num_experiments
    ax = nexttile(1 + (2 + num_subrows) * tot_subcols + (count - 1)*num_subcols, [num_subrows num_subcols]);
    hold;
    vec1 = reshape(corr_array_hom(:, :, count), [], 1);
    vec2 = reshape(corr_array_het(:, :, count), [], 1);
    vec1 = vec1(vec1 < Inf); vec2 = vec2(vec2 < Inf);
    xlim([0 1]);
    ylim([0 1]);    
    % Compute correlation between SCF's and annotate
    corr = 1 - pdist2(vec1', vec2', 'correlation');
    text(1, 0, append("r = $",num2str(corr, 2), "$"), 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    if count == 1
        text(corr_array_hom(position_b_grid(1), position_b_grid(2), 1), ...
            corr_array_het(position_b_grid(1), position_b_grid(2), 1), '$\gamma_\mathbf{p}(\mathbf{q})$', 'Interpreter', 'latex', 'FontSize', 15, 'Color', [0, 0.5, 0], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    end
    % Compute linear fit and annotate
    coefficients = polyfit(vec1, vec2, 1);
    plot(xlim, polyval(coefficients, xlim), 'Color', [0.6, 0.6, 0.6], 'LineStyle', '--', 'LineWidth', 1);
    % Plot SCF points on top of linear fit
    scatter(vec1, vec2, 2, 'black', 'filled');
    scatter(...
        corr_array_hom(position_b_grid(1), position_b_grid(2), count), ...
        corr_array_het(position_b_grid(1), position_b_grid(2), count), ...
        50, [0 0.5 0], 'filled');  
    if count == 1
        xticks([0 1]);
        yticks([0 1]);
        ax.TickLabelInterpreter = 'latex';
        xlabel('$\gamma_\mathbf{p}^\mathrm{geo}(\mathbf{r})$', 'Interpreter', 'latex', 'FontSize', 15);
        ylabel('$\gamma_\mathbf{p}^\mathrm{hyb}(\mathbf{r})$', 'Interpreter', 'latex', 'FontSize', 15);
        ax.XAxis.FontSize = 10;
        ax.YAxis.FontSize = 10;
        ax.LabelFontSizeMultiplier  = 1.5;
    else
        xticks([]);
        yticks([]);       
    end
    text(0.02, 0.99, sigma_n_titles(count), 'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 20);
    title(titles(count), 'Interpreter', 'latex', 'FontSize', 24, 'Units', 'normalized');
    ax.TitleHorizontalAlignment = 'left';
    hold off;
end

% save as scfvssigma_n.tiff
exportgraphics(gcf, '5_scfvssigma_n.tiff', 'Resolution', 300);
close(fig);

