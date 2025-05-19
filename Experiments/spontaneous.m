%% Create schematic of correlation plot

clear; clc; 
loadparam; 

f = figure;
f.Position = [100, 100, 520, 500];
ax = gca;
set(ax, 'color', 'white');
axis off;
% box on;
% ax.LineWidth = 3;
% ax.Color = 'k';
hold on;

dx = topology.L / topology.Nx;

xticks([]);
yticks([]);
xlim([0 topology.L]);
ylim([0 topology.L]);


%draw circle
% region_radius = 0.02;
xposition = hetparam_het.a(1);
yposition = hetparam_het.a(2);
region_radius = 0.04;
theta = linspace(-pi,pi,1000);
x = region_radius * cos(theta) + xposition;
y = region_radius * sin(theta) + yposition;
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
        stimulationinputspace(i, j) = exp(-0.5*(distx^2 + disty^2) * (dx^2) / (region_radius^2));
    end
end

stimulationinputspace = stimulationinputspace.* randn(topology.Nx);

imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), stimulationinputspace');
set(gca, 'YDir', 'normal');
colormap([linspace(0, 1, 256)' linspace(0, 1, 256)', linspace(1, 1, 256)';
    linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)']);
clim([-max(stimulationinputspace, [], 'all') max(stimulationinputspace, [], 'all')])

% % Draw computational grid
% grid_values = 0.0:0.01:0.4;
% for i = 1:length(grid_values)
%     xline = plot([grid_values(i), grid_values(i)], [0, topology.L], 'LineWidth', 0.1, 'Color', [0.4 0.4 0.4 0.5]);
%     % uistack(xline, 'bottom');
%     yline = plot([0, topology.L], [grid_values(i), grid_values(i)], 'LineWidth', 0.1, 'Color', [0.4 0.4 0.4 0.5]);
%     % uistack(yline, 'bottom');    
% end

scatter(xposition , yposition, 20, 'r', 'filled');
theta = linspace(0,2*pi,1000);
x = region_radius * cos(theta) + xposition;
y = region_radius * sin(theta) + yposition;
fill(x, y, 'r', 'FaceAlpha', 0., 'EdgeColor', 'r', 'LineWidth', 1);


quiver(xposition - region_radius, yposition - 2*region_radius, ...
2*region_radius, 0, ...
'Color', 'k', 'LineWidth', 1, ...
'MaxHeadSize', 0.3 / norm(hetparam_het.a - hetparam_het.b), ...
'Marker', '.', 'MarkerSize', 0.0001, ...
'AutoScale','off');
quiver(xposition + region_radius, yposition - 2*region_radius, ...
-2*region_radius, 0, ...
'Color', 'k', 'LineWidth', 1, ...
'MaxHeadSize', 0.3 / norm(hetparam_het.a - hetparam_het.b), ...
'Marker', '.', 'MarkerSize', 0.0001, ...
'AutoScale','off');
text(xposition, yposition - 2*region_radius - 0.01, "$2\sigma_n$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','top', 'FontSize', 14);




% Draw point on position_centre and annotate
scatter(hetparam_het.a(1) , hetparam_het.a(2), 40, 'k', 'filled');

text(hetparam_het.a(1), hetparam_het.a(2) - 0.015, "$\mathbf{p}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 18);


% Draw random point
randpt = [0.25, 0.3];
scatter(randpt(1), randpt(2), 20, 'k', 'filled');
text((hetparam_het.a(1) + randpt(1))/2 + 0.01, (hetparam_het.a(2) + randpt(2))/2 - 0.01, ["$\gamma_{\mathbf{p}}(\mathbf{r}) = \mathrm{corr}($","$\hspace{1em} \phi(\mathbf{r}, t), \phi(\mathbf{p}, t))$"], 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','left','VerticalAlignment','middle', 'FontSize', 15);
text(randpt(1), randpt(2) + 0.015, "$\mathbf{r}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 18);

text(xposition - 0.05, yposition + 0.02, "$\mathbf{i}\textbf{.}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 24);
text((hetparam_het.a(1) + randpt(1))/2, (hetparam_het.a(2) + randpt(2))/2 - 0.01, "$\mathbf{ii}\textbf{.}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontSize', 24);


plot([hetparam_het.a(1), randpt(1)], [hetparam_het.a(2), randpt(2)], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

text(topology.L - 0.01, 0, "$\Omega$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontSize', 24);

% Draw border
pos = get(ax, 'Position');
annotation('rectangle', [pos(1), pos(2), pos(1) + pos(3) - pos(1), pos(2) + pos(4) - pos(2)], ...
    'LineWidth', 1, 'Color', 'k');


% save as schematiccorrelation.svg
print(gcf, 'schematiccorrelation.svg', '-dsvg');
exportgraphics(gcf, 'schematiccorrelation.tiff', 'Resolution', 300);

%% Calculate SCF of heterogeneous model vs spatial width of noisy input
% All points in stimulated region are stimulated with white noise

% Preparation

clear; clc;

loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

% Set max_distance of stimulation point from Ri

position_centre = hetparam_het.a;
position_centre_grid = position_centre / dx;

num_experiments = 4;
sigma_n_array = [0.015, 0.05, 0.10, Inf];

corr_array_concat = zeros(topology.Nx, topology.Nx, length(sigma_n_array));

tol = 1e-7;

parfor count = 1:length(sigma_n_array)

    rng(0, "twister");

    % Set region radius

    sigma_n = sigma_n_array(count);

    % Create stim struct

    sum_xy = zeros(topology.Nx, topology.Nx);
    sum_x2 = zeros(topology.Nx, topology.Nx);
    sum_x = zeros(topology.Nx, topology.Nx);
    sum_y2 = 0;
    sum_y = 0;
    init = 0;

    cov_array = zeros(topology.Nx, topology.Nx);
    std_array = zeros(topology.Nx, topology.Nx);
    stdy = 0;
    corr_array = zeros(topology.Nx, topology.Nx);

    % simulate both models and calculate correlation of all points 
    % with centre of Ri, for both models

    delta = Inf;
    tottime = 0;

    stimulationinputspace = zeros(topology.Nx, topology.Nx);
    i0 = position_centre(1) / dx; j0 = position_centre(2) / dx;
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
      
    while(delta > tol || isnan(delta))
        tottime = tottime + topology.Nt;
        stiminput = (1/(sqrt(dt) * dx)) * randn(topology.Nx, topology.Nx, topology.Nt);
        if sigma_n < Inf
            stiminput = stiminput .* reshape(stimulationinputspace, [topology.Nx, topology.Nx, 1]);
        end
        ts = run_init_cstminput_periodic(topology,homparam,hetparam_het,stiminput,init);
        init = ts(:, :, end-1:end);
        ts_centre = ts(position_centre_grid(1), position_centre_grid(2), :);
        sum_y = sum_y + sum(ts_centre);
        sum_y2 = sum_y2 + sum(ts_centre.^2);
        for i = 1:topology.Nx
            for j = 1:topology.Nx
                ts_point = ts(i, j, :);
                sum_xy(i, j) = sum_xy(i, j) + sum(ts_centre.*ts_point);
                sum_x2(i, j) = sum_x2(i, j) + sum(ts_point.*ts_point);
                sum_x(i, j) = sum_x(i, j) + sum(ts_point);
            end
        end
        cov_array = (sum_xy / tottime) - (sum_x / tottime) * (sum_y / tottime);
        std_array = sqrt( (sum_x2 / tottime) - (sum_x / tottime).^2 );
        stdy = sqrt( (sum_y2 / tottime) - (sum_y / tottime).^2 );
        corr_array_old = corr_array;
        corr_array = cov_array ./ (std_array * stdy);
        if all(corr_array(:) == 0) || all(corr_array_old(:) == 0)
            delta = NaN;
        else
            delta = pdist2(corr_array_old(:)', corr_array(:)', 'cosine');
        end
        disp(num2str([sigma_n, delta]));
    end

    corr_array_concat(:, :, count) = corr_array;

end 

save('scfvssigma_n.mat', 'num_experiments', 'sigma_n_array', 'corr_array_concat', 'topology')

%% Calculate SCF of homogeneous model vs spatial width of noisy input
% All points in stimulated region are stimulated with white noise

% Preparation

clear; clc;

loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

% Set max_distance of stimulation point from Ri

position_centre = hetparam_het.a;
position_centre_grid = position_centre / dx;

num_experiments = 4;
sigma_n_array = [0.015, 0.05, 0.10, Inf];

corr_array_concat = zeros(topology.Nx, topology.Nx, length(sigma_n_array));

tol = 1e-7;

parfor count = 1:length(sigma_n_array)

    rng(0, "twister");

    % Set region radius

    sigma_n = sigma_n_array(count);

    % Create stim struct

    sum_xy = zeros(topology.Nx, topology.Nx);
    sum_x2 = zeros(topology.Nx, topology.Nx);
    sum_x = zeros(topology.Nx, topology.Nx);
    sum_y2 = 0;
    sum_y = 0;
    init = 0;

    cov_array = zeros(topology.Nx, topology.Nx);
    std_array = zeros(topology.Nx, topology.Nx);
    stdy = 0;
    corr_array = zeros(topology.Nx, topology.Nx);

    % simulate both models and calculate correlation of all points 
    % with centre of Ri, for both models

    delta = Inf;
    tottime = 0;

    stimulationinputspace = zeros(topology.Nx, topology.Nx);
    i0 = position_centre(1) / dx; j0 = position_centre(2) / dx;
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
      
    while(delta > tol || isnan(delta))
        tottime = tottime + topology.Nt;
        stiminput = (1/(sqrt(dt) * dx)) * randn(topology.Nx, topology.Nx, topology.Nt);
        if sigma_n < Inf
            stiminput = stiminput .* reshape(stimulationinputspace, [topology.Nx, topology.Nx, 1]);
        end
        ts = run_init_cstminput_periodic(topology,homparam,hetparam_hom,stiminput,init);
        init = ts(:, :, end-1:end);
        ts_centre = ts(position_centre_grid(1), position_centre_grid(2), :);
        sum_y = sum_y + sum(ts_centre);
        sum_y2 = sum_y2 + sum(ts_centre.^2);
        for i = 1:topology.Nx
            for j = 1:topology.Nx
                ts_point = ts(i, j, :);
                sum_xy(i, j) = sum_xy(i, j) + sum(ts_centre.*ts_point);
                sum_x2(i, j) = sum_x2(i, j) + sum(ts_point.*ts_point);
                sum_x(i, j) = sum_x(i, j) + sum(ts_point);
            end
        end
        cov_array = (sum_xy / tottime) - (sum_x / tottime) * (sum_y / tottime);
        std_array = sqrt( (sum_x2 / tottime) - (sum_x / tottime).^2 );
        stdy = sqrt( (sum_y2 / tottime) - (sum_y / tottime).^2 );
        corr_array_old = corr_array;
        corr_array = cov_array ./ (std_array * stdy);
        if all(corr_array(:) == 0) || all(corr_array_old(:) == 0)
            delta = NaN;
        else
            delta = pdist2(corr_array_old(:)', corr_array(:)', 'cosine');
        end
        disp(num2str([sigma_n, delta]));
    end

    corr_array_concat(:, :, count) = corr_array;

end 

save('scfvssigma_n_homogeneous.mat', 'num_experiments', 'sigma_n_array', 'corr_array_concat', 'topology')

%% Plot SCF of heterogeneous model vs input width
clear; clc;
loadparam;

dx = topology.L / topology.Nx;

load('scfvssigma_n.mat')
corr_array_concat_het = corr_array_concat;
load("scfvssigma_n_homogeneous.mat")
corr_array_concat_hom = corr_array_concat;
clear corr_array_concat;

position_centre = hetparam_het.a;
position_centre_grid = position_centre / dx;
for count = 1:length(sigma_n_array)
    corr_array_concat_het(position_centre_grid(1), position_centre_grid(2), count) = Inf;
    corr_array_concat_hom(position_centre_grid(1), position_centre_grid(2), count) = Inf;
end
        
f = figure;
f.Position = [100, 50, 1170, 910];
num_subcols = 8;
num_subrows = 4;
tot_subcols = num_subcols*num_experiments;
tot_subrows = 2 + 2*num_subrows;
t = tiledlayout(tot_subrows, tot_subcols, 'TileSpacing', 'loose', 'Padding', 'tight');

titles = "\textbf{" + ["i", "ii", "iii", "iv"] + ".}";
sigma_n_titles = "$\sigma_n=" + ["3L/80", "L/8", "L/4", "\infty"] + "$";

nexttile(1, [1, tot_subcols])
axis off;
text(0, 1, "$\gamma_\mathbf{p}^\mathrm{hyb}(\mathbf{r})$ over $\Omega$ - varying $\sigma_n$", 'Interpreter', 'latex', 'FontSize', 24, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')

for count = 1:num_experiments
    sigma_n = sigma_n_array(count);
    ax = nexttile(1 + tot_subcols + (count - 1)*num_subcols, [num_subrows num_subcols]);
    box on;
    ax.LineWidth = 1;
    ax.Color = 'k';
    hold on;
    imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), corr_array_concat_het(:, :, count)');
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

% nexttile(2*tot_subcols, [num_subrows, 1]);
% axis off;
% cb = colorbar;
% cb.FontSize = 15;
% cb.TickLabelInterpreter = 'latex';
% title(cb, "$\gamma_{\mathbf{p}}^{\mathrm{II}}(\mathbf{r})$", 'Interpreter', 'latex', 'Position', [7 195 0]);

nexttile(1 + tot_subcols * (1 + num_subrows), [1, tot_subcols])
axis off;
text(0, 1, "$(\gamma_\mathbf{p}^\mathrm{geo}(\mathbf{r}), \gamma_\mathbf{p}^\mathrm{hyb}(\mathbf{r}))$ - varying $\sigma_n$", 'Interpreter', 'latex', 'FontSize', 24, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')


position_rf_grid = hetparam_het.b / dx;

for count = 1:num_experiments
    ax = nexttile(1 + (2 + num_subrows) * tot_subcols + (count - 1)*num_subcols, [num_subrows num_subcols]);
    hold;
    vec1 = reshape(corr_array_concat_hom(:, :, count), [], 1);
    vec2 = reshape(corr_array_concat_het(:, :, count), [], 1);
    vec1 = vec1(vec1 < Inf); vec2 = vec2(vec2 < Inf);
    coefficients = polyfit(vec1, vec2, 1);
    corr = 1 - pdist2(vec1', vec2', 'correlation');
    text(1, 0, append("r = $",num2str(corr, 2), "$"), 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    if count == 1
        text(corr_array_concat_hom(position_rf_grid(1), position_rf_grid(2), 1), ...
            corr_array_concat_het(position_rf_grid(1), position_rf_grid(2), 1), '$\gamma_\mathbf{p}(\mathbf{q})$', 'Interpreter', 'latex', 'FontSize', 15, 'Color', [0, 0.5, 0], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    end
    xlim([0 1]);
    ylim([0 1]);
    plot(xlim, polyval(coefficients, xlim), 'Color', [0.6, 0.6, 0.6], 'LineStyle', '--', 'LineWidth', 1);
    scatter(vec1, vec2, 2, 'black', 'filled');
    scatter(...
        corr_array_concat_hom(position_rf_grid(1), position_rf_grid(2), count), ...
        corr_array_concat_het(position_rf_grid(1), position_rf_grid(2), count), ...
        50, [0 0.5 0], 'filled');
    % if count == 1
    %     scatter(...
    %         corr_array_concat_hom(position_centre_grid(1), position_centre_grid(2) + [-1 1], count), ...
    %         corr_array_concat_het(position_centre_grid(1), position_centre_grid(2) + [-1 1], count), ...
    %         50, 'b', 'filled');
    %     scatter(...
    %         corr_array_concat_hom(position_centre_grid(1) + [-1 1], position_centre_grid(2), count), ...
    %         corr_array_concat_het(position_centre_grid(1) + [-1 1], position_centre_grid(2), count), ...
    %         50, 'b', 'filled');
    % end    
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

% save as correlationvsstimradius.svg
print(gcf, 'scfvssigma_n.svg', '-dsvg');
exportgraphics(gcf, 'scfvssigma_n.tiff', 'Resolution', 300);

%% Plot SCF of homogeneous model vs input width
clear; clc;
loadparam;

dx = topology.L / topology.Nx;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

load("scfvssigma_n_homogeneous.mat")

position_centre = hetparam_het.a;
position_centre_grid = position_centre / dx;
for count = 1:length(sigma_n_array)
    corr_array_concat(position_centre_grid(1), position_centre_grid(2), count) = Inf;
end
        
f = figure;
f.Position = [100, 50, 1170, 465];
num_subrows = 4;
tot_subrows = 2 + num_subrows;
t = tiledlayout(tot_subrows, num_experiments, 'TileSpacing', 'compact', 'Padding', 'tight');

titles = "\textbf{" + ["i", "ii", "iii", "iv"] + ".}";
sigma_n_titles = "$\sigma_n=" + ["3L/80", "L/8", "L/4", "\infty"] + "$";

nexttile(1, [1, num_experiments])
axis off;
text(0, 1, "$\gamma_\mathbf{p}^\mathrm{geo}(\mathbf{r})$ over $\Omega$ - varying $\sigma_n$", 'Interpreter', 'latex', 'FontSize', 24, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')

for count = 1:num_experiments
    sigma_n = sigma_n_array(count);
    ax = nexttile(num_experiments + count, [num_subrows, 1]);
    box on;
    ax.LineWidth = 1;
    ax.Color = 'k';
    hold on;
    imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), corr_array_concat(:, :, count)');
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

% save as correlationvsstimradius_homogeneous.svg
print(gcf, 'scfvssigma_n_homogeneous.svg', '-dsvg');
exportgraphics(gcf, 'scfvssigma_n_homogeneous.tiff', 'Resolution', 300);