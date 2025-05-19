%% Animate geometric model

clear; clc;

loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

Colormap = [linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];

ts = run_periodic(topology,homparam,hetparam_hom,stim);

maxval = 650;
%
f = figure;
sgtitle('Evoked Response - Geometric Connectivity')
f.Position = [300 100 650 600];
set(gcf, 'color', 'white');

formatSpec = '%.0f';

h = imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), ts(:,:,1)');
hold;
cross = scatter3([stim.stimR(1)], [stim.stimR(2)], 0, 400, 'Marker', '.', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2);
xlabel(append('t = -5 ms'));
view(0,90);
colormap(Colormap);
shading flat;
clim([0 maxval]);
xlim([0 topology.L]);
ylim([0 topology.L]);
xticks([]);
yticks([]);

pause;

% gif('Homogeneous.gif','DelayTime',1/24)
while(1)
    for count = round((stim.stimt - 0.005) / dt):1:topology.Nt
        time = count*dt - stim.stimt;
        set(h, 'Cdata', ts(:,:,count)');
        if time < 0
            set(cross, 'MarkerEdgeColor', 'k')
        else
            set(cross, 'MarkerEdgeColor', 'none');
        end
        xlabel(append('t = ', num2str(1000*time, formatSpec), ' ms'));
        pause(0.02);
        % gif;
    end
end

%% Animate hybrid model (geometric connectivity + long-range connectivity)

clear; clc;

loadparam;

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
hetparam_het.m = 2*num_fnps;
hetparam_het.c = (homparam.r)^2 * ones(1, 2*num_fnps);
hetparam_het.tau = zeros(1, 2*num_fnps);
hetparam_het.a = [a_array, b_array];
hetparam_het.b = [b_array, a_array];

Colormap = [linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];

ts = run_periodic(topology,homparam,hetparam_het,stim);

maxval = 400;
%
f = figure;
sgtitle('Evoked Response - Hybrid Connectivity')
f.Position = [300 100 650 600];
set(gcf, 'color', 'white');

formatSpec = '%.0f';

h = imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), ts(:,:,1)');
hold;
cross = scatter3([stim.stimR(1)], [stim.stimR(2)], 0, 400, 'Marker', '.', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2);

% Draw FNPs
for k = 1:hetparam_het.m
    quiver(hetparam_het.a(1, k), hetparam_het.a(2, k), ...
    hetparam_het.b(1, k) - hetparam_het.a(1, k), ...
    hetparam_het.b(2, k) - hetparam_het.a(2, k),...
    'Color', 'k', 'LineWidth', 1, ...
    'MaxHeadSize', 0.0 / norm(hetparam_het.a(:, k) - hetparam_het.b(:, k)), ...
    'Marker', '.', 'MarkerSize', 0.0001, ...
    'AutoScale','off');
end


xlabel(append('t = -5 ms'));
view(0,90);
colormap(Colormap);
shading flat;
clim([0 maxval]);
xlim([0 topology.L]);
ylim([0 topology.L]);
xticks([]);
yticks([]);

pause;

% gif('Homogeneous.gif','DelayTime',1/24)
while(1)
    for count = round((stim.stimt - 0.005) / dt):1:topology.Nt
        time = count*dt - stim.stimt;
        set(h, 'Cdata', ts(:,:,count)');
        if time < 0
            set(cross, 'MarkerEdgeColor', 'k')
        else
            set(cross, 'MarkerEdgeColor','none');
        end
        xlabel(append('t = ', num2str(1000*time, formatSpec), ' ms'));
        pause(0.02);
        % gif;
    end
end

%% A. Animate underlying evoked response of homogeneous model, and heterogeneous model with 10 and 20 LRCs

clear; clc;

loadparam;

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

num_lrcs = 50;

% Create LRC topology
rng(1, "twister");
a_array = zeros(2, num_lrcs);
b_array = zeros(2, num_lrcs);

% Sample connectome

for j = 1:num_lrcs
    
    a = topology.L*rand(2, 1);
    b = topology.L*rand(2, 1);
    a_array(:, j) = a;
    b_array(:, j) = b;
    
end

hetparam_het.m = num_lrcs;
hetparam_het.c = (homparam.r)^2 * ones(1, num_lrcs);
hetparam_het.tau = zeros(1, num_lrcs);
hetparam_het.a = a_array;
hetparam_het.b = b_array;

hetparam_array = {hetparam_hom, hetparam_het};

ts_full = zeros(topology.Nx, topology.Nx, topology.Nt, 2);
ts_full(:, :, :, 1) = run_periodic(topology,homparam,hetparam_hom,stim);
ts_full(:, :, :, 2) = run_periodic(topology,homparam,hetparam_het,stim);
%%
%maxval = max(ts_full(:, :, end), [], 'all') * 3;
maxval = 600;
timepoints = 0:0.002:0.010;
num_snapshots = length(timepoints);


%load 'CustomColormap.mat'
%Colormap = CustomColormap;
Colormap = [linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

f = figure;
f.Position = [100, 100, 1400, 400];

num_subcols = 5;
num_subrows = 10;
tot_subcols = 3+num_subcols*(1+num_snapshots);
tot_subrows = 1+num_subrows*2;

t = tiledlayout(tot_subrows,tot_subcols,"TileSpacing", "tight", "Padding", "compact");

modeltitles = ["Model " + ["I", "II"]; ["Isotropic", "Perturbed"]];
stimulustitles = "\bf" + ["i", "ii"] + ".";
responsetitles = "\bf" + ["iii", "iv"] + ".";
coltitles = strings(1, 1+num_snapshots);
coltitles(1) = "\bf \rm Stimulus";
for i = 1:num_snapshots
    coltitles(i+1) = append("$t = ", num2str(1000*timepoints(i)), "\ \mathrm{ms}$");
end

for iter = 1:2
    hetparam_het = hetparam_array{iter};
    ts = squeeze(ts_full(:, :, :, iter));
    % Label Heterogeneous or Homogeneous Model
    nexttile(1 + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, 1]);
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(0.0, 0.5, modeltitles(1, iter), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'Units', 'normalized', 'Rotation', 90, 'FontWeight','bold');
    text(0.5, 0.5, modeltitles(2, iter), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'Units', 'normalized', 'Rotation', 90, 'FontWeight','bold');
    % Label Stimulus Subplot
    nexttile(2 + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, 1])
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(1, 0.5, stimulustitles(iter), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 16, 'Units', 'normalized', 'FontWeight','bold');
    % Plot Stimulus Subplot
    ax = nexttile(3 + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, num_subcols]);
    hold on;
    ax.Box = "on";
    ax.LineWidth = 1;
    set(gca, 'color', 'white');
    for k = 1:hetparam_het.m
        quiver(hetparam_het.a(1, k) + dx/2, hetparam_het.a(2, k) + dx/2, ...
        hetparam_het.b(1, k) - hetparam_het.a(1, k), ...
        hetparam_het.b(2, k) - hetparam_het.a(2, k),...
        'Color', 'k', 'LineWidth', 1, ...
        'MaxHeadSize', 0.05 / norm(hetparam_het.a(:, k) - hetparam_het.b(:, k)), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end
    scatter(stim.stimR(1) + dx/2 , stim.stimR(2) + dx/2, 15, 'r', 'filled');
    if iter == 1
        text(stim.stimR(1) + dx/2 + 0.03, stim.stimR(2) + dx/2 - 0.03, "$\mathbf{r}_0$", 'Interpreter', 'latex', 'Color', 'r', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 15);
    end
    view(0,90);
    xlim([0, topology.L + dx]);
    ylim([0, topology.L + dx]);
    xticks([]);
    yticks([]);
    if iter == 1
        title(coltitles(1), 'FontSize', 12);
    end
    hold off;
    % Label Response Subplot
    nexttile(3 + num_subcols + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, 1])
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(1, 0.5, responsetitles(iter), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 16, 'Units', 'normalized', 'FontWeight','bold');
    % Plot Response Subplot
    for i = 1:num_snapshots
        ax = nexttile(4 + num_subcols + (iter - 1)*(1+num_subrows)*tot_subcols + (i - 1)*num_subcols, [num_subrows, num_subcols]);
        hold on;
        ax.Box = "on";
        ax.LineWidth = 1;
        imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), ts(:,:,ceil((stim.stimt + timepoints(i)) / dt))');
        set(ax,'YDir','normal');
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
        if iter == 1
            title(coltitles(i+1), 'Interpreter', 'latex', 'FontSize', 12);
        end
        hold off;
    end
end

cb = colorbar; 
cb.Layout.Tile = 'east'; 


%% Animate underlying evoked response of 
% homogeneous model, and heterogeneous model with 10 and 50 recip LRCs

clear; clc;

loadparam;

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

num_lrcs_array = [50, 100];
max_num_lrcs = max(num_lrcs_array);

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

ts_full = zeros(topology.Nx, topology.Nx, topology.Nt, 1 + length(num_lrcs_array));

ts_full(:, :, :, 1) = run_periodic(topology,homparam,hetparam_hom,stim);

hetparam_het.sigmaeps = 0.01;
for i = 1:length(num_lrcs_array)
    num_lrcs = num_lrcs_array(i);
    hetparam_het.m = num_lrcs;
    hetparam_het.c = (homparam.r)^2 * ones(1, num_lrcs);
    hetparam_het.tau = zeros(1, num_lrcs);
    hetparam_het.a = a_array(:, 1:num_lrcs);
    hetparam_het.b = b_array(:, 1:num_lrcs);
 
    ts_full(:, :, :, 1 + i) = run_periodic(topology,homparam,hetparam_het,stim);
end
%%
%maxval = max(ts_full(:, :, end), [], 'all') * 3;
maxval = 600;
timepoints = [0, 0.005, 0.01, 0.015, 0.02, 0.05];
num_snapshots = length(timepoints);


%load 'CustomColormap.mat'
%Colormap = CustomColormap;
Colormap = [linspace(0, 1, 256)' linspace(0, 1, 256)', linspace(1, 1, 256)';
    linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];

f = figure;
f.Position = [100, 100, 1400, 400];

num_subcols = 5;
num_subrows = 10;
tot_subcols = 3+num_subcols*(1+num_snapshots);
tot_subrows = 2+num_subrows*3;

t = tiledlayout(tot_subrows,tot_subcols,"TileSpacing", "tight", "Padding", "compact");

modeltitles = ["Model " + ["I", "II", "III"]; ["Isotropic", "Perturbed - " + string(num_lrcs_array) + " tracts"]];
stimulustitles = "\bf" + ["i", "ii", "iii"] + ".";
responsetitles = "\bf" + ["iv", "v", "vi"] + ".";
coltitles = strings(1, 1+num_snapshots);
coltitles(1) = "\bf \rm Stimulus";
for i = 1:num_snapshots
    coltitles(i+1) = append("$t = ", num2str(1000*timepoints(i)), "\ \mathrm{ms}$");
end

for iter = 1:(1+length(num_lrcs_array))
    if iter == 1
        num_lrcs = 0;
    else
        num_lrcs = num_lrcs_array(iter-1);
    end
    ts = squeeze(ts_full(:, :, :, iter));
    % Label Heterogeneous or Homogeneous Model
    nexttile(1 + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, 1]);
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(0.0, 0.5, modeltitles(1, iter), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'Units', 'normalized', 'Rotation', 90, 'FontWeight','bold');
    text(0.5, 0.5, modeltitles(2, iter), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'Units', 'normalized', 'Rotation', 90, 'FontWeight','bold');
    % Label Stimulus Subplot
    nexttile(2 + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, 1])
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(1, 0.5, stimulustitles(iter), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 16, 'Units', 'normalized', 'FontWeight','bold');
    % Plot Stimulus Subplot
    ax = nexttile(3 + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, num_subcols]);
    hold on;
    ax.Box = "on";
    ax.LineWidth = 1;
    set(gca, 'color', 'white');
    for m = 1:num_lrcs
        quiver(dx/2 + a_array(1, m), dx/2 + a_array(2, m), ...
        b_array(1, m) - a_array(1, m), b_array(2, m) - a_array(2, m),...
        'Color', 'k', 'LineWidth', 1, ...
        'MaxHeadSize', 0.05 / norm(a_array(:, i) - b_array(:, i)), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end
    scatter(stim.stimR(1) + dx/2 , stim.stimR(2) + dx/2, 15, 'r', 'filled');
    if iter == 1
        text(stim.stimR(1) + dx/2 + 0.03, stim.stimR(2) + dx/2 - 0.03, "$\mathbf{r}_0$", 'Interpreter', 'latex', 'Color', 'r', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 15);
    end
    view(0,90);
    xlim([0, topology.L + dx]);
    ylim([0, topology.L + dx]);
    xticks([]);
    yticks([]);
    if iter == 1
        title(coltitles(1), 'FontSize', 12);
    end
    hold off;
    % Label Response Subplot
    nexttile(3 + num_subcols + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, 1])
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(1, 0.5, responsetitles(iter), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 16, 'Units', 'normalized', 'FontWeight','bold');
    % Plot Response Subplot
    for i = 1:num_snapshots
        ax = nexttile(4 + num_subcols + (iter - 1)*(1+num_subrows)*tot_subcols + (i - 1)*num_subcols, [num_subrows, num_subcols]);
        hold on;
        ax.Box = "on";
        ax.LineWidth = 1;
        imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), ts(:,:,ceil((stim.stimt + timepoints(i)) / dt))');
        set(ax,'YDir','normal');
        colormap(ax, Colormap);
        shading flat;
        clim([-maxval maxval]);
        view(0,90);
        xlim([0, topology.L + dx]);
        ylim([0, topology.L + dx]);
        xticks([]);
        yticks([]);
        if iter == 1
            title(coltitles(i+1), 'Interpreter', 'latex', 'FontSize', 12);
        end
        hold off;
    end
end

cb = colorbar; 
cb.Layout.Tile = 'east'; 

%% A. Visualize default model optimal stimulus with snapshots over time

clear; clc;
loadparam;

% Optimal stimulation
stim.stimR = hetparam_het.a;

%stim.sigma(2) = stim.sigma(2)*0.;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

ts_full = zeros(topology.Nx, topology.Nx, topology.Nt, 2);
ts_full(:, :, :, 1) = run_periodic(topology,homparam,hetparam_hom,stim);
ts_full(:, :, :, 2) = run_periodic(topology,homparam,hetparam_het,stim);

%
%maxval = max(ts_full(:, :, end), [], 'all') * 3;
maxval = 600;
timepoints = [-0.0008, 0.005, 0.01, 0.02, 0.03, 0.05];
num_snapshots = length(timepoints);

%load 'CustomColormap.mat'
%Colormap = CustomColormap;
Colormap = [linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

f = figure;
f.Position = [100, 100, 1400, 400];

num_subcols = 5;
num_subrows = 10;
tot_subcols = 3+num_subcols*(1+num_snapshots);
tot_subrows = 1+num_subrows*2;

t = tiledlayout(tot_subrows,tot_subcols,"TileSpacing", "tight", "Padding", "compact");

modeltitles = ["Model " + ["I", "II"]; ["Isotropic", "Perturbed"]];
stimulustitles = "\bf" + ["i", "ii"] + ".";
responsetitles = "\bf" + ["iii", "iv"] + ".";
coltitles(1) = "\bf \rm Stimulus";
coltitles(2) = append("$t \approx 0 \ \mathrm{ms}$");
for i = 2:num_snapshots
    coltitles(i+1) = append("$t = ", num2str(1000*timepoints(i)), "\ \mathrm{ms}$");
end

for iter = 1:2
    ts = squeeze(ts_full(:, :, :, iter));
    % Label Heterogeneous or Homogeneous Model
    nexttile(1 + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, 1]);
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(0.0, 0.5, modeltitles(1, iter), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'Units', 'normalized', 'Rotation', 90, 'FontWeight','bold');
    text(0.5, 0.5, modeltitles(2, iter), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'Units', 'normalized', 'Rotation', 90, 'FontWeight','bold');
    % Label Stimulus Subplot
    nexttile(2 + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, 1])
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(1, 0.5, stimulustitles(iter), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 16, 'Units', 'normalized', 'FontWeight','bold');
    % Plot Stimulus Subplot
    ax = nexttile(3 + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, num_subcols]);
    hold on;
    ax.Box = "on";
    ax.LineWidth = 1;
    set(gca, 'color', 'white');
    if iter == 2
        quiver(hetparam_het.a(1) + dx/2, hetparam_het.a(2) + dx/2, ...
        hetparam_het.b(1) - hetparam_het.a(1), ...
        hetparam_het.b(2) - hetparam_het.a(2),...
        'Color', 'k', 'LineWidth', 1, ...
        'MaxHeadSize', 0.05 / norm(hetparam_het.a - hetparam_het.b), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end
    scatter(stim.stimR(1) + dx/2 , stim.stimR(2) + dx/2, 15, 'r', 'filled');
    if iter == 2
        text(hetparam_het.a(1) + dx/2 - 0.02, hetparam_het.a(2) + dx/2 - 0.02, "$\mathbf{p}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 15);
        text(hetparam_het.b(1) + dx/2 + 0.02, hetparam_het.b(2) + dx/2 + 0.02, "$\mathbf{q}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 15);            
    end
    view(0,90);
    xlim([0, topology.L + dx]);
    ylim([0, topology.L + dx]);
    xticks([]);
    yticks([]);
    if iter == 1
        title(coltitles(1), 'FontSize', 12);
    end
    hold off;
    % Label Response Subplot
    nexttile(3 + num_subcols + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, 1])
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(1, 0.5, responsetitles(iter), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 16, 'Units', 'normalized', 'FontWeight','bold');
    % Plot Response Subplot
    for i = 1:num_snapshots
        ax = nexttile(4 + num_subcols + (iter - 1)*(1+num_subrows)*tot_subcols + (i - 1)*num_subcols, [num_subrows, num_subcols]);
        hold on;
        ax.Box = "on";
        ax.LineWidth = 1;
        imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), ts(:,:,ceil((stim.stimt + timepoints(i)) / dt))');
        set(ax,'YDir','normal');
        colormap(ax, Colormap);
        shading flat;
        clim([0 maxval]);
        view(0,90);
        xlim([0, topology.L + dx]);
        ylim([0, topology.L + dx]);
        xticks([]);
        yticks([]);
        if iter == 1
            title(coltitles(i+1), 'Interpreter', 'latex', 'FontSize', 12);
        end
        hold off;
    end
end

cb = colorbar; 
cb.Layout.Tile = 'east';

% save as snapshotsdefaultoptimal.svg

%% C. Plot dissimilarity of default model over time, optimal stimulus

clear; clc;
loadparam;

%Optimal Stimulation
stim.stimR = hetparam_het.a(:, 1);

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

time_array = dt * (1:topology.Nt);

ts_hom = run_periodic(topology,homparam,hetparam_hom,stim);
ts_het = run_periodic(topology,homparam,hetparam_het,stim);

dissim_time = zeros(1, topology.Nt);
for n = 1:topology.Nt
    dissim_time(n) = pdist2(...
        reshape(ts_hom(:, :, n), [], 1)', reshape(ts_het(:, :, n), [], 1)', 'cosine');
end


% Plot divergence of default model and random models on same plot

f = figure;
f.Position = [100, 100, 500, 250];
hold;

% Optimal Stimulation
plot(1000*(time_array - stim.stimt), dissim_time', 'LineWidth', 1);

ylim([-0.001, 0.09]);
xlim([-5, Inf]);

xlabel("$t \ (\mathrm{ms})$", 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$C_{\phi}$', 'Interpreter', 'latex', 'FontSize', 20)

% save as dissimcurveoptimal.svg

%% A. Visualize default model nonoptimal stimulus with snapshots over time

clear; clc;
loadparam;

% Nonoptimal stimulation
stim.stimR = [0.25; 0.15];

%stim.sigma(2) = stim.sigma(2)*0.;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

ts_full = zeros(topology.Nx, topology.Nx, topology.Nt, 2);
ts_full(:, :, :, 1) = run_periodic(topology,homparam,hetparam_hom,stim);
ts_full(:, :, :, 2) = run_periodic(topology,homparam,hetparam_het,stim);

%
%maxval = max(ts_full(:, :, end), [], 'all') * 3;
maxval = 600;
timepoints = [-0.0008, 0.005, 0.01, 0.02, 0.03, 0.05];
num_snapshots = length(timepoints);

%load 'CustomColormap.mat'
%Colormap = CustomColormap;
Colormap = [linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

f = figure;
f.Position = [100, 100, 1400, 400];

num_subcols = 5;
num_subrows = 10;
tot_subcols = 3+num_subcols*(1+num_snapshots);
tot_subrows = 1+num_subrows*2;

t = tiledlayout(tot_subrows,tot_subcols,"TileSpacing", "tight", "Padding", "compact");

modeltitles = ["Model " + ["I", "II"]; ["Isotropic", "Perturbed"]];
stimulustitles = "\bf" + ["i", "ii"] + ".";
responsetitles = "\bf" + ["iii", "iv"] + ".";
coltitles(1) = "\bf \rm Stimulus";
coltitles(2) = append("$t \approx 0 \ \mathrm{ms}$");
for i = 2:num_snapshots
    coltitles(i+1) = append("$t = ", num2str(1000*timepoints(i)), "\ \mathrm{ms}$");
end

for iter = 1:2
    ts = squeeze(ts_full(:, :, :, iter));
    % Label Heterogeneous or Homogeneous Model
    nexttile(1 + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, 1]);
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(0.0, 0.5, modeltitles(1, iter), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'Units', 'normalized', 'Rotation', 90, 'FontWeight','bold');
    text(0.5, 0.5, modeltitles(2, iter), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'Units', 'normalized', 'Rotation', 90, 'FontWeight','bold');
    % Label Stimulus Subplot
    nexttile(2 + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, 1])
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(1, 0.5, stimulustitles(iter), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 16, 'Units', 'normalized', 'FontWeight','bold');
    % Plot Stimulus Subplot
    ax = nexttile(3 + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, num_subcols]);
    hold on;
    ax.Box = "on";
    ax.LineWidth = 1;
    set(gca, 'color', 'white');
    if iter == 2
        quiver(hetparam_het.a(1) + dx/2, hetparam_het.a(2) + dx/2, ...
        hetparam_het.b(1) - hetparam_het.a(1), ...
        hetparam_het.b(2) - hetparam_het.a(2),...
        'Color', 'k', 'LineWidth', 1, ...
        'MaxHeadSize', 0.05 / norm(hetparam_het.a - hetparam_het.b), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end
    scatter(stim.stimR(1) + dx/2 , stim.stimR(2) + dx/2, 15, 'r', 'filled');
    if iter == 2
        text(hetparam_het.a(1) + dx/2 - 0.02, hetparam_het.a(2) + dx/2 - 0.02, "$\mathbf{p}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 15);
        text(hetparam_het.b(1) + dx/2 + 0.02, hetparam_het.b(2) + dx/2 + 0.02, "$\mathbf{q}$", 'Interpreter', 'latex', 'Color', 'k', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', 15);            
    end
    view(0,90);
    xlim([0, topology.L + dx]);
    ylim([0, topology.L + dx]);
    xticks([]);
    yticks([]);
    if iter == 1
        title(coltitles(1), 'FontSize', 12);
    end
    hold off;
    % Label Response Subplot
    nexttile(3 + num_subcols + (iter - 1)*(1+num_subrows)*tot_subcols, [num_subrows, 1])
    axis off;
    set(gca, 'Color', 'none');
    xticks([]); yticks([]);
    text(1, 0.5, responsetitles(iter), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 16, 'Units', 'normalized', 'FontWeight','bold');
    % Plot Response Subplot
    for i = 1:num_snapshots
        ax = nexttile(4 + num_subcols + (iter - 1)*(1+num_subrows)*tot_subcols + (i - 1)*num_subcols, [num_subrows, num_subcols]);
        hold on;
        ax.Box = "on";
        ax.LineWidth = 1;
        imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), ts(:,:,ceil((stim.stimt + timepoints(i)) / dt))');
        set(ax,'YDir','normal');
        colormap(ax, Colormap);
        shading flat;
        clim([0 maxval]);
        view(0,90);
        xlim([0, topology.L + dx]);
        ylim([0, topology.L + dx]);
        xticks([]);
        yticks([]);
        if iter == 1
            title(coltitles(i+1), 'Interpreter', 'latex', 'FontSize', 12);
        end
        hold off;
    end
end

cb = colorbar; 
cb.Layout.Tile = 'east';

% save as snapshotsdefaultoptimal.svg

%% C. Plot dissimilarity of default model over time, optimal stimulus

clear; clc;
loadparam;



dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

time_array = dt * (1:topology.Nt);

dissim_time = zeros(2, topology.Nt);

for iter = 1:2

    if iter == 1
        %Optimal Stimulation
        stim.stimR = hetparam_het.a(:, 1);
    else
        stim.stimR = [0.25; 0.15];
    end
    
    ts_hom = run_periodic(topology,homparam,hetparam_hom,stim);
    ts_het = run_periodic(topology,homparam,hetparam_het,stim);

    for n = 1:topology.Nt
        dissim_time(iter, n) = pdist2(...
            reshape(ts_hom(:, :, n), [], 1)', reshape(ts_het(:, :, n), [], 1)', 'cosine');
    end
end


% Plot divergence of default model and random models on same plot

f = figure;
f.Position = [100, 100, 500, 250];
hold;


for iter = 1:2
    plot(1000*(time_array - stim.stimt), dissim_time(iter, :)', 'LineWidth', 1);
end

ylim([-0.001, 0.09]);
xlim([-5, Inf]);

xlabel("$t \ (\mathrm{ms})$", 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$C_{\phi}$', 'Interpreter', 'latex', 'FontSize', 20)

legend({'Stimulus 1 (Proximate)', 'Stimulus 2 (Distant)'});

% save as dissimcurveoptimal.svg


%% Plot dissimilarity curve of rng models over time

clear; clc;
loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

time_array = dt * (1:topology.Nt);

% Plot divergence of default model and random models on same plot

f = figure;
f.Position = [100, 100, 600, 450];
hold;

load('dissimcurveensemble_multilrcs.mat');

meandissimcurve = squeeze(mean(dissimcurve_array, 3));

for i = 1:length(num_lrcs_array)
    for j = 1:num_samples
        currentColor = get(gca, 'ColorOrder');
        currentColor = currentColor(i, :);
        mixedColor = 0.3*currentColor + 0.7*[1 1 1];
        plot(1000*(time_array - stim.stimt), dissimcurve_array(:, i, j), 'Color', mixedColor, 'LineWidth', 0.1);
    end
end

set(gca,'ColorOrderIndex',1);
for i = 1:length(num_lrcs_array)
    plot(1000*(time_array - stim.stimt), meandissimcurve(:, i), 'LineWidth', 2);
end

ylim([-0.001, Inf]);
xlim([-5, Inf]);

xticks([0:10:100]);
yticks([0:0.1:1]);

xlabel("$t \ (\mathrm{ms})$", 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$ C_{\phi} (t) $', 'Interpreter', 'latex', 'FontSize', 16)

Legend = cell(1, (num_samples + 1)*length(num_lrcs_array));
for i = 1:num_samples*length(num_lrcs_array)
    Legend{i} = "";
end
for i = 1:length(num_lrcs_array)
    Legend{num_samples * length(num_lrcs_array) + i} = append('$', num2str(num_lrcs_array(i)), '$ tracts');
end
l = legend(Legend, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 12);

%% Animate underlying evoked response with varying position

clear; clc;

loadparam;

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

ts_full = zeros(topology.Nx, topology.Nx, topology.Nt, 2);

% hetparam_het.sigmaeps = 0.004;
stim.stimR = hetparam_het.a;
ts_full(:, :, :, 1) = run_periodic(topology,homparam,hetparam_het,stim);
stim.stimR = [0.25; 0.15];
ts_full(:, :, :, 2) = run_periodic(topology,homparam,hetparam_het,stim);

%%
timepoints = [-0.0008, 0.005, 0.01, 0.015, 0.02, 0.05];
num_snapshots = length(timepoints);

maxval = 600;
titles = ["Proximate Stimulus", " Distant Stimulus"];

num_subrows = 10;

Colormap = [linspace(0, 1, 256)' linspace(0, 1, 256)', linspace(1, 1, 256)';
    linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];
f = figure;
f.Position = [100, 100, 890, 310];
tiledlayout(2*num_subrows,num_snapshots,"TileSpacing","tight","Padding","compact");
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

for iter = 1:2
    ts = squeeze(ts_full(:, :, :, iter));
    for i = 1:num_snapshots
        ax = nexttile;
        hold on;
        ax.Box = "on";
        ax.LineWidth = 3;
        ax.XColor = 'k';
        ax.YColor = 'k';
        imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), ts(:,:,ceil((stim.stimt + timepoints(i)) / dt))');
        view(0,90);
        colormap(ax, Colormap);
        shading flat;
        clim([-maxval maxval]);
        xlim([0, topology.L]);
        ylim([0, topology.L]);
        xticks([]);
        yticks([]);
        if i == 1
            quiver(dx/2 + hetparam_het.a(1), dx/2 + hetparam_het.a(2), ...
            hetparam_het.b(1) - hetparam_het.a(1), hetparam_het.b(2) - hetparam_het.a(2),...
            'Color', 'k', 'LineWidth', 1, ...
            'MaxHeadSize', 0.05 / norm(hetparam_het.b - hetparam_het.a), ...
            'Marker', '.', 'MarkerSize', 0.0001, ...
            'AutoScale','off');
            ylabel(titles(iter), 'Color', 'k', 'FontSize', 10)
        end
        if iter == 1
            if i > 1
                title(append(num2str(1000*timepoints(i)), " ms"), 'Interpreter', 'latex', 'FontSize', 10);
            else
                title(append('t = ',num2str(0), " ms"), 'Interpreter', 'latex', 'FontSize', 10);
            end
        end
        hold off;
    end
end


%% Create schematics of topological constraints

clear; clc;

loadparam;

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

f = figure;
f.Position = [100, 100, 1400, 155];
set(gca,'Color','white')
box on;
t = tiledlayout(1, 28, 'TileSpacing', 'none', 'Padding', 'compact');

num_lrcs = 50;

% Create Ensemble of Connectomes

rng(1, "twister");

a_array = zeros(2, num_lrcs);
b_array = zeros(2, num_lrcs);

for j = 1:num_lrcs
    a = topology.L * rand(2, 1);
    b = topology.L * rand(2, 1);
    a_array(:, j) = a;
    b_array(:, j) = b;
end

% Create schematics of default topology with EDR

lambda_array = [0, 0.5, 1];

dist_array = zeros(1, num_lrcs);

for i = 1:length(lambda_array)

    rng(1, "twister");

    lambda = lambda_array(i);

    a_array1 = a_array; b_array1 = b_array;
    
    dist_array = zeros(1, num_lrcs);

    % Create connectome

    for j = 1:num_lrcs

        a = a_array1(:, j); b = b_array1(:, j);
        test = 0;
        while (test == 0)
            distx = abs(a(1) - b(1));
            disty = abs(a(2) - b(2));
            dist = norm([distx; disty], 2);
            dist_array(j) = dist;
            test = 1;
            if rand() > (exp(-dist*(100*lambda)))
                test = 0;
                a = topology.L * rand(2, 1);
                b = topology.L * rand(2, 1);                    
            end
        end
        a_array1(:, j) = a; b_array1(:, j) = b;
        dist_array(j) = dist;
    end

    ax = nexttile(3*i - 2, [1 2]);
    hold;
    set(gca, 'Color', 'white');
    ax.Box = "on";
    ax.LineWidth = 1;

    xlim([0 topology.L])
    ylim([0 topology.L])
    xticks([]);
    yticks([]);
    for k = 1:num_lrcs
        quiver(a_array1(1, k), a_array1(2, k), ...
        b_array1(1, k) - a_array1(1, k), ...
        b_array1(2, k) - a_array1(2, k),...
        'Color', 'k', 'LineWidth', 0.5, ...
        'MaxHeadSize', 0.05 / norm(b_array1(:, k) - a_array1(:, k)), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end

    xlabel(append('$\lambda_l = ', num2str(lambda), '$'), 'Interpreter', 'latex', 'FontSize', 12);

    if i == 2
        title('Exponential Distance Rule (EDR)', 'Interpreter', 'latex', 'FontSize', 12);
    end

    if i == 1
        title('i.', 'FontSize', 16, 'HorizontalAlignment', 'left');
        ax.TitleHorizontalAlignment = 'left';
    end


    hold off;

end

for i = 1:2
    nexttile(3*i);
    hold;
    quiver(0.35, 0.5, 0.3, 0, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 1, 'Marker', '.', 'MarkerSize', 0.0001, 'AutoScale', 'off');
    xlim([0 1])
    ylim([0 1])
    xticks([]);
    yticks([]);
    axis off;
end


% Create schematics of default topology with hubs

% Set hub region position
num_hubs = 4;
hub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
hublength = topology.L / sqrt(34);

lambda_array = [0, 0.5, 1];

for i = 1:length(lambda_array)

    rng(1, "twister");

    lambda = lambda_array(i);

    a_array1 = a_array; b_array1 = b_array;

    % Create connectome

    for j = 1:num_lrcs

        a = a_array1(:, j); b = b_array1(:, j);
        test = 0;
        reject_rv = rand();
        while (test == 0)
            test = 1;
            if reject_rv < lambda
                % Check if LRC is not incident with any hub region
                % Calculate L-inf distance of Ri from hub centre
                dist_a = zeros(1, num_hubs);
                for iter0 = 1:num_hubs
                    distx = abs(a(1) - hub_centres(iter0, 1));
                    disty = abs(a(2) - hub_centres(iter0, 2));
                    dist_a(iter0) = norm([distx; disty], Inf);
                end
                % Calculate L-inf distance of Rf from hub centre
                dist_b = zeros(1, num_hubs);
                for iter0 = 1:num_hubs
                    distx = abs(b(1) - hub_centres(iter0, 1));
                    disty = abs(b(2) - hub_centres(iter0, 2));
                    dist_b(iter0) = norm([distx; disty], Inf);
                end
                % Rejection 1: a and b are both inside a hub
                test1 = any(dist_a <= hublength/2) & any(dist_b <= hublength/2);
                % Rejection 2: a and b are both outside a hub
                test2 = all(dist_a > hublength/2) & all(dist_b > hublength/2);
                if test1 || test2
                    test = 0;
                    a = topology.L * rand(2, 1);
                    b = topology.L * rand(2, 1);   
                end
            end
        end
        a_array1(:, j) = a; b_array1(:, j) = b;
    end

    ax = nexttile(10 + 3*i - 2, [1 2]);
    hold;
    set(gca, 'Color', 'white');
    ax.Box = "on";
    ax.LineWidth = 1;

    xlim([0 topology.L])
    ylim([0 topology.L])
    xticks([]);
    yticks([]);
    for k = 1:num_lrcs
        quiver(a_array1(1, k), a_array1(2, k), ...
        b_array1(1, k) - a_array1(1, k), ...
        b_array1(2, k) - a_array1(2, k),...
        'Color', 'k', 'LineWidth', 0.5, ...
        'MaxHeadSize', 0.05 / norm(b_array1(:, k) - a_array1(:, k)), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end

    xlabel(append('$\lambda_h = ', num2str(lambda), '$'), 'Interpreter', 'latex', 'FontSize', 12);

    % Draw hub region
    for k = 1:num_hubs
        fill(hub_centres(k, 1) + hublength/2 * [-1 -1 1 1], hub_centres(k, 2) + hublength/2 * [-1 1 1 -1], 'k', 'FaceColor', '#D95319', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    end
    
    if i == 2
        title('Hub Specificity', 'Interpreter', 'latex', 'FontSize', 12);
    end
    
    if i == 1
        title('ii.', 'FontSize', 16, 'HorizontalAlignment', 'left');
        ax.TitleHorizontalAlignment = 'left';
    end
   
    hold off;

end



for i = 1:2
    nexttile(10 + 3*i);
    hold;
    quiver(0.35, 0.5, 0.3, 0, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 1, 'Marker', '.', 'MarkerSize', 0.0001, 'AutoScale', 'off');
    xlim([0 1])
    ylim([0 1])
    xticks([]);
    yticks([]);
    axis off;
end

% Create schematics of default topology with cores

% Set number of rich club
% Set rich club region position
num_richclubs = 4;
richclub_centres = topology.L * [0.25 0.25; 0.25 0.75; 0.75 0.75; 0.75 0.25];
richclublength = topology.L / sqrt(34);

lambda_array = [0, 0.5, 1];

for i = 1:length(lambda_array)

    rng(1, "twister");

    lambda = lambda_array(i);

    a_array1 = a_array; b_array1 = b_array;

    % Create connectome

    for j = 1:num_lrcs

        a = a_array1(:, j); b = b_array1(:, j);
        test = 0;
        reject_rv = rand();
        while (test == 0)
            distx = abs(a(1) - b(1));
            disty = abs(a(2) - b(2));
            dist = norm([distx; disty], 2);
            test = 1;
            if reject_rv < lambda
                % Check if LRC is not incident with any hub region
                % Calculate L-inf distance of Ri from hub centre
                dist_a = zeros(1, num_richclubs);
                for iter0 = 1:num_richclubs
                    distx = abs(a(1) - richclub_centres(iter0, 1));
                    disty = abs(a(2) - richclub_centres(iter0, 2));
                    dist_a(iter0) = norm([distx; disty], Inf);
                end
                % Calculate L-inf distance of Rf from hub centre
                dist_b = zeros(1, num_richclubs);
                for iter0 = 1:num_richclubs
                    distx = abs(b(1) - richclub_centres(iter0, 1));
                    disty = abs(b(2) - richclub_centres(iter0, 2));
                    dist_b(iter0) = norm([distx; disty], Inf);
                end
                % Rejection 1: either a or b are outside a hub
                test1 = all(dist_a > richclublength/2) | all(dist_b > richclublength/2);
                % Rejection 2: a and b are both inside the same hub
                test2 = any(dist_a <= richclublength/2 & dist_b <= richclublength/2);
                if test1 || test2
                    test = 0;
                    a = topology.L * rand(2, 1);
                    b = topology.L * rand(2, 1);     
                end
            end
        end
        a_array1(:, j) = a; b_array1(:, j) = b;
    end

    ax = nexttile(20 + 3*i - 2, [1 2]);
    hold;
    set(gca, 'Color', 'white');
    ax.Box = "on";
    ax.LineWidth = 1;

    xlim([0 topology.L])
    ylim([0 topology.L])
    xticks([]);
    yticks([]);
    for k = 1:num_lrcs
        quiver(a_array1(1, k), a_array1(2, k), ...
        b_array1(1, k) - a_array1(1, k), ...
        b_array1(2, k) - a_array1(2, k),...
        'Color', 'k', 'LineWidth', 0.5, ...
        'MaxHeadSize', 0.05 / norm(b_array1(:, k) - a_array1(:, k)), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    end

    xlabel(append('$\lambda_c = ', num2str(lambda), '$'), 'Interpreter', 'latex', 'FontSize', 12);

    % Draw hub region
    for k = 1:num_richclubs
        fill(richclub_centres(k, 1) + richclublength/2 * [-1 -1 1 1], richclub_centres(k, 2) + richclublength/2 * [-1 1 1 -1], 'k', 'FaceColor', '#D95319', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    end
    
    if i == 2
        title('Core Specificity', 'Interpreter', 'latex', 'FontSize', 12);
    end
    
    if i == 1
        title('iii.', 'FontSize', 16, 'HorizontalAlignment', 'left');
        ax.TitleHorizontalAlignment = 'left';
    end
       
    hold off;

end



for i = 1:2
    nexttile(20 + 3*i);
    hold;
    quiver(0.35, 0.5, 0.3, 0, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 1, 'Marker', '.', 'MarkerSize', 0.0001, 'AutoScale', 'off');
    xlim([0 1])
    ylim([0 1])
    xticks([]);
    yticks([]);
    axis off;
end

%% Plot mean BOLD dissimilarity versus three topological control parameters

clear; clc;
loadparam;

f = figure;
f.Position = [100, 100, 1400, 400];
hold;
title(['Cosine Distance vs. Probability of Hub Incidence'])
set(gcf, 'Color', 'white')

% Load random connectome statistics
load('dissimboldensemble_multilrcs.mat', 'md_array', 'num_lrcs_array')
% Keep only 10, 20, 50 and 100 LRCs
md_array0 = md_array;
clear md_array;
md_array0 = reshape(md_array0, [1 size(md_array0)]);

lambda = 188;
reciprobability = 0.31;
hubprobability = 1 - sqrt(12 / 68); % 12 hub regions identified from 68 across both hemispheres in VDH 2011

t = tiledlayout(1, 28, 'TileSpacing', 'None', 'Padding', 'Compact');

ax = nexttile(10*1 - 9, [1 8]); hold;

load('dissimboldensemble_edr_old.mat');

% Append random connectome statistics
md_array = cat(1, md_array0, md_array);
lambda_array = [0 lambda_array];

mean_dissimilarity = squeeze(mean(md_array, 3));
std_dissimilarity = squeeze(std(md_array, 1, 3));

set(gca,'ColorOrderIndex',1);
for i = 1:length(num_lrcs_array)
    h = errorbar(lambda_array, mean_dissimilarity(:, i), std_dissimilarity(:, i), 'o', 'MarkerSize', 0.000001, 'LineWidth', 1);
    alpha = 0.2;   
    % Set transparency (undocumented)
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
end
set(gca,'ColorOrderIndex',1);
for i = 1:length(num_lrcs_array)
    scatter(lambda_array, mean_dissimilarity(:, i), 50, 'filled');
end
set(gca,'ColorOrderIndex',1);
for i = 1:length(num_lrcs_array)
    plot(lambda_array, mean_dissimilarity(:, i));
end

% % Plot double exponential of best fit
% set(gca,'ColorOrderIndex',1);
% xi = linspace(min(reciprobability_array), max(reciprobability_array), 100);
% for i = 1:length(num_lrcs_array)
%     ft = fit(reciprobability_array', mean_dissimilarity(:, i), 'exp2', 'Upper', [Inf, 0, Inf, 0], 'Lower', [-Inf, 0, -Inf, -Inf]);
%     a = ft.a;
%     b = ft.b;
%     c = ft.c;    
%     d = ft.d;
%     plot(xi, a*exp(b*xi) + c*exp(d*xi));
% end

% text(0, 1, '$\langle C_{\mathrm{BOLD}} \rangle_{\lambda}$', 'Interpreter', 'latex', 'VerticalAlignment', 'top', 'FontSize', 20, 'Units', 'normalized')


xlabel('$\lambda_l$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$C_{\mathrm{BOLD}}$', 'Interpreter', 'latex', 'FontSize', 20, 'Units', 'normalized')

xlim([-0.05*max(xlim) 1.05*max(xlim)])
ylim([0, Inf])
yticks([0:0.02:0.1])
xticks(lambda_array)

% title('i.', 'FontSize', 16, 'HorizontalAlignment', 'left');
% ax.TitleHorizontalAlignment = 'left';

hold off;

ax = nexttile(10*2 - 9, [1 8]); hold;

load('dissimboldensemble_hub.mat');


% Append random connectome statistics
md_array = cat(1, md_array0, md_array);
lambda_array = [0 lambda_array];

mean_dissimilarity = squeeze(mean(md_array, 3));
std_dissimilarity = squeeze(std(md_array, 1, 3));

set(gca,'ColorOrderIndex',1);
for i = 1:length(num_lrcs_array)
    h = errorbar(lambda_array, mean_dissimilarity(:, i), std_dissimilarity(:, i), 'o', 'MarkerSize', 0.000001, 'LineWidth', 1);
    alpha = 0.2;   
    % Set transparency (undocumented)
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
end
set(gca,'ColorOrderIndex',1);
for i = 1:length(num_lrcs_array)
    scatter(lambda_array, mean_dissimilarity(:, i), 50, 'filled');
end
set(gca,'ColorOrderIndex',1);
for i = 1:length(num_lrcs_array)
    plot(lambda_array, mean_dissimilarity(:, i));
end

xlabel('$\lambda_h$', 'Interpreter', 'latex', 'FontSize', 16);

% text(0, 1, '$\langle C_{\mathrm{BOLD}} \rangle_{p_r}$', 'Interpreter', 'latex', 'VerticalAlignment', 'top', 'FontSize', 20, 'Units', 'normalized')

xlim([-0.05*max(xlim) 1.05*max(xlim)])
ylim([0, Inf])
yticks([0:0.02:0.1])
xticks(lambda_array);

% title('ii.', 'FontSize', 16, 'HorizontalAlignment', 'left');
% ax.TitleHorizontalAlignment = 'left';

hold off;

ax = nexttile(10*3 - 9, [1 8]); hold;

load('dissimboldensemble_core.mat');

% Append random connectome statistics
md_array = cat(1, md_array0, md_array);
lambda_array = [0 lambda_array];

mean_dissimilarity = squeeze(mean(md_array, 3));
std_dissimilarity = squeeze(std(md_array, 1, 3));

set(gca,'ColorOrderIndex',1);
for i = 1:length(num_lrcs_array)
    h = errorbar(lambda_array, mean_dissimilarity(:, i), std_dissimilarity(:, i), 'o', 'MarkerSize', 0.000001, 'LineWidth', 1);
    alpha = 0.2;   
    % Set transparency (undocumented)
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
end
set(gca,'ColorOrderIndex',1);
for i = 1:length(num_lrcs_array)
    scatter(lambda_array, mean_dissimilarity(:, i), 50, 'filled');
end
set(gca,'ColorOrderIndex',1);
for i = 1:length(num_lrcs_array)
    plot(lambda_array, mean_dissimilarity(:, i));
end

xlabel('$\lambda_c$', 'Interpreter', 'latex', 'FontSize', 16);

% text(0, 1, '$\langle C_{\mathrm{BOLD}} \rangle_{\rho_h}$', 'Interpreter', 'latex', 'VerticalAlignment', 'top', 'FontSize', 20, 'Units', 'normalized')

xlim([-0.05*max(xlim) 1.05*max(xlim)])
ylim([0, Inf])
yticks([0:0.02:0.1])
xticks(lambda_array);

% title('iii.', 'FontSize', 16, 'HorizontalAlignment', 'left');
% ax.TitleHorizontalAlignment = 'left';

hold off;

Legend = cell(1, 2*length(num_lrcs_array));
for i = 1:length(num_lrcs_array)
    Legend{i} = "";
    Legend{length(num_lrcs_array) + i} = append('$', num2str(num_lrcs_array(i)), '$ tracts');
end

l = legend(Legend, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 12);
l.Layout.Tile = 'south';
l.Orientation = 'horizontal';

% save as dissimboldensemble_topologyconstraints.svg

%% Animate BOLD evoked response of 
% homogeneous model, and heterogeneous model with 10 and 50 recip LRCs

clear; clc;

loadparam;

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

% Reduce the intensity of the stimulus (seems to cause unwanted
% nonlinearities otherwise)

stim.stimI = 0.01;

timepoints = round((0:1:5) / topology.T);

num_snapshots = length(timepoints);

num_shortcuts_array = [0, 10, 50];
max_num_shortcuts = max(num_shortcuts_array);

region_radius = hetparam_uni.Ri(3, 1);
min_distance = 0.09;
max_distance = topology.L / 2;

% Create LRC topology
rng(1, "twister");
r_centre = zeros(2*max_num_shortcuts, 2);
for j = 1:max_num_shortcuts
    
    % Add ith shortcut

    dist = Inf;
    while (dist < min_distance || dist > max_distance)
        ri_centre = dx*randi(topology.Nx, 1, 2);
        rf_centre = dx*randi(topology.Nx, 1, 2);
        distx = abs(ri_centre(1) - rf_centre(1));
        if distx > topology.L / 2
            distx = topology.L - distx;
        end
        disty = abs(ri_centre(2) - rf_centre(2));
        if disty > topology.L / 2
            disty = topology.L - disty;
        end
        dist = norm([distx; disty], 2);
        unif_rv = rand();
        if unif_rv > (1/dist) / (1/min_distance)
            dist = Inf;
        end         
    end
    r_centre(2*j-1, :) = ri_centre;
    r_centre(2*j, :) = rf_centre;
end

hetparam_het_array = cell(1, length(num_shortcuts_array));
bold_full = zeros(topology.Nx, topology.Nx, num_snapshots, length(num_shortcuts_array));

for iter = 1:length(num_shortcuts_array)
    num_shortcuts = num_shortcuts_array(iter);
    if num_shortcuts == 0
        hetparam_het = hetparam_hom;
    else
        hetparam_het = hetparam_uni;
        hetparam_het.m = 2*num_shortcuts;
        hetparam_het.c = ones(1, 2*num_shortcuts);
        hetparam_het.tau = zeros(1, 2*num_shortcuts);
        hetparam_het.Ri = zeros(3, 2*num_shortcuts);
        hetparam_het.Rf = zeros(3, 2*num_shortcuts);
        hetparam_het.Ri(3, :) = region_radius;
        hetparam_het.Rf(3, :) = region_radius;
    
        for j2 = 1:num_shortcuts
            hetparam_het.Ri(1:2, j2) = r_centre(2*j2-1, :);
            hetparam_het.Rf(1:2, j2) = r_centre(2*j2, :);
            hetparam_het.Ri(1:2, j2 + num_shortcuts) = r_centre(2*j2, :);
            hetparam_het.Rf(1:2, j2 + num_shortcuts) = r_centre(2*j2-1, :);
        end
    end

    hetparam_het_array{iter} = hetparam_het;

    init = 0;
    s_init = 0; f_init = 1; v_init = 1; q_init = 1;

    %Skip the first snapshot (leave it zero) as its when t = 0
    snapshot_count = 2;

    for i = 1:max(timepoints)
        if (i > 1)
            stim.stimnum = 0;
        else
            stim.stimnum = 1;
        end
        ts = run_absorblayer(topology,homparam,hetparam_het,stim);
        init = ts(:, :, end-1:end);
        [bold, s_init, f_init, v_init, q_init] = ...
            balloonmodel(ts, dt, s_init, f_init, v_init, q_init);
        bold = bold(:, :, end);
        disp(append('iter ', num2str(i)));
        if i == timepoints(snapshot_count)
            bold_full(:, :, snapshot_count, iter) = bold;
            disp(append('added snapshot ', num2str(snapshot_count), ' to array'));
            snapshot_count = snapshot_count + 1;
        end
    end

end

% DATA SAVED AS BOLDresponse.mat

%%

clear; 
clc; 
loadparam;
load("BOLDresponse.mat");
dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

maxval = max(bold_full(:, :, :, :), [], 'all');
titles = ["Homogeneous", "10 LRCs", "50 LRCs"];

load 'CustomColormap.mat'
f = figure;
f.Position = [100, 100, 850, 450];
tiledlayout(length(num_shortcuts_array),num_snapshots,"TileSpacing","tight","Padding","tight");
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

for iter = 1:length(num_shortcuts_array)
    bold = squeeze(bold_full(:, :, :, iter));
    hetparam_het = hetparam_het_array{iter};
    for i = 1:num_snapshots
        ax = nexttile;
        hold on;
        ax.Box = "on";
        ax.LineWidth = 3;
        ax.XColor = 0.5*[1 1 1];
        ax.YColor = 0.5*[1 1 1];
        surf(X, Y, bold(:, :, i));
        view(0,90);
        colormap(ax, "bone");
        shading flat;
        clim([0 maxval]);
        xlim([0, topology.L]);
        ylim([0, topology.L]);
        xticks([]);
        yticks([]);
        if i == 1
            for m = 1:hetparam_het.m
                quiver3(dx/2 + hetparam_het.Ri(1, m), dx/2 + hetparam_het.Ri(2, m), maxval, ...
                hetparam_het.Rf(1, m) - hetparam_het.Ri(1, m), ...
                hetparam_het.Rf(2, m) - hetparam_het.Ri(2, m), 0,...
                'Color', 'g', 'LineWidth', 0.5, ...
                'MaxHeadSize', 0.0 / norm(hetparam_uni.Ri(:, 1) - hetparam_uni.Rf(:, 1)), ...
                'Marker', '.', 'MarkerSize', 5, ...
                'AutoScale','off');
            end
            plot3(stim.stimR(1), stim.stimR(2), maxval, 'x', 'MarkerSize', 15, 'Color', 'white', 'LineWidth', 2);
            ylabel(titles(iter), 'Color', 'k', 'FontSize', 15)
        end
        if iter == 1
            if i > 1
                title(append(num2str(i - 1), " s"), 'FontSize', 15);
            else
                title(append('t = ',num2str(i - 1), " s"), 'FontSize', 15);
            end
        end
        hold off;
    end
end

%% Plot divergence of random generated models over time

clear; clc;
loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

time_array = dt * (1:topology.Nt);

% Plot divergence of default model and random models on same plot

f = figure;
f.Position = [100, 100, 500, 300];
hold;

load('dissimovertime.mat');
dt = topology.T / topology.Nt;
time_array = dt*(1:topology.Nt);

mean_dissimilarity = zeros(topology.Nt, length(num_shortcuts_array));
for i = 1:length(num_shortcuts_array)
    mean_dissimilarity(:, i) = squeeze(mean(dissim_array(:, i, :), 3));
end

for i = 1:length(num_shortcuts_array)
    for j = 1:num_samples
        currentColor = get(gca, 'ColorOrder');
        currentColor = currentColor(i, :);
        mixedColor = 0.1*currentColor + 0.9*[1 1 1];
        plot(1000*(time_array - stim.stimt), dissim_array(:, i, j), 'Color', mixedColor);
    end
end

set(gca,'ColorOrderIndex',1);
for i = 1:length(num_shortcuts_array)
    plot(1000*(time_array - stim.stimt), mean_dissimilarity(:, i), 'LineWidth', 2);
end

Legend = cell(1, (num_samples+1)*length(num_shortcuts_array));
for i = 1:num_samples*length(num_shortcuts_array)
    Legend{i} = "";
end
for i = 1:length(num_shortcuts_array)
    Legend{num_samples*length(num_shortcuts_array) + i} = append(num2str(num_shortcuts_array(i)), ' LRCs');
end

ylim([-0.005, Inf]);
xlim([-10, Inf]);

xlabel('Time (ms)');
ylabel('Divergence')
legend(Legend);

% save as divergencedemo.svg

%% Animation homogeneous model under BALLOON FILTER transform

clear; clc;

loadparam;

load 'CustomColormap.mat'

topology.T = 0.01;
topology.Nx = 100;
dx = topology.L / topology.Nx;
dt = dx / (homparam.r * homparam.gamma * sqrt(2) * 2);
topology.Nt = ceil(topology.T / dt);
clear dx dt;

stim.stimI = 0.01;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

boldfull = [];
init_hom = 0;
init_het = 0;
s_hom_init = 0;
f_hom_init = 1;
v_hom_init = 1;
q_hom_init = 1;
s_het_init = 0;
f_het_init = 1;
v_het_init = 1;
q_het_init = 1;

for iter = 1:700
    ts_hom = run_init_periodic(topology,homparam,hetparam_hom,stim,init_hom);
    init_hom = ts_hom(:, :, end-1:end);
    [bold_hom, s_hom_init, f_hom_init, v_hom_init, q_hom_init] = ...
        balloonmodel(ts_hom, dt, s_hom_init, f_hom_init, v_hom_init, q_hom_init);
    if mod(iter, 50) == 0
        boldfull = cat(3, boldfull, bold_hom(:, :, end));
    end
    disp(iter);
    stim.stimnum = 0;
end

%%

maxval = max(abs(boldfull(:, :, :, :)), [], "all");

f = figure;
sgtitle('Observed fMRI Response')
f.Position = [10 100 425 400];
set(gcf, 'color', 'white');

subplot(1, 1, 1);
h = surf(X, Y, 0*boldfull(:, :, 1));
hold on;
hx = scatter3([dx/2 + stim.stimR(1)], [dx/2 + stim.stimR(2)], 0, 400, 'Marker', 'x', ...
    'MarkerEdgeColor', 'white', 'LineWidth', 2);
view(0,90);
colormap(CustomColormap);
shading flat;
clim([-maxval maxval]);
xlim([dx topology.L]);
ylim([dx topology.L]);
xticks([]);
yticks([]);
xlabel('t = -2.0 s');

gif('OriginalBOLD.gif', 'DelayTime', 1/2)
for count = -3:1:size(boldfull,3)
    if count > 0
        set(h, 'Zdata', boldfull(:,:,count));
        set(hx, 'Zdata', -maxval)
    end
    xlabel(append('t = ', sprintf('%.1f', count*0.5), ' s'));
    pause(0.5);
    gif;
end


%% Animate bidirectional model

clear; clc;

loadparam;

homparam.nu0 = 1;
hetparam_bi.Ri(1:2, 1) = [0.15; 0.2];
hetparam_bi.Rf(1:2, 1) = [0.25; 0.3];
hetparam_bi.Ri(1:2, 2) = hetparam_bi.Rf(1:2, 1);
hetparam_bi.Rf(1:2, 2) = hetparam_bi.Ri(1:2, 1);

% Increase simulation time
dt = topology.T / topology.Nt;
topology.T = 0.04;
topology.Nt = ceil(topology.T / dt);

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

titles = "1 shortcut";

ts = run_absorblayer(topology,homparam,hetparam_bi,stim);

maxval = 1*max(ts(:, :, end), [], 'all');
maxvalall = 1*max(ts, [], 'all');
%%
f = figure;
sgtitle('Response to Impulse Stimulation')
f.Position = [10 100 400 400];
set(gcf, 'color', 'white');

formatSpec = '%.0f';

h = surf(X, Y, 0*ts(:, :, 1));
hold on;
hx = scatter3([stim.stimR(1)], [stim.stimR(2)], 0, 400, 'Marker', 'x', ...
    'MarkerEdgeColor', 'white', 'LineWidth', 2);
plot3(dx/2 + [hetparam_bi.Ri(1, 1),hetparam_bi.Rf(1, 1)], ...
    dx/2 + [hetparam_bi.Ri(2, 1),hetparam_bi.Rf(2, 1)], ...
    [maxvalall, maxvalall], ...
    'Color', 'g', 'LineWidth', 1);
xlabel(append('t = ', num2str(1000*-200*dt, formatSpec), ' ms'));
%title(titles);
view(0,90);
colormap(CustomColormap);
shading flat;
clim([-maxval maxval]);
xlim([0 topology.L]);
ylim([0 topology.L]);
xticks([]);
yticks([]);

gif('HeterogeneousOriginalModified.gif','DelayTime',1/24)

for count = -200:5:topology.Nt
    if count > 0
        set(h, 'Zdata', ts(:,:,count))
        set(hx, 'Zdata', -maxvalall);
    end
    xlabel(append('t = ', num2str(1000*count*dt, formatSpec), ' ms'));
    pause(0);
    gif;
end

%% Animate different numbers of shortcuts

clear; clc;

loadparam;

homparam.nu0 = 1;

% Increase simulation time
dt = topology.T / topology.Nt;
topology.T = 0.04;
topology.Nt = ceil(topology.T / dt);

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

region_radius = hetparam_uni.Ri(3, 1);

num_shortcuts_array = [10, 50];

titles = ["0 shortcuts (Homogenous)" "10 shortcuts" "50 shortcuts"];

ts = zeros(topology.Nx, topology.Nx, topology.Nt, 3);
ts(:, :, :, 1) = run_absorblayer(topology,homparam,hetparam_hom,stim);
for i = 1:2
    num_pairs = num_shortcuts_array(i);
    hetparam_bi.m = 2*num_pairs;
    hetparam_bi.c = ones(1, 2*num_pairs);
    hetparam_bi.tau = zeros(1, 2*num_pairs);
    hetparam_bi.Ri = zeros(3, 2*num_pairs);
    hetparam_bi.Rf = zeros(3, 2*num_pairs);
    hetparam_bi.Ri(3, :) = region_radius;
    hetparam_bi.Rf(3, :) = region_radius;
        
    load(append('r_centre-',num2str(num_pairs), 'pairs.mat'));
    
    r_centre = dx * ceil(r_centre * topology.L / dx);
    
    for k = 1:num_pairs
        hetparam_bi.Ri(1:2, 2*k-1) = r_centre(2*k-1, :);
        hetparam_bi.Rf(1:2, 2*k-1) = r_centre(2*k, :);
        hetparam_bi.Ri(1:2, 2*k) = r_centre(2*k, :);
        hetparam_bi.Rf(1:2, 2*k) = r_centre(2*k-1, :);
    end
    hetparam_array(i) = hetparam_bi;
    ts(:, :, :, i+1) = run_absorblayer(topology,homparam,hetparam_bi,stim);
end
%%
maxval = 1.5*max(ts(:, :, end, :), [], 'all');
maxvalall = 1.5*max(ts(:, :, :, :), [], 'all');

f = figure;
sgtitle('Response to Impulse Stimulation')
f.Position = [10 100 1400 400];
set(gcf, 'color', 'white');

formatSpec = '%.0f';

for i = 1:3
    subplot(1, 3, i);
    h(i) = surf(X, Y, 0*ts(:, :, 1, 1));
    hold on;
    hx(i) = scatter3([stim.stimR(1)], [stim.stimR(2)], 0, 400, 'Marker', 'x', ...
        'MarkerEdgeColor', 'white', 'LineWidth', 2);
    if i > 1
        hetparam_bi = hetparam_array(i-1);
        for k = 1:num_shortcuts_array(i-1)
            plot3(dx/2 + [hetparam_bi.Ri(1, 2*k-1),hetparam_bi.Rf(1, 2*k-1)], ...
                dx/2 + [hetparam_bi.Ri(2, 2*k-1),hetparam_bi.Rf(2, 2*k-1)], ...
                [maxvalall, maxvalall], ...
                'Color', 'g', 'LineWidth', 1);
        end
    end
    title(titles(i));
    xlabel(append('t = ', num2str(1000*-200*dt, formatSpec), ' ms'));
    view(0,90);
    colormap(CustomColormap);
    shading flat;
    clim([-maxval maxval]);
    xlim([0 topology.L]);
    ylim([0 topology.L]);
    xticks([]);
    yticks([]);
end

gif('HeterogeneousLargeNumShortcutsBidirectional.gif','DelayTime',1/24)
for count = -200:5:topology.Nt
    for i = 1:3
        if count > 0
            set(h(i), 'Zdata', ts(:,:,count,i))
            set(hx(i), 'Zdata', -maxvalall);
        end
        xlabel(subplot(1, 3, i), append('t = ', num2str(1000*count*dt, formatSpec), ' ms'));
    end
    pause(0);
    gif
end

%% Animate unidirectional model
% Display with annotated INTERMEDIATE POINTS

clear; clc;

loadparam;

% Increase simulation time
dt = topology.T / topology.Nt;
topology.T = 0.04;
topology.Nt = ceil(topology.T / dt);

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

num_positions = 5;
regionicentre = hetparam_uni.Ri(1:2,1)';
regionfcentre = hetparam_uni.Rf(1:2,1)';
positions = linspace(regionicentre(1), regionfcentre(1), num_positions);
positions_grid = round(positions / dx);

ts = zeros(topology.Nx, topology.Nx, topology.Nt, 2);
for i = 1:2
    if i == 1
        ts_values = run_absorblayer(topology,homparam,hetparam_hom,stim);
    else
        ts_values = run_absorblayer(topology,homparam,hetparam_uni,stim);
    end
    ts(:, :, :, i) = ts_values;
end
%%
maxval = 3*max(ts(:, :, end, :), [], 'all');
maxvalall = 3*max(ts(:, :, :, :), [], 'all');

f = figure;
f.Position = [50 100 1000 400];
set(gcf, 'color', 'white');
titles = ["Homogeneous (0 shortcuts)" "1 shortcut"];

formatSpec = '%.0f';

for i = 1:2
    subplot(1, 2, i);
    hold on;
    xlabel(append('t = ', num2str(1000*-200*dt, formatSpec), ' ms'));
    h(i) = surf(X, Y, ts(:, :, 1, i));
    scatter3([dx/2 + stim.stimR(1)], [dx/2 + stim.stimR(2)], 0, 400, 'Marker', 'x', ...
        'MarkerEdgeColor', 'white', 'LineWidth', 2);
    scatter3(positions, positions, ...
        maxvalall, 600, 'Marker', '.', 'MarkerFaceColor', 'white');
    if i == 2
        for k = 1:hetparam_uni.m
            quiver3(dx/2 + hetparam_uni.Ri(1, k), dx/2 + hetparam_uni.Ri(2, k), maxvalall, ...
                hetparam_uni.Rf(1, k) - hetparam_uni.Ri(1, k), ...
                hetparam_uni.Rf(2, k) - hetparam_uni.Ri(2, k), 0, ...
                'Color', 'g', 'LineWidth', 2, ...
                'MaxHeadSize', 0.1 / norm(hetparam_uni.Ri(:, k) - hetparam_uni.Rf(:, k)), ...
                'Marker', '.', 'MarkerSize', 5, ...
                'AutoScale','off');
        end
    end
    title(titles(i));
    view(0,90);
    colormap(CustomColormap);
    shading flat;
    clim([-maxval maxval]);
    xlim([0 topology.L]);
    ylim([0 topology.L]);
    xticks([]);
    yticks([]);
end

gif('OriginalModelIntermediate.gif','DelayTime',1/24)

for count = -200:5:topology.Nt
    for i = 1:2
        xlabel(subplot(1, 2, i), append('t = ', num2str(1000*count*dt, formatSpec), ' ms'));
        set(h(i), 'Zdata', ts(:,:,max(count,1),i))
    end
    pause(0);
    gif
end

%% Plot response of homogeneous model only

Legend = cell(2,1);
f = figure;
f.Position = [0 100 1500 350];
tiledlayout('flow', 'TileSpacing', 'None','Padding', 'None');
set(gcf, 'color', 'white');

colorOrder = get(gca, 'ColorOrder');
for iter = 1:num_positions
    nexttile; 
    hold;
    title(['at Position ', num2str(iter)], 'Color', colorOrder(iter+1, :)); 
    xlabel('t (ms)'); 
    if iter == 1
        ylabel('Activity');
    end
    g = plot(1000*dt*[1:topology.Nt], reshape(ts(positions_grid(iter), positions_grid(iter), :, 1), [], 1), 'Color', 'Red', 'LineWidth', 2, 'LineStyle', '--');
    h(iter) = scatter(1000*dt*1, ts(positions_grid(iter), positions_grid(iter), 1, 1), 'r', 'filled');
    xlim(1000*[0, 0.03])
    ylim([0, max(ts(positions_grid(iter), positions_grid(iter), :, :), [], 'all')]);
    hold off;
end

sgtitle('Temporal Response');
hl = legend({"Homogeneneous"});
hl.Layout.Tile = 'south';

% save as temprespRiRf.svg

gif('TemporalResponseHomogeneousModel.gif','DelayTime',1/24)

for count = -200:5:topology.Nt
    if count > 0
        for iter = 1:num_positions
            set(h(iter), 'Xdata', 1000*dt*[count])
            set(h(iter), 'Ydata', ts(positions_grid(iter), positions_grid(iter),count,1))
        end
    end
    pause(0);
    gif
end

%% Plot response of heterogeneous model on top

Legend = cell(2,1);
f = figure;
f.Position = [0 100 1500 350];
tiledlayout('flow', 'TileSpacing', 'None','Padding', 'None');
set(gcf, 'color', 'white');

colorOrder = get(gca, 'ColorOrder');
for iter = 1:num_positions
    nexttile(iter); 
    hold;
    title(['at Position ', num2str(iter)], 'Color', colorOrder(iter+1, :)); 
    xlabel('t (ms)'); 
    if iter == 1
        ylabel('Activity');
    end
    plot(1000*dt*[1:topology.Nt], reshape(ts(positions_grid(iter), positions_grid(iter), :, 1), [], 1), 'Color', 'Red', 'LineWidth', 2, 'LineStyle', '--');
    plot(1000*dt*[1:topology.Nt], reshape(ts(positions_grid(iter), positions_grid(iter), :, 2), [], 1), 'Color', 'Red', 'LineWidth', 1);
    h(iter) = scatter(1000*dt*1, ts(positions_grid(iter), positions_grid(iter), 1, 2), 'r', 'filled');
    xlim(1000*[0, 0.03])
    ylim([0, max(ts(positions_grid(iter), positions_grid(iter), :, :), [], 'all')]);
    hold off;
end
Legend{1} = "homogeneous";
Legend{2} = "1 shortcut";
Legend{3} = "";

hl = legend(Legend);
hl.Layout.Tile = 'south';
sgtitle('Temporal Response');

% save as temprespRiRf.svg

gif('TemporalResponseOriginalModel.gif','DelayTime',1/24)

for count = -200:5:topology.Nt
    if count > 0
        for iter = 1:num_positions
            for i = 2
                set(h(iter), 'Xdata', 1000*dt*[count])
                set(h(iter), 'Ydata', ts(positions_grid(iter), positions_grid(iter),count,i))
            end
        end
    end
    pause(0);
    gif
end

%% Plot dissimilarity of original heterogeneous model over time

clear; clc;

loadparam;

% Increase simulation time
dt = topology.T / topology.Nt;
topology.T = 0.05;
topology.Nt = ceil(topology.T / dt);

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

num_positions = 5;
regionicentre = hetparam_uni.Ri(1:2,1)';
regionfcentre = hetparam_uni.Rf(1:2,1)';
positions = linspace(regionicentre(1), regionfcentre(1), num_positions);
positions_grid = round(positions / dx);

ts = zeros(topology.Nx, topology.Nx, topology.Nt, 2);
for i = 1:2
    if i == 1
        ts_values = run_absorblayer(topology,homparam,hetparam_hom,stim);
    else
        ts_values = run_absorblayer(topology,homparam,hetparam_uni,stim);
    end
    ts(:, :, :, i) = ts_values;
end

norm1 = sqrt(sum((ts(:, :, :, 1)).^2, [1 2]));
norm2 = sqrt(sum((ts(:, :, :, 2)).^2, [1 2]));
dissim = 1 - (sum(ts(:, :, :, 2) .* ts(:, :, :, 1), [1 2])) ./ (norm1 .* norm2);
dissim = squeeze(dissim);
%%
maxval = 1.5*max(ts(:, :, end, :), [], 'all');
maxvalall = 1.5*max(ts(:, :, :, :), [], 'all');

f = figure;
sgtitle('Dissimilarity over Time')
f.Position = [10 100 1400 400];
titles = ["Homogeneous (0 shortcuts)", "1 shortcut", "Dissimilarity"];
set(gcf, 'Color', 'White')

formatSpec = '%.0f';

for i = 1:2
    subplot(1, 3, i);
    hold on;
    xlabel(append('t = ', num2str(1000*-200*dt, formatSpec), ' ms'));
    h(i) = surf(X, Y, 0*ts(:, :, 1, i));
    scatter3([dx/2 + stim.stimR(1)], [dx/2 + stim.stimR(2)], 0, 400, 'Marker', 'x', ...
        'MarkerEdgeColor', 'white', 'LineWidth', 2);
    if i == 2
        for k = 1:hetparam_uni.m
            quiver3(dx/2 + hetparam_uni.Ri(1, k), dx/2 + hetparam_uni.Ri(2, k), maxvalall, ...
                hetparam_uni.Rf(1, k) - hetparam_uni.Ri(1, k), ...
                hetparam_uni.Rf(2, k) - hetparam_uni.Ri(2, k), 0, ...
                'Color', 'g', 'LineWidth', 2, ...
                'MaxHeadSize', 0.1 / norm(hetparam_uni.Ri(:, k) - hetparam_uni.Rf(:, k)), ...
                'Marker', '.', 'MarkerSize', 5, ...
                'AutoScale','off');
        end
    end
    title(titles(i));
    view(0,90);
    colormap(CustomColormap);
    shading flat;
    clim([-maxval maxval]);
    xlim([0 topology.L]);
    ylim([0 topology.L]);
    xticks([]);
    yticks([]);
    hold off;
end

subplot(1, 3, 3);
hold;
h(3) = plot(1000*dt*[1:1], dissim(1:1));
title(titles(3));
xlabel('t (ms)');
ylabel('D(t)');
xlim([0, 1000*topology.T]);
ylim([0, 0.27]);
%xticks([]);
%yticks([]);

gifname = append('DissimilarityoverTimeOriginalModel.gif');
gif(gifname,'DelayTime',1/24)

for count = -200:5:topology.Nt
    for i = 1:2
        xlabel(subplot(1, 3, i), append('t = ', num2str(1000*count*dt, formatSpec), ' ms'));
    end
    if count > 0
        for i = 1:2
            set(h(i), 'Zdata', ts(:,:,count,i))
        end
        set(h(3), 'XData', 1000*dt*(1:count));
        set(h(3), 'YData', dissim(1:count));
    end
    pause(0);
    gif
end

%% Plot dissimilarity of 50 shortcut heterogeneous model over time

clear; clc;
loadparam;

% Increase simulation time
dt = topology.T / topology.Nt;
topology.T = 0.05;
topology.Nt = ceil(topology.T / dt);

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

num_shortcuts = 50;
hetparam_uni.m = num_shortcuts;
hetparam_uni.c = ones(1, num_shortcuts);
hetparam_uni.tau = zeros(1, num_shortcuts);
region_radius = hetparam_uni.Ri(3, 1);
hetparam_uni.Ri = zeros(3, num_shortcuts);
hetparam_uni.Rf = zeros(3, num_shortcuts);
hetparam_uni.Ri(3, :) = region_radius;
hetparam_uni.Rf(3, :) = region_radius;
    
load(append('r_centre-',num2str(num_shortcuts), 'pairs.mat'));

r_centre = dx * ceil(r_centre * topology.L / dx);

for k = 1:num_shortcuts
    hetparam_uni.Ri(1:2, k) = r_centre(2*k-1, :);
    hetparam_uni.Rf(1:2, k) = r_centre(2*k, :);
end

ts = zeros(topology.Nx, topology.Nx, topology.Nt, 2);
for i = 1:2
    if i == 1
        ts_values = run_absorblayer(topology,homparam,hetparam_hom,stim);
    else
        ts_values = run_absorblayer(topology,homparam,hetparam_uni,stim);
    end
    ts(:, :, :, i) = ts_values;
end


norm1 = sqrt(sum((ts(:, :, :, 1)).^2, [1 2]));
norm2 = sqrt(sum((ts(:, :, :, 2)).^2, [1 2]));
dissim = 1 - (sum(ts(:, :, :, 2) .* ts(:, :, :, 1), [1 2])) ./ (norm1 .* norm2);
dissim = squeeze(dissim);
max_dissim = max(dissim);
maxval = 1.5*max(ts(:, :, end, :), [], 'all');
maxvalall = 1.5*max(ts(:, :, :, :), [], 'all');
%%
f = figure;
sgtitle('Dissimilarity over Time')
f.Position = [10 100 1400 400];
titles = ["Homogeneous (0 shortcuts)" append(num2str(num_shortcuts), " shortcuts") "Dissimilarity"];
set(gcf, 'Color', 'White')

formatSpec = '%.0f';

for i = 1:2
    subplot(1, 3, i);
    hold on;
    xlabel(append('t = ', num2str(1000*-200*dt, formatSpec), ' ms'));
    h(i) = surf(X, Y, 0*ts(:, :, 1, i));
    scatter3([dx/2 + stim.stimR(1)], [dx/2 + stim.stimR(2)], 0, 400, 'Marker', 'x', ...
        'MarkerEdgeColor', 'white', 'LineWidth', 2);
    if i == 2
        for k = 1:hetparam_uni.m
            quiver3(dx/2 + hetparam_uni.Ri(1, k), dx/2 + hetparam_uni.Ri(2, k), maxvalall, ...
                hetparam_uni.Rf(1, k) - hetparam_uni.Ri(1, k), ...
                hetparam_uni.Rf(2, k) - hetparam_uni.Ri(2, k), 0, ...
                'Color', 'g', 'LineWidth', 0.5, ...
                'MaxHeadSize', 0.1 / norm(hetparam_uni.Ri(:, k) - hetparam_uni.Rf(:, k)), ...
                'Marker', '.', 'MarkerSize', 5, ...
                'AutoScale','off');
        end
    end
    title(titles(i));
    view(0,90);
    colormap(CustomColormap);
    shading flat;
    clim([-maxval maxval]);
    xlim([0 topology.L]);
    ylim([0 topology.L]);
    xticks([]);
    yticks([]);
    hold off;
end

subplot(1, 3, 3);
hold;
h(3) = plot(dt*[1:1], dissim(1:1));
title(titles(3));
xlabel('t (ms)');
ylabel('D(t)');
xlim([0, 1000*topology.T]);
ylim([0, 0.27]);
%xticks([]);
%yticks([]);

gifname = append('DissimilarityoverTime', num2str(num_shortcuts), 'shortcuts.gif');
gif(gifname,'DelayTime',1/24)

for count = -200:5:topology.Nt
    for i = 1:2
        xlabel(subplot(1, 3, i), append('t = ', num2str(1000*count*dt, formatSpec), ' ms'));
    end
    if count > 0
        for i = 1:2
            set(h(i), 'Zdata', ts(:,:,count,i))
        end
        set(h(3), 'XData', 1000*dt*(1:count));
        set(h(3), 'YData', dissim(1:count));
    end
    pause(0);
    gif
end

%% Plot dissimilarity of 10 bidirectional shortcut heterogeneous model over time

clear; clc;
loadparam;

homparam.nu0 = 1;

% Increase simulation time
dt = topology.T / topology.Nt;
topology.T = 0.05;
topology.Nt = ceil(topology.T / dt);

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

num_shortcuts = 10;
hetparam_bi.m = 2*num_shortcuts;
hetparam_bi.c = ones(1, 2*num_shortcuts);
hetparam_bi.tau = zeros(1, 2*num_shortcuts);
region_radius = hetparam_uni.Ri(3, 1);
hetparam_bi.Ri = zeros(3, 2*num_shortcuts);
hetparam_bi.Rf = zeros(3, 2*num_shortcuts);
hetparam_bi.Ri(3, :) = region_radius;
hetparam_bi.Rf(3, :) = region_radius;
    
load(append('r_centre-',num2str(num_shortcuts), 'pairs.mat'));

r_centre = dx * ceil(r_centre * topology.L / dx);

for k = 1:num_shortcuts
    hetparam_bi.Ri(1:2, 2*k-1) = r_centre(2*k-1, :);
    hetparam_bi.Rf(1:2, 2*k-1) = r_centre(2*k, :);
    hetparam_bi.Ri(1:2, 2*k) = r_centre(2*k, :);
    hetparam_bi.Rf(1:2, 2*k) = r_centre(2*k-1, :);
end

ts = zeros(topology.Nx, topology.Nx, topology.Nt, 2);
for i = 1:2
    if i == 1
        ts_values = run_absorblayer(topology,homparam,hetparam_hom,stim);
    else
        ts_values = run_absorblayer(topology,homparam,hetparam_bi,stim);
    end
    ts(:, :, :, i) = ts_values;
end


norm1 = sqrt(sum((ts(:, :, :, 1)).^2, [1 2]));
norm2 = sqrt(sum((ts(:, :, :, 2)).^2, [1 2]));
dissim = 1 - (sum(ts(:, :, :, 2) .* ts(:, :, :, 1), [1 2])) ./ (norm1 .* norm2);
dissim = squeeze(dissim);
%%
max_dissim = max(dissim);
maxval = 1.5*max(ts(:, :, end, :), [], 'all');
maxvalall = 1.5*max(ts(:, :, :, :), [], 'all');

f = figure;
sgtitle('Dissimilarity over Time')
f.Position = [10 100 1400 400];
titles = ["0 shortcuts (Homogeneous)" append(num2str(num_shortcuts), "  shortcuts") "Dissimilarity between Models"];
set(gcf, 'Color', 'White')

formatSpec = '%.0f';

for i = 1:2
    subplot(1, 3, i);
    hold on;
    xlabel(append('t = ', num2str(1000*-150*dt, formatSpec), ' ms'));
    h(i) = surf(X, Y, 0*ts(:, :, 1, i));
    hx(i) = scatter3([dx/2 + stim.stimR(1)], [dx/2 + stim.stimR(2)], 0, 400, 'Marker', 'x', ...
        'MarkerEdgeColor', 'white', 'LineWidth', 2);
    if i == 2
        for k = 1:num_shortcuts
            plot3(dx/2 + [hetparam_bi.Ri(1, 2*k-1),hetparam_bi.Rf(1, 2*k-1)], ...
                dx/2 + [hetparam_bi.Ri(2, 2*k-1),hetparam_bi.Rf(2, 2*k-1)], ...
                [maxvalall, maxvalall], ...
                'Color', 'g', 'LineWidth', 1);
        end
    end
    title(titles(i));
    view(0,90);
    colormap(CustomColormap);
    shading flat;
    clim([-maxval maxval]);
    xlim([0 topology.L]);
    ylim([0 topology.L]);
    xticks([]);
    yticks([]);
    hold off;
end

subplot(1, 3, 3);
hold;
h(3) = plot(1000*dt*[1:1], dissim(1:1));
title(titles(3));
xlabel('t (ms)');
ylabel('Dissimilarity');
xlim([0, 1000*topology.T]);
ylim([0, 0.27]);
%xticks([]);
%yticks([]);

gifname = append('DissimilarityoverTime', num2str(num_shortcuts), 'reciprocalshortcuts.gif');
gif(gifname,'DelayTime',1/24)

for count = -150:3:topology.Nt
    for i = 1:2
        xlabel(subplot(1, 3, i), append('t = ', num2str(1000*count*dt, formatSpec), ' ms'));
    end
    if count > 0
        for i = 1:2
            set(h(i), 'Zdata', ts(:,:,count,i))
            set(hx(i), 'Zdata', -maxvalall)
        end
        set(h(3), 'XData', 1000*dt*(1:count));
        set(h(3), 'YData', dissim(1:count));
    end
    pause(0);
    gif
end




%% Plot dissimilarity of 50 bidirectional shortcut heterogeneous model over time

clear; clc;
loadparam;

homparam.nu0 = 1;

% Increase simulation time
dt = topology.T / topology.Nt;
topology.T = 0.05;
topology.Nt = ceil(topology.T / dt);

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

num_shortcuts = 50;
hetparam_bi.m = 2*num_shortcuts;
hetparam_bi.c = ones(1, 2*num_shortcuts);
hetparam_bi.tau = zeros(1, 2*num_shortcuts);
region_radius = hetparam_uni.Ri(3, 1);
hetparam_bi.Ri = zeros(3, 2*num_shortcuts);
hetparam_bi.Rf = zeros(3, 2*num_shortcuts);
hetparam_bi.Ri(3, :) = region_radius;
hetparam_bi.Rf(3, :) = region_radius;
    
load(append('r_centre-',num2str(num_shortcuts), 'pairs.mat'));

r_centre = dx * ceil(r_centre * topology.L / dx);

for k = 1:num_shortcuts
    hetparam_bi.Ri(1:2, 2*k-1) = r_centre(2*k-1, :);
    hetparam_bi.Rf(1:2, 2*k-1) = r_centre(2*k, :);
    hetparam_bi.Ri(1:2, 2*k) = r_centre(2*k, :);
    hetparam_bi.Rf(1:2, 2*k) = r_centre(2*k-1, :);
end

ts = zeros(topology.Nx, topology.Nx, topology.Nt, 2);
for i = 1:2
    if i == 1
        ts_values = run_absorblayer(topology,homparam,hetparam_hom,stim);
    else
        ts_values = run_absorblayer(topology,homparam,hetparam_bi,stim);
    end
    ts(:, :, :, i) = ts_values;
end


norm1 = sqrt(sum((ts(:, :, :, 1)).^2, [1 2]));
norm2 = sqrt(sum((ts(:, :, :, 2)).^2, [1 2]));
dissim = 1 - (sum(ts(:, :, :, 2) .* ts(:, :, :, 1), [1 2])) ./ (norm1 .* norm2);
dissim = squeeze(dissim);
%%
max_dissim = max(dissim);
maxval = 1.5*max(ts(:, :, end, :), [], 'all');
maxvalall = 1.5*max(ts(:, :, :, :), [], 'all');

f = figure;
sgtitle('Dissimilarity over Time')
f.Position = [10 100 1400 400];
titles = ["0 shortcuts (Homogeneous)" append(num2str(num_shortcuts), "  shortcuts") "Dissimilarity between Models"];
set(gcf, 'Color', 'White')

formatSpec = '%.0f';

for i = 1:2
    subplot(1, 3, i);
    hold on;
    xlabel(append('t = ', num2str(1000*-150*dt, formatSpec), ' ms'));
    h(i) = surf(X, Y, 0*ts(:, :, 1, i));
    hx(i) = scatter3([dx/2 + stim.stimR(1)], [dx/2 + stim.stimR(2)], 0, 400, 'Marker', 'x', ...
        'MarkerEdgeColor', 'white', 'LineWidth', 2);
    if i == 2
        for k = 1:num_shortcuts
            plot3(dx/2 + [hetparam_bi.Ri(1, 2*k-1),hetparam_bi.Rf(1, 2*k-1)], ...
                dx/2 + [hetparam_bi.Ri(2, 2*k-1),hetparam_bi.Rf(2, 2*k-1)], ...
                [maxvalall, maxvalall], ...
                'Color', 'g', 'LineWidth', 1);
        end
    end
    title(titles(i));
    view(0,90);
    colormap(CustomColormap);
    shading flat;
    clim([-maxval maxval]);
    xlim([0 topology.L]);
    ylim([0 topology.L]);
    xticks([]);
    yticks([]);
    hold off;
end

subplot(1, 3, 3);
hold;
h(3) = plot(1000*dt*[1:1], dissim(1:1));
title(titles(3));
xlabel('t (ms)');
ylabel('Dissimilarity');
xlim([0, 1000*topology.T]);
ylim([0, 0.27]);
%xticks([]);
%yticks([]);

gifname = append('DissimilarityoverTime', num2str(num_shortcuts), 'reciprocalshortcuts.gif');
gif(gifname,'DelayTime',1/24)

for count = -150:3:topology.Nt
    for i = 1:2
        xlabel(subplot(1, 3, i), append('t = ', num2str(1000*count*dt, formatSpec), ' ms'));
    end
    if count > 0
        for i = 1:2
            set(h(i), 'Zdata', ts(:,:,count,i))
            set(hx(i), 'Zdata', -maxvalall)
        end
        set(h(3), 'XData', 1000*dt*(1:count));
        set(h(3), 'YData', dissim(1:count));
    end
    pause(0);
    gif
end





%% Demonstrate dissimilarity measure with unidirectional model

clear; clc;

loadparam;

load 'CustomColormap.mat'

topology.T = 0.07;
topology.Nx = 300;
dx = topology.L / topology.Nx;
dt = dx / (homparam.r * homparam.gamma * sqrt(2) * 2);
topology.Nt = ceil(topology.T / dt);
clear dx dt;

stim1 = stim;
stim1.stimR = [0.1; 0.15];

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

titles = ["Homogeneous" append("Heterogeneous - ",num2str(hetparam_uni.m), " shortcut"), "Difference"];

ts = zeros(topology.Nx, topology.Nx, topology.Nt, 3, 2);
for i = 1:2
    if i == 1
        ts_values = run_periodic(topology,homparam,hetparam_hom,stim);
    else
        ts_values = run_periodic(topology,homparam,hetparam_uni,stim);
    end
    ts(:, :, :, i, 1) = ts_values;
end

for i = 1:2
    if i == 1
        ts_values = run_periodic(topology,homparam,hetparam_hom,stim1);
    else
        ts_values = run_periodic(topology,homparam,hetparam_uni,stim1);
    end
    ts(:, :, :, i, 2) = ts_values;
end


u_cum = zeros(topology.Nx, topology.Nx, 2);
v_cum = zeros(topology.Nx, topology.Nx, 2);

for i = 1:topology.Nt
    u_cum = u_cum + squeeze(ts(:, :, i, 1, :));
    v_cum = v_cum + squeeze(ts(:, :, i, 2, :));
    ts(:, :, i, 1, :) = u_cum;
    ts(:, :, i, 2, :) = v_cum;
    ts(:, :, i, 3, :) = v_cum - u_cum;
end

%%
maxval = max(ts(:, :, end, 1:2, :), [], 'all');
maxdiff = max(ts(:, :, end, 3, :), [], 'all');

f = figure;
sgtitle('Time-Integrated Response to Impulse Stimulation')
set(gcf, 'color', 'white');
f.Position = [10 10 1250 1100];

for j = 1:2
    for i = 1:3
        subplot(2, 3, (j-1)*3 + i);
        h(i, j) = surf(X, Y, ts(:, :, 1, i, j));
        hold on;
        if j == 1
            scatter3([stim.stimR(1)], [stim.stimR(2)], max(ts, [], 'all'), 100, 'Marker', 'x');
        else
            scatter3([stim1.stimR(1)], [stim1.stimR(2)], max(ts, [], 'all'), 100, 'Marker', 'x');
        end
        if i == 2
            for k = 1:hetparam_uni.m
                quiver3(dx/2 + hetparam_uni.Ri(1, k), dx/2 + hetparam_uni.Ri(2, k), max(ts, [], 'all'), ...
                    hetparam_uni.Rf(1, k) - hetparam_uni.Ri(1, k), ...
                    hetparam_uni.Rf(2, k) - hetparam_uni.Ri(2, k), 0, ...
                    'Color', 'g', 'LineWidth', 0.5, ...
                    'MaxHeadSize', 0.05 / norm(hetparam_uni.Ri(:, k) - hetparam_uni.Rf(:, k)), ...
                    'Marker', '.', 'MarkerSize', 5, ...
                    'AutoScale','off');
            end
        end
        if (j == 1)
            title(titles(i));
        end
        if (i == 1)
            if j == 1
                ylabel("Stimulate on shortcut");
            else
                ylabel("Stimulate off shortcut");
            end
        end
        view(0,90);
        colormap(CustomColormap);
        shading flat;
        if i <= 2
            clim([-maxval maxval]);
        else
            clim([-maxdiff maxdiff]);
        end
        xlim([0 topology.L]);
        ylim([0 topology.L]);
        xticks([]);
        yticks([]);
    end
end

gif('HeterogeneousDissimilarity.gif','DelayTime',1/24)
for count = 1:10:topology.Nt
    for j = 1:2
        for i = 1:3
            set(h(i, j), 'Zdata', ts(:,:,count,i, j))
        end
    end
    pause(0);
    gif
end

%% Animation dissimilarity of original heterogenous model under 
% AFTER A REAL BALLOON FILTER transform

clear; clc;

loadparam;

load 'CustomColormap.mat'

topology.T = 0.01;
topology.Nx = 100;
dx = topology.L / topology.Nx;
dt = dx / (homparam.r * homparam.gamma * sqrt(2) * 2);
topology.Nt = ceil(topology.T / dt);
clear dx dt;

stim.stimI = 0.01;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

boldfull = [];
init_hom = 0;
init_het = 0;
s_hom_init = 0;
f_hom_init = 1;
v_hom_init = 1;
q_hom_init = 1;
s_het_init = 0;
f_het_init = 1;
v_het_init = 1;
q_het_init = 1;

for iter = 1:700
    ts_hom = run_init_periodic(topology,homparam,hetparam_hom,stim,init_hom);
    init_hom = ts_hom(:, :, end-1:end);
    ts_het = run_init_periodic(topology,homparam,hetparam_uni,stim,init_het);
    init_het = ts_het(:, :, end-1:end);
    [bold_hom, s_hom_init, f_hom_init, v_hom_init, q_hom_init] = ...
        balloonmodel(ts_hom, dt, s_hom_init, f_hom_init, v_hom_init, q_hom_init);
    [bold_het, s_het_init, f_het_init, v_het_init, q_het_init] = ...
        balloonmodel(ts_het, dt, s_het_init, f_het_init, v_het_init, q_het_init);
    if mod(iter, 50) == 0
        boldfull = cat(3, boldfull, ...
            cat(4, bold_het(:, :, end), bold_hom(:, :, end)));
    end
    disp(iter);
    stim.stimnum = 0;
end

%%

maxval = max(abs(boldfull(:, :, :, :)), [], "all");

titles = [append("Heterogeneous - 1 shortcut"), "Homogeneous Approximation"];

f = figure;
sgtitle('BOLD Response to Impulse Stimulation')
f.Position = [10 100 1400 600];
set(gcf, 'color', 'white');

for i = 1:2
    subplot(1, 2, i);
    h(i) = surf(X, Y, 0*boldfull(:, :, 1, i));
    hold on;
    scatter3([dx/2 + stim.stimR(1)], [dx/2 + stim.stimR(2)], 0, 400, 'Marker', 'x', ...
        'MarkerEdgeColor', 'white', 'LineWidth', 2);
    if i == 1
        for k = 1:hetparam_uni.m
            quiver3(dx/2 + hetparam_uni.Ri(1, k), dx/2 + hetparam_uni.Ri(2, k), 0, ...
                hetparam_uni.Rf(1, k) - hetparam_uni.Ri(1, k), ...
                hetparam_uni.Rf(2, k) - hetparam_uni.Ri(2, k), 0, ...
                'Color', 'g', 'LineWidth', 0.5, ...
                'MaxHeadSize', 0.05 / norm(hetparam_uni.Ri(:, k) - hetparam_uni.Rf(:, k)), ...
                'Marker', '.', 'MarkerSize', 5, ...
                'AutoScale','off');
        end
    end
    title(titles(i));
    view(0,90);
    colormap(CustomColormap);
    shading flat;
    clim([-maxval maxval]);
    xlim([0 topology.L]);
    ylim([0 topology.L]);
    xticks([]);
    yticks([]);
    xlabel('t = -2.0');
end

gif('OriginalBOLD.gif', 'DelayTime', 1/2)
for count = -3:1:size(boldfull,3)
    for i = 1:2
        if count > 0
            set(h(i), 'Zdata', boldfull(:,:,max(count,1),i));
        end
        xlabel(subplot(1, 2, i), append('t = ', sprintf('%.1f', count*0.5)));
    end
    pause(0.5);
    gif;
end


%% Animation dissimilarity of 50 shortcut model under 
% AFTER A REAL BALLOON FILTER transform

clear; clc;

loadparam;

load 'CustomColormap.mat'

topology.T = 0.01;
topology.Nx = 100;
dx = topology.L / topology.Nx;
dt = dx / (homparam.r * homparam.gamma * sqrt(2) * 2);
topology.Nt = ceil(topology.T / dt);
clear dx dt;

stim.stimI = 0.01;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

num_shortcuts = 50;
hetparam_uni.m = num_shortcuts;
hetparam_uni.c = ones(1, num_shortcuts);
hetparam_uni.tau = zeros(1, num_shortcuts);
region_radius = hetparam_uni.Ri(3, 1);
hetparam_uni.Ri = zeros(3, num_shortcuts);
hetparam_uni.Rf = zeros(3, num_shortcuts);
hetparam_uni.Ri(3, :) = region_radius;
hetparam_uni.Rf(3, :) = region_radius;
    
load(append('r_centre-',num2str(num_shortcuts), 'pairs.mat'));

r_centre = dx * ceil(r_centre * topology.L / dx);

for k = 1:num_shortcuts
    hetparam_uni.Ri(1:2, k) = r_centre(2*k-1, :);
    hetparam_uni.Rf(1:2, k) = r_centre(2*k, :);
end

boldfull = [];
init_hom = 0;
init_het = 0;
s_hom_init = 0;
f_hom_init = 1;
v_hom_init = 1;
q_hom_init = 1;
s_het_init = 0;
f_het_init = 1;
v_het_init = 1;
q_het_init = 1;

for iter = 1:700
    ts_hom = run_init_periodic(topology,homparam,hetparam_hom,stim,init_hom);
    init_hom = ts_hom(:, :, end-1:end);
    ts_het = run_init_periodic(topology,homparam,hetparam_uni,stim,init_het);
    init_het = ts_het(:, :, end-1:end);
    [bold_hom, s_hom_init, f_hom_init, v_hom_init, q_hom_init] = ...
        balloonmodel(ts_hom, dt, s_hom_init, f_hom_init, v_hom_init, q_hom_init);
    [bold_het, s_het_init, f_het_init, v_het_init, q_het_init] = ...
        balloonmodel(ts_het, dt, s_het_init, f_het_init, v_het_init, q_het_init);
    if mod(iter, 50) == 0
        boldfull = cat(3, boldfull, ...
            cat(4, bold_het(:, :, end), bold_hom(:, :, end)));
    end
    disp(iter);
    stim.stimnum = 0;
end

%%

maxval = max(abs(boldfull(:, :, :, :)), [], "all");

titles = [append("Heterogeneous - 50 shortcuts"), "Homogeneous Approximation"];

f = figure;
sgtitle('BOLD Response to Impulse Stimulation')
f.Position = [10 100 1400 600];
set(gcf, 'color', 'white');

for i = 1:2
    subplot(1, 2, i);
    h(i) = surf(X, Y, boldfull(:, :, 1, i));
    hold on;
    scatter3([dx/2 + stim.stimR(1)], [dx/2 + stim.stimR(2)], maxval, 400, 'Marker', 'x', ...
        'MarkerEdgeColor', 'white', 'LineWidth', 2);
    % if i == 1
    %     for k = 1:hetparam_uni.m
    %         quiver3(dx/2 + hetparam_uni.Ri(1, k), dx/2 + hetparam_uni.Ri(2, k), 0.01, ...
    %             hetparam_uni.Rf(1, k) - hetparam_uni.Ri(1, k), ...
    %             hetparam_uni.Rf(2, k) - hetparam_uni.Ri(2, k), 0, ...
    %             'Color', 'g', 'LineWidth', 0.5, ...
    %             'MaxHeadSize', 0.05 / norm(hetparam_uni.Ri(:, k) - hetparam_uni.Rf(:, k)), ...
    %             'Marker', '.', 'MarkerSize', 5, ...
    %             'AutoScale','off');
    %     end
    % end
    title(titles(i));
    view(0,90);
    colormap(CustomColormap);
    shading flat;
    clim([-maxval maxval]);
    xlim([0 topology.L]);
    ylim([0 topology.L]);
    xticks([]);
    yticks([]);
    xlabel('t = -2.0');
end

gif('50shortcutsBOLD.gif', 'DelayTime', 1/2)
for count = -3:1:size(boldfull,3)
    for i = 1:2
        set(h(i), 'Zdata', boldfull(:,:,max(count,1),i));
        xlabel(subplot(1, 2, i), append('t = ', sprintf('%.1f', count*0.5)));
    end
    pause(0.5);
    gif;
end

%% Animate different numbers of BIDIRECTIONAL shortcuts

clear; clc;

loadparam;

% Increase simulation time
dt = topology.T / topology.Nt;
topology.T = 0.1;
topology.Nt = ceil(topology.T / dt);

load 'CustomColormap.mat'

dx = topology.L / topology.Nx;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

region_radius = hetparam_uni.Ri(3, 1);

num_shortcuts_array = [10, 20, 50];

titles = ["10 shortcuts" "20 shortcuts" "50 shortcuts"];

ts = zeros(topology.Nx, topology.Nx, topology.Nt, 3);
for i = 1:3
    num_pairs = num_shortcuts_array(i);
    hetparam_bi.m = 2*num_pairs;
    hetparam_bi.c = ones(1, 2*num_pairs);
    hetparam_bi.tau = zeros(1, 2*num_pairs);
    hetparam_bi.Ri = zeros(3, 2*num_pairs);
    hetparam_bi.Rf = zeros(3, 2*num_pairs);
    hetparam_bi.Ri(3, :) = region_radius;
    hetparam_bi.Rf(3, :) = region_radius;
        
    load(append('r_centre-',num2str(num_pairs), 'pairs.mat'));
    
    r_centre = dx * ceil(r_centre * topology.L / dx);
    
    for k = 1:num_pairs
        hetparam_bi.Ri(1:2, 2*k-1) = r_centre(2*k-1, :);
        hetparam_bi.Rf(1:2, 2*k-1) = r_centre(2*k, :);
        hetparam_bi.Ri(1:2, 2*k) = r_centre(2*k, :);
        hetparam_bi.Rf(1:2, 2*k) = r_centre(2*k-1, :);
    end
    hetparam_array(i) = hetparam_bi;
    ts_values = run_periodic(topology,homparam,hetparam_bi,stim);
    ts(:, :, :, i) = ts_values;
end

%%
maxval = 5*max(ts(:, :, end, :), [], 'all');
maxvalall = 5*max(ts(:, :, :, :), [], 'all');

f = figure;
sgtitle('Response to Impulse Stimulation')
f.Position = [10 100 1400 400];
set(gcf, 'color', 'white');

formatSpec = '%.3f';

for i = 1:3
    subplot(1, 3, i);
    h(i) = surf(X, Y, 0*ts(:, :, 1, i));
    hold on;
    scatter3([stim.stimR(1)], [stim.stimR(2)], 0, 400, 'Marker', 'x', ...
        'MarkerEdgeColor', 'white', 'LineWidth', 2);
    scatter3(topology.L/2, topology.L/2, ...
        maxvalall, 600, 'Marker', '.', 'MarkerFaceColor', 'white');
    hetparam_bi = hetparam_array(i);
    for k = 1:hetparam_bi.m
        quiver3(dx/2 + hetparam_bi.Ri(1, k), dx/2 + hetparam_bi.Ri(2, k), 0, ...
            hetparam_bi.Rf(1, k) - hetparam_bi.Ri(1, k), ...
            hetparam_bi.Rf(2, k) - hetparam_bi.Ri(2, k), 0, ...
            'Color', 'g', 'LineWidth', 0.5, ...
            'MaxHeadSize', 0.05 / norm(hetparam_bi.Ri(:, k) - hetparam_bi.Rf(:, k)), ...
            'Marker', '.', 'MarkerSize', 5, ...
            'AutoScale','off');
    end
    title(titles(i));
    xlabel(append('t = ', num2str(-200*dt, formatSpec)));
    view(0,90);
    colormap(CustomColormap);
    shading flat;
    clim([-maxval maxval]);
    xlim([0 topology.L]);
    ylim([0 topology.L]);
    xticks([]);
    yticks([]);
end

gif('HeterogeneousBidirectional.gif','DelayTime',1/24)
for count = -200:5:ceil(topology.Nt/2)
    for i = 1:3
        if count > 0
            set(h(i), 'Zdata', ts(:,:,count,i))
        end
        xlabel(subplot(1, 3, i), append('t = ', num2str(count*dt, formatSpec)));
    end
    pause(0);
    gif
end

%% Animate BOLD response with different numbers of UNIDIRECTIONAL shortcuts
% AFTER A REAL BALLOON FILTER transform

clear; clc;

loadparam;

load 'CustomColormap.mat'

topology.T = 0.01;
topology.Nx = 100;
dx = topology.L / topology.Nx;
dt = dx / (homparam.r * homparam.gamma * sqrt(2) * 2);
topology.Nt = ceil(topology.T / dt);
clear dx dt;

stim.stimI = 0.01;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

region_radius = hetparam_uni.Ri(3, 1);

num_shortcuts_array = [10, 20, 50];

titles = ["10 shortcuts" "20 shortcuts" "50 shortcuts"];

boldfull = zeros(topology.Nx, topology.Nx, 14, 3);

for i = 1:3
    num_pairs = num_shortcuts_array(i);
    hetparam_uni.m = 2*num_pairs;
    hetparam_uni.c = ones(1, 2*num_pairs);
    hetparam_uni.tau = zeros(1, 2*num_pairs);
    hetparam_uni.Ri = zeros(3, 2*num_pairs);
    hetparam_uni.Rf = zeros(3, 2*num_pairs);
    hetparam_uni.Ri(3, :) = region_radius;
    hetparam_uni.Rf(3, :) = region_radius;

    load(append('r_centre-',num2str(num_pairs), 'pairs.mat'));
    r_centre = dx * ceil(r_centre * topology.L / dx);
    for k = 1:num_pairs
        hetparam_uni.Ri(1:2, k) = r_centre(2*k-1, :);
        hetparam_uni.Rf(1:2, k) = r_centre(2*k, :);
    end
    hetparam_array(i) = hetparam_uni;

    bold_values = [];
    init_het = 0;
    s_het_init = 0;
    f_het_init = 1;
    v_het_init = 1;
    q_het_init = 1;

    stim.stimnum = 1;
    
    for iter = 1:700
        ts_het = run_init_periodic(topology,homparam,hetparam_uni,stim,init_het);
        init_het = ts_het(:, :, end-1:end);
        [bold_het, s_het_init, f_het_init, v_het_init, q_het_init] = ...
            balloonmodel(ts_het, dt, s_het_init, f_het_init, v_het_init, q_het_init);
        if mod(iter, 50) == 0
            bold_values = cat(3, bold_values, ...
                cat(4, bold_het(:, :, end)));
        end
        disp(iter);
        stim.stimnum = 0;
    end

    boldfull(:, :, :, i) = bold_values;
end

%%
maxval = max(abs(boldfull(:, :, :, :)), [], "all");

f = figure;
sgtitle('BOLD Response to Impulse Stimulation')
f.Position = [10 100 1400 400];
set(gcf, 'color', 'white');

formatSpec = '%.3f';

for i = 1:3
    subplot(1, 3, i);
    h(i) = surf(X, Y, 0*boldfull(:, :, 1, i));
    hold on;
    scatter3([dx/2 + stim.stimR(1)], [dx/2 + stim.stimR(2)], 0, 400, 'Marker', 'x', ...
        'MarkerEdgeColor', 'white', 'LineWidth', 2);
    hetparam_uni = hetparam_array(i);
    for k = 1:hetparam_uni.m
        quiver3(dx/2 + hetparam_uni.Ri(1, k), dx/2 + hetparam_uni.Ri(2, k), 0, ...
            hetparam_uni.Rf(1, k) - hetparam_uni.Ri(1, k), ...
            hetparam_uni.Rf(2, k) - hetparam_uni.Ri(2, k), 0, ...
            'Color', 'g', 'LineWidth', 0.5, ...
            'MaxHeadSize', 0.05 / norm(hetparam_uni.Ri(:, k) - hetparam_uni.Rf(:, k)), ...
            'Marker', '.', 'MarkerSize', 5, ...
            'AutoScale','off');
    end
    title(titles(i));
    view(0,90);
    colormap(CustomColormap);
    shading flat;
    clim([-maxval maxval]);
    xlim([0 topology.L]);
    ylim([0 topology.L]);
    xticks([]);
    yticks([]);
    xlabel('t = -2.0');
end

gif('UnidirectionalBOLD.gif','DelayTime',1/2)
for count = -3:1:size(boldfull,3)
    for i = 1:3
        if count > 0
            set(h(i), 'Zdata', boldfull(:,:,count,i));
        end
        xlabel(subplot(1, 3, i), append('t = ', sprintf('%.1f', count*0.5)));
    end
    pause(0.5);
    gif;
end



%% Animate BOLD response with different numbers of BIDIRECTIONAL shortcuts
% AFTER A REAL BALLOON FILTER transform

clear; clc;

loadparam;

load 'CustomColormap.mat'

topology.T = 0.01;
topology.Nx = 100;
dx = topology.L / topology.Nx;
dt = dx / (homparam.r * homparam.gamma * sqrt(2) * 2);
topology.Nt = ceil(topology.T / dt);
clear dx dt;

stim.stimI = 0.01;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

region_radius = hetparam_uni.Ri(3, 1);

num_shortcuts_array = [0, 10, 50];

titles = ["0 shortcuts (homogeneous)" "10 shortcuts" "50 shortcuts"];

boldfull = zeros(topology.Nx, topology.Nx, 14, 4);

for i = 1:3

    if i == 1
        hetparam_bi2 = hetparam_hom;
    else
        num_pairs = num_shortcuts_array(i);
        hetparam_bi2 = hetparam_bi;
        hetparam_bi2.m = 2*num_pairs;
        hetparam_bi2.c = ones(1, 2*num_pairs);
        hetparam_bi2.tau = zeros(1, 2*num_pairs);
        hetparam_bi2.Ri = zeros(3, 2*num_pairs);
        hetparam_bi2.Rf = zeros(3, 2*num_pairs);
        hetparam_bi2.Ri(3, :) = region_radius;
        hetparam_bi2.Rf(3, :) = region_radius;
    
        load(append('r_centre-',num2str(num_pairs), 'pairs.mat'));
        r_centre = dx * ceil(r_centre * topology.L / dx);
        for k = 1:num_pairs
            hetparam_bi2.Ri(1:2, 2*k-1) = r_centre(2*k-1, :);
            hetparam_bi2.Rf(1:2, 2*k-1) = r_centre(2*k, :);
            hetparam_bi2.Ri(1:2, 2*k) = r_centre(2*k, :);
            hetparam_bi2.Rf(1:2, 2*k) = r_centre(2*k-1, :);
        end
    end

    hetparam_array(i) = hetparam_bi2;

    bold_values = [];
    init_het = 0;
    s_het_init = 0;
    f_het_init = 1;
    v_het_init = 1;
    q_het_init = 1;

    stim.stimnum = 1;
    
    for iter = 1:700
        ts_het = run_init_periodic(topology,homparam,hetparam_bi2,stim,init_het);
        init_het = ts_het(:, :, end-1:end);
        [bold_het, s_het_init, f_het_init, v_het_init, q_het_init] = ...
            balloonmodel(ts_het, dt, s_het_init, f_het_init, v_het_init, q_het_init);
        if mod(iter, 50) == 0
            bold_values = cat(3, bold_values, ...
                cat(4, bold_het(:, :, end)));
        end
        disp(iter);
        stim.stimnum = 0;
    end

    boldfull(:, :, :, i) = bold_values;
end

%%
maxvalall = max(boldfull(:, :, :, :), [], 'all');

f = figure;
sgtitle('Observed fMRI Response to Impulse Stimulation')
f.Position = [10 100 1400 400];
set(gcf, 'color', 'white');

formatSpec = '%.0f';

for i = 1:3
    subplot(1, 3, i);
    h(i) = surf(X, Y, 0*boldfull(:, :, 1, 1));
    hold on;
    scatter3([stim.stimR(1)], [stim.stimR(2)], 0, 400, 'Marker', 'x', ...
        'MarkerEdgeColor', 'white', 'LineWidth', 2);
    if i > 1
        hetparam_bi = hetparam_array(i);
        for k = 1:num_shortcuts_array(i)
            plot3(dx/2 + [hetparam_bi.Ri(1, 2*k-1),hetparam_bi.Rf(1, 2*k-1)], ...
                dx/2 + [hetparam_bi.Ri(2, 2*k-1),hetparam_bi.Rf(2, 2*k-1)], ...
                [0, 0], ...
                'Color', 'g', 'LineWidth', 1);
        end
    end
    title(titles(i));
    xlabel('t = -2.0 s');
    view(0,90);
    colormap(CustomColormap);
    shading flat;
    clim([-maxvalall maxvalall]);
    xlim([dx topology.L]);
    ylim([dx topology.L]);
    xticks([]);
    yticks([]);
end

gif('BidirectionalBOLD.gif','DelayTime',1/2)
for count = -3:1:size(boldfull,3)
    for i = 1:3
        if count > 0
            set(h(i), 'Zdata', boldfull(:,:,count,i))
        end
        xlabel(subplot(1, 3, i), append('t = ', sprintf('%.1f', count*0.5), ' s'));
    end
    pause(0.5);
    gif
end


%% Animation of RANDOM 50 shortcut model under 
% AFTER A REAL BALLOON FILTER transform

clear; clc;

loadparam;

load 'CustomColormap.mat'

topology.T = 0.01;
topology.Nx = 50;
dx = topology.L / topology.Nx;
dt = dx / (homparam.r * homparam.gamma * sqrt(2) * 2);
topology.Nt = ceil(topology.T / dt);
clear dx dt;

stim.stimI = 0.01;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

region_radius = hetparam_uni.Ri(3, 1);
min_distance = 0.09;
max_distance = topology.L / 2;

lambda = 188;

num_shortcuts = 50;
hetparam_uni.m = num_shortcuts;
hetparam_uni.c = ones(1, num_shortcuts);
hetparam_uni.tau = zeros(1, num_shortcuts);
hetparam_uni.Ri = zeros(3, num_shortcuts);
hetparam_uni.Rf = zeros(3, num_shortcuts);
hetparam_uni.Ri(3, :) = region_radius;
hetparam_uni.Rf(3, :) = region_radius;

rmid_centre = dx*randi(topology.Nx, 1, 2);
    
r_centre = zeros(2*num_shortcuts, 2);

for k = 1:num_shortcuts
    
    % Add ith shortcut

    dist = Inf;
    while (dist < min_distance || dist > max_distance)
        % ri_centre = dx*randi(topology.Nx, 1, 2);
        % rf_centre = dx*randi(topology.Nx, 1, 2);
        if k <= num_shortcuts / 3
            ri_centre = rmid_centre;
            rf_centre = dx*randi(topology.Nx, 1, 2);
        else
            ri_centre = dx*randi(topology.Nx, 1, 2);
            rf_centre = rmid_centre;
        end
        distx = abs(ri_centre(1) - rf_centre(1));
        if distx > topology.L / 2
            distx = topology.L - distx;
        end
        disty = abs(ri_centre(2) - rf_centre(2));
        if disty > topology.L / 2
            disty = topology.L - disty;
        end
        dist = norm([distx; disty], 2);
        unif_rv = rand();
        if unif_rv > (exp(-dist/lambda)/dist) / (exp(-min_distance/lambda)/min_distance)
            dist = Inf;
        end
    end
    hetparam_uni.Ri(1:2, k) = ri_centre;
    hetparam_uni.Rf(1:2, k) = rf_centre;
end

boldfull = [];
init_hom = 0;
init_het = 0;
s_hom_init = 0;
f_hom_init = 1;
v_hom_init = 1;
q_hom_init = 1;
s_het_init = 0;
f_het_init = 1;
v_het_init = 1;
q_het_init = 1;

for iter = 1:150
    ts_hom = run_init_periodic(topology,homparam,hetparam_hom,stim,init_hom);
    init_hom = ts_hom(:, :, end-1:end);
    ts_het = run_init_periodic(topology,homparam,hetparam_uni,stim,init_het);
    init_het = ts_het(:, :, end-1:end);
    [bold_hom, s_hom_init, f_hom_init, v_hom_init, q_hom_init] = ...
        balloonmodel(ts_hom, dt, s_hom_init, f_hom_init, v_hom_init, q_hom_init);
    [bold_het, s_het_init, f_het_init, v_het_init, q_het_init] = ...
        balloonmodel(ts_het, dt, s_het_init, f_het_init, v_het_init, q_het_init);
    if mod(iter, 50) == 0
        boldfull = cat(3, boldfull, ...
            cat(4, bold_het(:, :, end), bold_hom(:, :, end)));
        disp(iter);
    end
    stim.stimnum = 0;
end

S_u = squeeze(sum(boldfull(:, :, end, 1), [1 2]));
S_v = squeeze(sum(boldfull(:, :, end, 2), [1 2]));
S_u2 = squeeze(sum(boldfull(:, :, end, 1).^2, [1 2]));
S_v2 = squeeze(sum(boldfull(:, :, end, 2).^2, [1 2]));
S_uv = squeeze(sum(boldfull(:, :, end, 1).*boldfull(:, :, end, 2), [1 2]));

num = (topology.Nx)^4 * (S_uv).^2 - (S_u).^2 .* (S_v).^2;
denom = sqrt((topology.Nx)^4 * (S_u2).^2 - (S_u).^4) .* sqrt((topology.Nx)^4 * (S_v2).^2 - (S_v).^4);

dissim = 1 - num ./ denom;
disp(dissim);

disp(pdist2(reshape(boldfull(:, :, end, 1), [], 1)', reshape(boldfull(:, :, end, 2), [], 1)', 'correlation'))

%%

maxval = max(abs(boldfull(:, :, :, :)), [], "all");

titles = [append("Heterogeneous - 50 shortcuts"), "Homogeneous Approximation"];

f = figure;
sgtitle('BOLD Response to Impulse Stimulation')
f.Position = [10 100 1400 600];
set(gcf, 'color', 'white');

for i = 1:2
    subplot(1, 2, i);
    h(i) = surf(X, Y, 0*boldfull(:, :, 1, i));
    hold on;
    scatter3([dx/2 + stim.stimR(1)], [dx/2 + stim.stimR(2)], maxval, 400, 'Marker', 'x', ...
        'MarkerEdgeColor', 'white', 'LineWidth', 2);
    if i == 1
        for k = 1:hetparam_uni.m
            quiver3(dx/2 + hetparam_uni.Ri(1, k), dx/2 + hetparam_uni.Ri(2, k), 0, ...
                hetparam_uni.Rf(1, k) - hetparam_uni.Ri(1, k), ...
                hetparam_uni.Rf(2, k) - hetparam_uni.Ri(2, k), 0, ...
                'Color', 'g', 'LineWidth', 0.5, ...
                'MaxHeadSize', 0.05 / norm(hetparam_uni.Ri(:, k) - hetparam_uni.Rf(:, k)), ...
                'Marker', '.', 'MarkerSize', 5, ...
                'AutoScale','off');
        end
    end
    title(titles(i));
    view(0,90);
    colormap(CustomColormap);
    shading flat;
    clim([-maxval maxval]);
    xlim([0 topology.L]);
    ylim([0 topology.L]);
    xticks([]);
    yticks([]);
    xlabel('t = -2.0');
end

%gif('50shortcutsBOLD.gif', 'DelayTime', 1/2)
for count = -3:1:size(boldfull,3)
    for i = 1:2
        if count > 0
            set(h(i), 'Zdata', boldfull(:,:,max(count,1),i));
        end
        xlabel(subplot(1, 2, i), append('t = ', sprintf('%.1f', count*0.5)));
    end
    pause(0.5);
    %gif;
end

%% Animate heterogeneous model with central hub

clear; clc;

loadparam;

load 'CustomColormap.mat'

% Increase simulation time
topology.Nx = 50;
topology.T = 0.1;
dx = topology.L / topology.Nx;
dt = dx / (homparam.r * homparam.gamma * sqrt(2) * 10);
topology.Nt = ceil(topology.T / dt);
dt = topology.T / topology.Nt;

dx = topology.L / topology.Nx;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

hub_centre = dx*randi(topology.Nx, 1, 2);

region_radius = hetparam_uni.Ri(3, 1);
min_distance = 0.09;
max_distance = topology.L / 2;
lambda = 188;
num_shortcuts = 50;
in_degree = 10;

r_centre = zeros(num_shortcuts, 2);

for j = 1:num_shortcuts
    
    % Add ith shortcut

    dist = Inf;
    while (dist < min_distance || dist > max_distance)
        local_centre = dx*randi(topology.Nx, 1, 2);
        distx = abs(local_centre(1) - hub_centre(1));
        if distx > topology.L / 2
            distx = topology.L - distx;
        end
        disty = abs(local_centre(2) - hub_centre(2));
        if disty > topology.L / 2
            disty = topology.L - disty;
        end
        dist = norm([distx; disty], 2);
        unif_rv = rand();
        if unif_rv > (exp(-dist/lambda)/dist) / (exp(-min_distance/lambda)/min_distance)
            dist = Inf;
        end
    end
    r_centre(j, :) = local_centre;
end

hetparam_uni2 = hetparam_uni;
hetparam_uni2.m = num_shortcuts;
hetparam_uni2.c = ones(1, num_shortcuts);
hetparam_uni2.tau = zeros(1, num_shortcuts);
hetparam_uni2.Ri = zeros(3, num_shortcuts);
hetparam_uni2.Rf = zeros(3, num_shortcuts);
hetparam_uni2.Ri(3, :) = region_radius;
hetparam_uni2.Rf(3, :) = region_radius;

for j2 = 1:num_shortcuts
    if j2 <= in_degree
        hetparam_uni2.Ri(1:2, j2) = r_centre(j2, :);
        hetparam_uni2.Rf(1:2, j2) = hub_centre;
    else
        hetparam_uni2.Ri(1:2, j2) = hub_centre;
        hetparam_uni2.Rf(1:2, j2) = r_centre(j2, :);
    end
end

ts = run_periodic(topology,homparam,hetparam_uni2,stim);

maxval = max(ts(:, :, end), [], 'all');
maxvalall = max(ts(:, :, :), [], 'all');

f = figure;
sgtitle('Response to Impulse Stimulation')
f.Position = [10 100 425 400];
set(gcf, 'color', 'white');

formatSpec = '%.3f';

subplot(1, 1, 1);
h(1) = surf(X, Y, 0*ts(:, :, 1));
hold;
scatter3([stim.stimR(1)], [stim.stimR(2)], 0, 400, 'Marker', 'x', ...
        'MarkerEdgeColor', 'white', 'LineWidth', 2);
xlabel(append('t = ', num2str(-200*dt, formatSpec)));
view(0,90);
colormap(CustomColormap);
shading flat;
clim([-maxval maxval]);
xlim([0 topology.L]);
ylim([0 topology.L]);
xticks([]);
yticks([]);

%gif('Homogeneous.gif','DelayTime',1/24)

for count = -200:5:topology.Nt
    if count > 0
        set(h(1), 'Zdata', ts(:,:,count));
    end
    xlabel(subplot(1, 1, 1), append('t = ', num2str(count*dt, formatSpec)));
    pause(0.01);
    %gif;
end

%% Stimulate heterogeneous model with 50 shortcuts with custom input

clear; clc;

loadparam;

load 'CustomColormap.mat'

topology.Nx = 100;
topology.T = 1;

dx = topology.L / topology.Nx;
dt = dx / (homparam.r * homparam.gamma * sqrt(2) * 2);
topology.Nt = ceil(topology.T / dt);
dt = topology.T / topology.Nt;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

% Create 50 shortcut struct
num_shortcuts = 50;
hetparam_bi.m = 2*num_shortcuts;
hetparam_bi.c = ones(1, 2*num_shortcuts);
hetparam_bi.tau = zeros(1, 2*num_shortcuts);
region_radius = hetparam_bi.Ri(3, 1);
hetparam_bi.Ri = zeros(3, 2*num_shortcuts);
hetparam_bi.Rf = zeros(3, 2*num_shortcuts);
hetparam_bi.Ri(3, :) = 2*region_radius;
hetparam_bi.Rf(3, :) = 2*region_radius;
    
load(append('r_centre-',num2str(num_shortcuts), 'pairs.mat'));

r_centre = dx * ceil(r_centre * topology.L / dx);

for k = 1:num_shortcuts
    hetparam_bi.Ri(1:2, 2*k-1) = r_centre(2*k-1, :);
    hetparam_bi.Rf(1:2, 2*k-1) = r_centre(2*k, :);
    hetparam_bi.Ri(1:2, 2*k) = r_centre(2*k, :);
    hetparam_bi.Rf(1:2, 2*k) = r_centre(2*k-1, :);
end

position_centre = [round(topology.Nx / 2), round(topology.Nx / 2)];

region_radius = 0.005;

stiminput = zeros(topology.Nx, topology.Nx, topology.Nt);
for i = 1:topology.Nx
    for j = 1:topology.Nx
        distx = abs(i - position_centre(1));
        if distx > topology.Nx / 2
            distx = topology.Nx - distx;
        end
        disty = abs(j - position_centre(2));
        if disty > topology.Nx / 2
            disty = topology.Nx - disty;
        end
        dist = sqrt(distx^2 + disty^2);
        if dx * dist < region_radius
            stiminput(i, j, :) = randn(1, topology.Nt);
        end
    end
end
ts = run_init_cstminput_periodic(topology,homparam,hetparam_bi,stiminput,0);

maxval = max(ts(:, :, end), [], 'all');
maxvalall = max(ts(:, :, :), [], 'all');

f = figure;
sgtitle('Response to Impulse Stimulation')
f.Position = [10 100 425 400];
set(gcf, 'color', 'white');

formatSpec = '%.3f';

h(1) = surf(X, Y, ts(:, :, 1));
hold;
xlabel(append('t = ', num2str(-200*dt, formatSpec)));
view(0,90);
colormap(CustomColormap);
shading flat;
clim([-maxval maxval]);
xlim([0 topology.L]);
ylim([0 topology.L]);
xticks([]);
yticks([]);

%gif('Homogeneous.gif','DelayTime',1/24)

for count = -0:1:topology.Nt
    if count > 0
        set(h(1), 'Zdata', ts(:,:,count));
    end
    xlabel(append('t = ', num2str(count*dt, formatSpec)));
    pause(0.01);
    %gif;
end

%% Calculate SCF of heterogeneous model vs input radius
% All points in stimulated region are stimulated with white noise

% Preparation

clear; clc;

loadparam;

topology.T = topology.T * 2;
topology.Nt = topology.Nt * 2;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

% Set max_distance of stimulation point from Ri

position_centre = hetparam_het.a;
position_centre_grid = position_centre / dx;

hetparam_het.c = hetparam_het.c * 3;
hetparam_het.sigmaeps = 0.002;

region_radius_array = [0.005, 0.05, Inf];
num_experiments = length(region_radius_array);

ts_full = zeros(topology.Nx, topology.Nx, topology.Nt, num_experiments);

% num_lrcs = 100;
% 
% % Create LRC topology
% rng(1, "twister");
% a_array = zeros(2, num_lrcs);
% b_array = zeros(2, num_lrcs);
% 
% % Sample connectome
% 
% for j = 1:num_lrcs
% 
%     a = topology.L*rand(2, 1);
%     b = topology.L*rand(2, 1);
%     a_array(:, j) = a;
%     b_array(:, j) = b;
% 
% end
% 
% hetparam_het.m = num_lrcs;
% hetparam_het.c = (homparam.r)^2 * ones(1, num_lrcs);
% hetparam_het.tau = zeros(1, num_lrcs);
% hetparam_het.a = a_array;
% hetparam_het.b = b_array;

sigma = 1;
rng(0, "twister");
stiminput0 = sigma * sqrt(dt) * randn(topology.Nx, topology.Nx, topology.Nt);


for count = 1:length(region_radius_array)

    % Set region radius
    
    region_radius = region_radius_array(count);
    
    init = 0;

    if region_radius < Inf
    
        stimulationinputspace = zeros(topology.Nx, topology.Nx);
        i0 = position_centre(1) / dx; j0 = position_centre(2) / dx;
        for i = 1:topology.Nx
            for j = 1:topology.Nx
                distx = dx*abs(i - i0);
                distx = min(distx, topology.L - distx);
                disty = dx*abs(j - j0);
                disty = min(disty, topology.L - disty);
                stimulationinputspace(i, j) = exp(-0.5*(distx^2 + disty^2) / (region_radius^2));
            end
        end
        stiminput = stiminput0 .* reshape(stimulationinputspace, [topology.Nx, topology.Nx, 1]);
    else
        stiminput = stiminput0;
    end
    ts = run_init_cstminput_periodic(topology,homparam,hetparam_het,stiminput,init);
    ts_full(:, :, :, count) = ts;
    disp(num2str([region_radius]));

end

%%

maxval = max(ts_full, [], [1 2 3]) / 8;
Colormap = [linspace(0, 1, 256)' linspace(0, 1, 256)', linspace(1, 1, 256)';
    linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];

f = figure;
sgtitle('Noise-Input Driven Activity (Circle = location of input)', 'Interpreter', 'Latex', 'FontSize', 20)
f.Position = [10 100 1400 410];
set(gcf, 'color', 'white');

min_n = round(topology.Nt / 2);

for iter = 1:num_experiments
    region_radius = region_radius_array(iter);
    ax = subplot(1, num_experiments, iter);
    hold on;
    ax.Box = "on";
    ax.LineWidth = 2;
    h(iter) = imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), ts_full(:, :, min_n, iter));
    clim([-maxval(iter), maxval(iter)] * (exp( 2 * iter / num_experiments)));
    %draw circle
    if region_radius < Inf
        theta = linspace(0,2*pi,1000);
        x = region_radius * cos(theta) + position_centre(1) + dx/2;
        y = region_radius * sin(theta) + position_centre(2) + dx/2;
        fill(x, y, 'k', 'FaceAlpha', 0., 'EdgeColor', 'k', 'LineWidth', 1);
    else
        fill([dx dx topology.L topology.L], [dx topology.L topology.L dx], 'k', 'FaceAlpha', 0., 'EdgeColor', 'k', 'LineWidth', 1)
    end
    quiver(hetparam_het.a(1), hetparam_het.a(2), ...
        hetparam_het.b(1) - hetparam_het.a(1), hetparam_het.b(2) - hetparam_het.a(2),...
        'Color', 'k', 'LineWidth', 1, ...
        'MaxHeadSize', 0.05 / norm(hetparam_het.b - hetparam_het.a), ...
        'Marker', '.', 'MarkerSize', 0.0001, ...
        'AutoScale','off');
    colormap(Colormap);
    xlim([0, topology.L + dx]);
    ylim([0, topology.L + dx]);
    xticks([]);
    yticks([]);
    if iter == 1
        title('Controlled Input', 'Interpreter', 'latex', 'FontSize', 15);
    end
    if iter == 3
        title('Spontaneous (Uncontrolled)', 'Interpreter', 'latex', 'FontSize', 15);
    end
end

gif('spontaneousposition.gif', 'DelayTime', 1/12)
for count = min_n:2:topology.Nt
    for iter = 1:num_experiments
        set(h(iter), 'Cdata', ts_full(:,:,count,iter));
    end
    % xlabel(append('t = ', num2str(count*dt, formatSpec)));
    pause(0.01);
    gif;
end


%% Plot dissim vs stim position
clear; clc; 
loadparam; 

load('divergencevsposition_default.mat');
dissim_array = reshape(dissim_array, [topology.Nx, topology.Nx]);
maxval = max(dissim_array, [], 'all');

dx = topology.L / topology.Nx;
[Y, X] = meshgrid(dx * (1:topology.Nx), dx * (1:topology.Nx));

f = figure;
f.Position = [100, 100, 360, 315];
ax = gca;
ax.Box = "on";
ax.LineWidth = 1;
hold;

imagesc(dx*(1:topology.Nx), dx*(1:topology.Nx), dissim_array');
set(ax, 'YDir', 'normal');
view(0, 90);
colorbar;
colormap([linspace(1, 0, 256); linspace(1, 0, 256); linspace(1, 1, 256)]')
clim([0 maxval]);
shading flat;
xlim([0, topology.L + dx]);
ylim([0, topology.L + dx]);
xticks([]);
yticks([]);
quiver(hetparam_het.a(1), hetparam_het.a(2), ...
    hetparam_het.b(1) - hetparam_het.a(1), ...
    hetparam_het.b(2) - hetparam_het.a(2),...
    'Color', 'k', 'LineWidth', 1, ...
    'MaxHeadSize', 0.05 / norm(hetparam_het.a - hetparam_het.b), ...
    'Marker', '.', 'MarkerSize', 0.0001, ...
    'AutoScale','off');
title("$\max C_\phi$", 'Interpreter', 'latex', 'FontSize', 20);          




