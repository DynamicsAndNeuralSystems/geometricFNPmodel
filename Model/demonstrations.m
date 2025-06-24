%% Animate geometric model

clear; clc;
% Load model parameters
loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

ts = run_periodic(topology,homparam,hetparam_hom,stim);

% Create figure
f = figure;
sgtitle('Evoked Response - Geometric Connectivity')
f.Position = [300 100 650 600];
set(gcf, 'color', 'white');
maxval = 650;
formatSpec = '%.0f';

Colormap = [linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];

% Set heatmap
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

% Animate
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
% Load model parameters
loadparam;

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

% Create discrete connectome of 50 FNPs
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

% Create FNP struct
hetparam_het.m = num_fnps;
hetparam_het.c = (homparam.r)^2 * ones(1, num_fnps);
hetparam_het.tau = zeros(1, num_fnps);
hetparam_het.a = a_array;
hetparam_het.b = b_array;

ts = run_periodic(topology,homparam,hetparam_het,stim);

%
f = figure;
sgtitle('Evoked Response - Hybrid Connectivity')
f.Position = [300 100 650 600];
set(gcf, 'color', 'white');
maxval = 400;
formatSpec = '%.0f';

Colormap = [linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];

% Set heatmap
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
    'MaxHeadSize', 0.05 / norm(hetparam_het.a(:, k) - hetparam_het.b(:, k)), ...
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

% Animate
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

