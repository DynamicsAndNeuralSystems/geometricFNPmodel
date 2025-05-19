%% Draw neural field model

clear; clc;
loadparam;

Nx = 41;

f = @(x, y) 1-2*(x.^2+y.^2);

% Step 1: Create the surface mesh
x_coord = linspace(-0.5,0.5,Nx);
y_coord = linspace(-0.5,0.5,Nx);
[X, Y] = meshgrid(x_coord, y_coord); % Generate grid data
Z = f(X,Y); % Calculate Z values

% Create a new figure
fig = figure;
hold;
fig.Position = [10 10 450 300];
view(35, 40);
axis off;

% Plot the surface mesh
surf(X, Y, Z, 'FaceColor', 'none', 'EdgeColor', [0.5 0.5 1], 'FaceAlpha', 0.1);
plot3(x_coord, -0.5*ones(1, Nx), f(x_coord, -0.5), 'k', 'LineWidth', 2);
plot3(x_coord, 0.5*ones(1, Nx), f(x_coord, 0.5), 'k', 'LineWidth', 2);
plot3(-0.5*ones(1, Nx), y_coord, f(-0.5, y_coord), 'k', 'LineWidth', 2);
plot3(0.5*ones(1, Nx), y_coord, f(0.5, y_coord), 'k', 'LineWidth', 2);

Ri = [-0.25, -0.25];
Rf = [0.25, 0.25];
LineZ = -1.2;
radius = 0.08;

% Draw Laplacian
%lap_x = 0.2; lap_y = -0.2; lap_z = f(lap_x, lap_y);
%dzx = -4*lap_x; dzy = -4*lap_y;
%quiver3(lap_x, lap_y, lap_z, 0.15, 0, 0.15*dzx, 'Color', [0 0 1], 'LineWidth', 2, 'MaxHeadSize', 5)
%quiver3(lap_x, lap_y, lap_z, -0.15, 0, -0.15*dzx, 'Color', [0 0 1], 'LineWidth', 2, 'MaxHeadSize', 5)
%quiver3(lap_x, lap_y, lap_z, 0, 0.15, 0.15*dzy, 'Color', [0 0 1], 'LineWidth', 2, 'MaxHeadSize', 5)
%quiver3(lap_x, lap_y, lap_z, 0, -0.15, -0.15*dzy, 'Color', [0 0 1], 'LineWidth', 2, 'MaxHeadSize', 5)

zlim([-2.7, Inf])

% Step 2: Set the renderer to 'painters'
set(gcf, 'Renderer', 'painters');

% Step 3: Export the plot as SVG
print(gcf, 'schematicfield.svg', '-dsvg');

% save as surface_mesh_plot.svg

%% Draw neural mass model

clear; clc;
loadparam;

Nx = 41;

f = @(x, y) 1-2*(x.^2+y.^2);

% Step 1: Create the surface mesh
x_coord = linspace(-0.5,0.5,Nx);
y_coord = linspace(-0.5,0.5,Nx);
[X, Y] = meshgrid(x_coord, y_coord); % Generate grid data
Z = f(X,Y); % Calculate Z values

% Create a new figure
fig = figure;
hold;
fig.Position = [10 10 450 300];
view(35, 40);
axis off;

% Plot the surface mesh
surf(X, Y, Z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1], 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1);
plot3(x_coord, -0.5*ones(1, Nx), f(x_coord, -0.5), 'k', 'LineWidth', 2);
plot3(x_coord, 0.5*ones(1, Nx), f(x_coord, 0.5), 'k', 'LineWidth', 2);
plot3(-0.5*ones(1, Nx), y_coord, f(-0.5, y_coord), 'k', 'LineWidth', 2);
plot3(0.5*ones(1, Nx), y_coord, f(0.5, y_coord), 'k', 'LineWidth', 2);

Ri_array = [-0.25, -0.25; -0.3, 0.2];
Rf_array = [0.25, 0.25; 0.3, -0.2];

depth_array = [0.7, 1];

for iter = 1:2
    Ri = Ri_array(iter, :);
    Rf = Rf_array(iter, :);
    LineZ = 2;
    radius = 0.08;
    
    % Draw shortcuts
    % plot3([x_coord(i), Ri(1)], [y_coord(j), Ri(2)], ...
    %     [f(x_coord(i), y_coord(j)), -LineZ], ...
    %     'Color', [0 0.5 0 exp(-4000*(norm([x_coord(i), y_coord(j)] - Ri))^2)], 'LineWidth', 2);
    % plot3([x_coord(i), Rf(1)], [y_coord(j), Rf(2)], ...
    %     [f(x_coord(i), y_coord(j)), -LineZ], ...
    %     'Color', [0 0.5 0 exp(-4000*(norm([x_coord(i), y_coord(j)] - Rf))^2)], 'LineWidth', 2);
    plot3(linspace(Ri(1), Rf(1), 1000), linspace(Ri(2), Rf(2), 1000), f(Ri(1),Ri(2)) - depth_array(iter)*(1-linspace(-1, 1, 1000).^2).^(1/8), ...
        'Color', [0 0.5 0], 'LineWidth', 2)
    % x_circle = Ri(1) + radius*cos(linspace(0, 2*pi, 40));
    % y_circle = Ri(2) + radius*sin(linspace(0, 2*pi, 40));
    % plot3(x_circle, y_circle, f(x_circle, y_circle), 'k', 'LineWidth', 2);
    % x_circle = Rf(1) + radius*cos(linspace(0, 2*pi, 40));
    % y_circle = Rf(2) + radius*sin(linspace(0, 2*pi, 40));
    % plot3(x_circle, y_circle, f(x_circle, y_circle), 'k', 'LineWidth', 2);

end

zlim([-2.7, Inf])

% Step 2: Set the renderer to 'painters'
set(gcf, 'Renderer', 'painters');

% Step 3: Export the plot as SVG
print(gcf, 'schematicmass.svg', '-dsvg');

% save as schematicmass.svg

%% Draw mesh of square with a long range connection

clear; clc;
loadparam;

Nx = 41;

f = @(x, y) 1-2*(x.^2+y.^2);

% Step 1: Create the surface mesh
x_coord = linspace(-0.5,0.5,Nx);
y_coord = linspace(-0.5,0.5,Nx);
[X, Y] = meshgrid(x_coord, y_coord); % Generate grid data
Z = f(X,Y); % Calculate Z values

% Create a new figure
fig = figure;
hold;
fig.Position = [10 10 900 600];
view(35, 40);
axis off;

% Plot the surface mesh
surf(X, Y, Z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

Ri = [-0.25, -0.25];
Rf = [0.25, 0.25];
depth = 2;
C_arrow_pos = 0.08;

% Draw shortcuts
plot3(linspace(Ri(1), Rf(1), 1000), linspace(Ri(2), Rf(2), 1000), f(Ri(1),Ri(2)) - depth*(1-linspace(-1, 1, 1000).^2).^(1/8), ...
    'Color', 'k', 'LineWidth', 1)
plot3(C_arrow_pos*[-1 1], C_arrow_pos*[-1 1], (f(Ri(1),Ri(2))-1.1*depth) + [0 0], 'Color', [0 0.5 0], 'LineWidth', 2)
plot3(C_arrow_pos*[1 2/3], C_arrow_pos*[1 2/3], (f(Ri(1),Ri(2))-1.1*depth) + [0 0.075], 'Color', [0 0.5 0], 'LineWidth', 2)
plot3(C_arrow_pos*[1 2/3], C_arrow_pos*[1 2/3], (f(Ri(1),Ri(2))-1.1*depth) + [0 -0.075], 'Color', [0 0.5 0], 'LineWidth', 2)

% Draw circles around endpoints
theta = linspace(-pi,pi,1000);
x = 0.005 * cos(theta) + Ri(1);
y = 0.005 * sin(theta) + Ri(2);
plot3(x, y, f(x, y), 'Color', 'k', 'LineWidth', 2);
x = 0.005 * cos(theta) + Rf(1);
y = 0.005 * sin(theta) + Rf(2);
plot3(x, y, f(x, y), 'Color', 'k', 'LineWidth', 2);


% Draw Laplacian
lap_x = 0.; lap_y = -0.; lap_z = f(lap_x, lap_y);
dzx = -4*lap_x; dzy = -4*lap_y;
quiver3(lap_x, lap_y, lap_z, 0.15, 0, 0.15*dzx, 'Color', [0 0 1], 'LineWidth', 2, 'MaxHeadSize', 5)
quiver3(lap_x, lap_y, lap_z, -0.15, 0, -0.15*dzx, 'Color', [0 0 1], 'LineWidth', 2, 'MaxHeadSize', 5)
quiver3(lap_x, lap_y, lap_z, 0, 0.15, 0.15*dzy, 'Color', [0 0 1], 'LineWidth', 2, 'MaxHeadSize', 5)
quiver3(lap_x, lap_y, lap_z, 0, -0.15, -0.15*dzy, 'Color', [0 0 1], 'LineWidth', 2, 'MaxHeadSize', 5)

zlim([-2.7, Inf])

% Step 2: Set the renderer to 'painters'
set(gcf, 'Renderer', 'painters');

% Step 3: Export the plot as SVG
print(gcf, 'schematicfieldmass.svg', '-dsvg');

% save as surface_mesh_plot.svg