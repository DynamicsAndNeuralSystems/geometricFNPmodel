%% Draw neural field model

clear; clc;
loadparam;


% Create function for surface
f = @(x, y) 1-2*(x.^2+y.^2);

% Create surface mesh
Nx = 41; % 41 x 41 grid
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

% Set limits
zlim([-2.7, Inf])

% Export the plot as schematicfield.svg
set(gcf, 'Renderer', 'painters');
print(gcf, '1_schematicfield.svg', '-dsvg');
close(fig);

%% Draw neural mass model with FNPs (fast-conducting non-local projections)

clear; clc;
loadparam;

% Create function for surface
f = @(x, y) 1-2*(x.^2+y.^2);

% Create the surface mesh
Nx = 41; % 41 x 41 grid
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

% Specify the endpoints (a,b) of two FNPs
a_array = [-0.25, -0.25; -0.3, 0.2];
b_array = [0.25, 0.25; 0.3, -0.2];

% Specify the depth of each of the FNPs
depth_array = [0.7, 1];

% Draw FNPs
for iter = 1:size(a_array, 1)
    a = a_array(iter, :);
    b = b_array(iter, :);
    plot3(linspace(a(1), b(1), 1000), linspace(a(2), b(2), 1000), f(a(1),a(2)) - depth_array(iter)*(1-linspace(-1, 1, 1000).^2).^(1/8), ...
        'Color', [0 0.5 0], 'LineWidth', 2);
end

% Set limits
zlim([-2.7, Inf])

% Export the plot as schematicmass.svg
set(gcf, 'Renderer', 'painters');
print(gcf, '1_schematicmass.svg', '-dsvg');
close(fig);

%% Draw mesh of square with a long range connection

clear; clc;
loadparam;

% Define function for surface
f = @(x, y) 1-2*(x.^2+y.^2);

% Create the surface mesh
Nx = 41; % 41 x 41 grid
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

% Specify endpoints and depth of a single FNP
a = [-0.25, -0.25];
b = [0.25, 0.25];
depth = 2;

% Specify location of arrow
C_arrow_pos = 0.08;

% Draw FNP
plot3(linspace(a(1), b(1), 1000), linspace(a(2), b(2), 1000), f(a(1),a(2)) - depth*(1-linspace(-1, 1, 1000).^2).^(1/8), ...
    'Color', 'k', 'LineWidth', 1)
plot3(C_arrow_pos*[-1 1], C_arrow_pos*[-1 1], (f(a(1),a(2))-1.1*depth) + [0 0], 'Color', [0 0.5 0], 'LineWidth', 2)
plot3(C_arrow_pos*[1 2/3], C_arrow_pos*[1 2/3], (f(a(1),a(2))-1.1*depth) + [0 0.075], 'Color', [0 0.5 0], 'LineWidth', 2)
plot3(C_arrow_pos*[1 2/3], C_arrow_pos*[1 2/3], (f(a(1),a(2))-1.1*depth) + [0 -0.075], 'Color', [0 0.5 0], 'LineWidth', 2)

% Draw circles on endpoints
theta = linspace(-pi,pi,1000);
x = 0.005 * cos(theta) + a(1);
y = 0.005 * sin(theta) + a(2);
plot3(x, y, f(x, y), 'Color', 'k', 'LineWidth', 2);
x = 0.005 * cos(theta) + b(1);
y = 0.005 * sin(theta) + b(2);
plot3(x, y, f(x, y), 'Color', 'k', 'LineWidth', 2);

% Draw Laplacian
lap_x = 0.; lap_y = -0.; lap_z = f(lap_x, lap_y);
dzx = -4*lap_x; dzy = -4*lap_y;
quiver3(lap_x, lap_y, lap_z, 0.15, 0, 0.15*dzx, 'Color', [0 0 1], 'LineWidth', 2, 'MaxHeadSize', 5)
quiver3(lap_x, lap_y, lap_z, -0.15, 0, -0.15*dzx, 'Color', [0 0 1], 'LineWidth', 2, 'MaxHeadSize', 5)
quiver3(lap_x, lap_y, lap_z, 0, 0.15, 0.15*dzy, 'Color', [0 0 1], 'LineWidth', 2, 'MaxHeadSize', 5)
quiver3(lap_x, lap_y, lap_z, 0, -0.15, -0.15*dzy, 'Color', [0 0 1], 'LineWidth', 2, 'MaxHeadSize', 5)

% Set limits
zlim([-2.7, Inf])

% Export the plot as schematicfieldmass.svg
set(gcf, 'Renderer', 'painters');
print(gcf, '1_schematicfieldmass.svg', '-dsvg');
close(fig);