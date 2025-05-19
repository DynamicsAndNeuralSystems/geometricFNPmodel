function ts = run_periodic(topology,homparam,hetparam,stim)
%RUN with ZERO TRIVIAL CONDITIONS: phi(r,0) = phi_t(r,0) = 0;

% Takes in input of model parameters, timesteps and lengthsteps.
% Returns timeseries at each point.
%   List of step-size inputs (topology)
    %   T: time domain [0, T]
    %   L: space domain [0, L] x [0, L]
    %   Nt: number of timesteps
    %   Nx: number of grid points
%   List of homogeneous model parameters inputs (homparam)
    %   gamma: timescale/decay - check that gamma > 0
    %   r: characteristic axonal length 
        %   check r > 0 && r*gamma < dx/dt (stability)
    %   nu0: dissipative rate
        %   check that 0 <= nu0 < 1
%   List of heterogeneous model parameters (hetparam)
    %   G: global coupling parameter
        %   check G > 0
    %   m: number of unidirectional pipe
        %   check m > 0
    %   cj = c1,..,cm: connectivity strength of jth pipe
        %   check cj > 0
    %   tauj = tau1,...,taum: travel time of jth current (as multiples of dt)
        %   check tauj < |Rji-Rjf|/(r*gamma)
    %   Rji = R1i,...,Rmi: regions at which pipe starts
        %   check Rjis are in 1:N1x1:N2
    %   Rjf = R1f,...,Rmf: regions at which pipe ends
        %   check Rjfs are in 1:N1x1:N2
        %   check mutually exclusive: does not exist any j1,j2 with
        %   overlapping Rji and overlapping Rjf
%   List of stimulation inputs (stim)
    %   stimnum: number of times of stimulation
    %   stimtk = stimt1,...,stimtn = times of stimulation - check less than Nt
    %   stimIk = stimI1,...,stimIn = stimulation intensity
    %   stimRk = stimR1,...,stimRn = stimulation positions
    %   sigma = [sigma_x, sigma_t] = Gaussian approximation of delta impulse
% Outputs
    %   ts: Nx x Ny x Nt timeseries array

dx = topology.L / topology.Nx;
dt = topology.T / topology.Nt;

if (homparam.r * homparam.gamma * sqrt(2) > dx / dt)
    error('Error, Unstable choice of dx and dt')
end

% Initiate output timeseries array  
ts = zeros(topology.Nx, topology.Nx, topology.Nt);
% Initiate temporary recursive timeseries array
% Need at least tauj previous timesteps for convolution 
% and at least 1 previous timestep for wave equation
% Therefore we track 1 + max(tau, 1) timesteps for each iteration
if hetparam.m > 0
    phi = zeros(topology.Nx, topology.Nx, 1 + max([ceil(hetparam.tau/dt), 1]));
else
    phi = zeros(topology.Nx, topology.Nx, 2);
end

if hetparam.m > 0

    % Change hetparam.Ri and hetparam.Rf to 3D array if they are 2D
    
    Ra = zeros(topology.Nx, topology.Nx, hetparam.m);
    Rb = zeros(topology.Nx, topology.Nx, hetparam.m);
    
    for k = 1:hetparam.m
        a_x = hetparam.a(1, k); a_y = hetparam.a(2, k);
    
        for i = 1:topology.Nx
            for j = 1:topology.Nx
                distx = abs(i*dx - a_x);
                if distx > topology.L / 2
                    distx = topology.L - distx;
                end
                disty = abs(j*dx - a_y);
                if disty > topology.L / 2
                    disty = topology.L - disty;
                end
                Ra(i, j, k) = exp(-0.5*(distx^2 + disty^2) / (hetparam.sigmaeps)^2);
            end
        end
    end
    
    for k = 1:hetparam.m
        b_x = hetparam.b(1, k); b_y = hetparam.b(2, k);
    
        for i = 1:topology.Nx
            for j = 1:topology.Nx
                distx = abs(i*dx - b_x);
                if distx > topology.L / 2
                    distx = topology.L - distx;
                end
                disty = abs(j*dx - b_y);
                if disty > topology.L / 2
                    disty = topology.L - disty;
                end
                Rb(i, j, k) = exp(-0.5*(distx^2 + disty^2) / (hetparam.sigmaeps)^2);
            end
        end
    end
    
    for k = 1:hetparam.m
        Ra(:, :, k) = Ra(:, :, k)/sum(Ra(:, :, k),'all');
        Rb(:, :, k) = Rb(:, :, k)/sum(Rb(:, :, k),'all');
    end
    
    hetparam.Ra = Ra;
    hetparam.Rb = Rb;
    
    % Turn hetparam.c into c / dx^2, so I dont have to keep multiplying
    hetparam.c = hetparam.c / dx^2;

end
    
% Create array of stimulation inputs by space

stimulationinputspace = zeros(topology.Nx, topology.Nx, stim.stimnum);
for n = 1:stim.stimnum
    i0 = stim.stimR(1, n) / dx; j0 = stim.stimR(2, n) / dx;
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
            stimulationinputspace(i, j, n) = exp(-0.5*(distx^2 + disty^2) * (dx^2) / (stim.sigma(1)^2));
        end
    end
end

stimulationinputspace = stimulationinputspace * 1/(2*pi * stim.sigma(1)^2);

% Create array of stimulation inputs by time

stimulationinputtime = zeros(topology.Nt, stim.stimnum);
for k = 1:stim.stimnum
    t0 = stim.stimt(k) / dt;
    for n = 1:topology.Nt
        distt = abs((n - 1) - t0);
        stimulationinputtime(n, k) = exp(-0.5 *(distt^2) * (dt^2) / (stim.sigma(2)^2));
    end
end

stimulationinputtime = stimulationinputtime * 1/(sqrt(2*pi) * stim.sigma(2));

% Normalise stimulationinputspace and stimulationinputtime so that weights
% add to 1/dx^2 and 1/dt respectively

for k = 1:stim.stimnum
    stimulationinputspace(:, :, k) = (1/dx^2) * ...
        stimulationinputspace(:, :, k) / sum(stimulationinputspace(:, :, k), "all");
    stimulationinputtime(:, k) = (1/dt) * ...
        stimulationinputtime(:, k) / sum(stimulationinputtime(:, k), 'all');
end


for count = 1:topology.Nt
    
    % Calculate input for wave equation: nu0*phi + r^2*nabla^2(phi) + C(phi) + P
    recurrentinput = homparam.nu0 * phi(:, :, end);
    laplacianinput = (homparam.r / dx)^2 * Laplacian_2D(phi(:, :, end));
    if hetparam.m > 0
        convolutioninput = Convolution(hetparam, phi, dt);
    else
        convolutioninput = 0;
    end
    stimulationinput = zeros(topology.Nx, topology.Nx);
    for k = 1:stim.stimnum
        stimulationinput = stimulationinput + ...
            stim.stimI(k) * stimulationinputspace(:, :, k) * stimulationinputtime(count, k);
    end
    input = recurrentinput + laplacianinput + convolutioninput + stimulationinput;
    if count == 1
        % Initial velocity conditions imply that phi_1 = phi_{-1}, hence...
        phinew = 0.5 * input * (homparam.gamma * dt)^2;
    else
        phinew = wave_eq_2D(phi, input, homparam.gamma * dt);
    end
    phi(:, :, 1:end-1) = phi(:, :, 2:end);
    phi(:, :, end) = phinew;
    ts(:, :, count) = phinew;
end


end

