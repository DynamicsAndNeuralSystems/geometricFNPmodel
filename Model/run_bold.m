function [ts_bold] = run_bold(topology, homparam, hetparam, stim, tol)
% Compute the BOLD response, which is approximated as the zero frequency component of the underlying stimulus evoked response
%   ts_bold(r) = sum_n=0^\infty phi(r, t = ndt).

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
    %   ts_bold: Nx x Ny array containing ts_bold(r)

if nargin < 5 || isempty(tol)
    tol = 1e-5;
end

dt = topology.T / topology.Nt;
dx = topology.L / topology.Nx;

iter = 0;

init = 0;

ts_bold = zeros(topology.Nx, topology.Nx);

delta = Inf;

while(delta > tol || isnan(delta))

    if (iter > 0)
        stim.stimnum = 0;
    else
        stim.stimnum = 1;
    end

    ts_het = run_init_periodic(topology, homparam, hetparam, stim, init);

    init = ts_het(:, :, end-1:end);

    ts_bold_prev = ts_bold;
    ts_bold = ts_bold + dt*squeeze(sum(ts_het, 3));

    if (iter > 0)
        delta = pdist2(ts_bold_prev(:)', ts_bold(:)', 'cosine');
    else
        delta = NaN;
    end

    iter = iter + 1;

end

ts_bold = ts_bold/sum(ts_bold, 'all') * (1/(dx^2 * (1-homparam.nu0)));

end