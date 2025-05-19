function [phi_new] = wave_eq_2D(phi,input,gamma)
%WAVE_EQ Simulate the wave equation of model using finite difference scheme
%   Wave equation can be written as:
%       phi_n 
%           + (2/gamma)*(phi_(n+1)-phi_(n-1))/(2*dt) 
%           + (1/gamma^2)*(phi_(n+1)+phi_(n-1)-2*phi_n)/(dt^2)
%       = input = r^2*Laplacian(phi) + C(phi) + nu0*phi + P
%   Therefore:
%       phi_(n+1)
%       = [input - (1/gamma^2-1/gamma)*phi_(n-1) - (1-2/gamma^2)*phi_n]/[1/gamma + 1/gamma^2]
%   List of activity inputs (phi)
	%	phi: Nx x Ny x 2 array containing phi(r, t=(n-1)dt) and phi(r, t = ndt)
%	List of driving inputs (input)
	%	input: Nx x Ny array containing r^2*Laplacian(phi)(r, t = ndt) + C(phi)(r, t = ndt) + nu0*phi(r, t = ndt) + P(r, t = ndt)
%	List of homogeneous parameters (gamma)
	%	gamma = damping rate
%	List of outputs (phi_new)
	%	phi_new = Nx x Ny array contianin phi(r, t = (n+1)dt)

phi_new = (input - (1/gamma^2 - 1/gamma)*phi(:, :, end-1) - (1 - 2/gamma^2)*phi(:, :, end)) / ...
    (1/gamma + 1/gamma^2);
end