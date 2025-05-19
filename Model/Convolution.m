function Cphi = Convolution(hetparam,phi,dt)
%CONVOLUTION Calculates the Convolution term
%   for the jth shortcut, C is calculated as
%   C(r,t) = sum_(all points r') (cj/|Rji||Rjf|) I(r in Rjf) I(r' in Rji) phi(r',t - tauj) dx^2
%       - (cj/|Rji||Rjf|) I(r' in Rjf) I(r in Rji) phi(r,t) dx^2
%   We can calculate this using the following steps
%       1. Calculate normalized matrices Ranorm_m = K(r-a_m), Rbnorm_m = K(r-b_m);
%       2. Calculate average phi at a: avgphiRa_m = sum(phi(r')*K(r'-a_m)); avgphiRadelay_m = sum(phi(r',t - tau_m)*K(r'-a_m))
%       3. Then C = sum_m c_m*(Rbnorm.*avgphiRadelay - Ranorm.*phi)

% Compute avgphiRi
avgphiRa = zeros(1, 1, hetparam.m);
phi_end = phi(:, :, end);
Ra_reshaped = reshape(hetparam.Ra, [], hetparam.m);
avgphiRa = phi_end(:)' * Ra_reshaped;
if ~all(hetparam.tau == 0)
    avgphiRadelay = zeros(1, 1, hetparam.m);
    for j = 1:hetparam.m
        avgphiRadelay(j) = sum(phi(:, :, end - ceil(hetparam.tau(j)/dt)).*hetparam.Ra(:, :, j), 'all');
    end
end
c = reshape(hetparam.c, [1, 1, hetparam.m]);
avgphiRa = reshape(avgphiRa, [1, 1, hetparam.m]);
if all(hetparam.tau == 0)
    Cphi = (hetparam.Rb - hetparam.Ra) .* avgphiRa;
else
    Cphi = hetparam.Rb .* avgphiRadelay - hetparam.Ra .* avgphiRa;
end
Cphi = sum(c .* Cphi, 3);



end

