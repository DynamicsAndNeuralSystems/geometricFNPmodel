function L = Laplacian_2D(phi)
%LAPLACIAN Calculate the Discrete Laplacian of the current phi variable
%   For each (i,j) element of phi:
%   Laplacian(phi(i)) =
%   (phi(i+1,j)+phi(i-1,j)+phi(i,j+1),phi(i,j-1)-4*phi(i,j))/dx^2

%L = zeros(size(phi));
% L(2:end-1,2:end-1) = (phi(3:end,2:end-1) + phi(1:end-2,2:end-1) + ...
%     phi(2:end-1,3:end) + phi(2:end-1,1:end-2) - ...
%     4*phi(2:end-1,2:end-1));
% L(1,2:end-1) = phi(1,1:end-2) + phi(2,2:end-1) + phi(1,3:end) + phi(end,2:end-1) - 4*phi(1,2:end-1);
% L(end,2:end-1) = phi(end,1:end-2) + phi(end-1,2:end-1) + phi(end,3:end) + phi(1,2:end-1) - 4*phi(end,2:end-1); 
% L(2:end-1,1) = phi(1:end-2,1) + phi(2:end-1,2) + phi(3:end,1) + phi(2:end-1,end) - 4*phi(2:end-1,1);
% L(2:end-1,end) = phi(1:end-2,end) + phi(2:end-1,end-1) + phi(3:end,end) + phi(2:end-1,1) - 4*phi(2:end-1,end);
% L(1,1) = phi(1,2) + phi(2,1) + phi(end,1) + phi(1,end) - 4*phi(1,1);
% L(1,end) = phi(1,end-1) + phi(2,end) + phi(end,end) + phi(1,1) - 4*phi(1,end);
% L(end,1) = phi(end,2) + phi(end-1,1) + phi(1,1) + phi(end,end) - 4*phi(end,1);
% L(end,end) = phi(end,end-1) + phi(end-1,end) + phi(1,end) + phi(end,1) - 4*phi(end,end);

L = phi([2:end 1], :) + phi([end 1:end-1], :) + phi(:, [2:end 1]) + phi(:, [end 1:end-1]) - 4*phi;


end

