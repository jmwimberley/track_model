% % Determine the displacement within the trajectory

% DPX = xt_in(2:end) - xt_in(1:end-1) + ...
%     (alpha(2:end)).*(xt_out(2:end) - xt_in(2:end)) - ...
%     (alpha(1:end-1)).*(xt_out(1:end-1) - xt_in(1:end-1));
% DPY = yt_in(2:end) - yt_in(1:end-1) + ...
%     (alpha(2:end)).*(yt_out(2:end) - yt_in(2:end)) - ...
%     (alpha(1:end-1)).*(yt_out(1:end-1) - yt_in(1:end-1));
% 
% % Determine the vectors for the position
% PX = xt_in + alpha.*(xt_out - xt_in);
% PY = yt_in + alpha.*(yt_out - yt_in);

% 
% % Total distance
% S = sum(sqrt(DPX*transpose(DPX) + DPY*transpose(DPY)));
%Minimise S for shortest distance

%% TRACKLINE CURVATURE 
% 
% % Distance increment between points
% ds = S_out(2:end) - S_out(1:end-1);
% 
% % Curvilinear abscissa increments
% bx_j = (Rx(2:end) - Rx(1:end-1))./ds;
% by_j = (Ry(2:end) - Ry(1:end-1))./ds;
% 
% % Normalised curvilinear abscissa
% % t = sqrt(bx_j.^2 + by_j.^2);
% 
% % Trajectory matrix
% n = length(ds);          
% v = ones(n,1)*2*(ds(1) + ds(1));    
% u = ones(n-1,1)*ds(1); 
% Hs_v = diag(v) + diag(u,1) + diag(u,-1); 
% 
% % Apply the boundary conditions
% bx_v = zeros(1,length(bx_j));
% by_v = zeros(1,length(by_j));
% bx_v(2:end) = 6*(bx_j(2:end) - bx_j(1:end-1));
% by_v(2:end) = 6*(by_j(2:end) - by_j(1:end-1));
% 
% % Curvilinear abscissa
% bx = zeros(length(bx_v)+2,1);
% by = zeros(length(by_v)+2,1);
% bx(2:end-1) = bx_v;
% by(2:end-1) = by_v;
% % Resulting trajectory matrix
% Hs = zeros(length(bx),length(bx));
% Hs(1,1) = 1;
% Hs(end,end) = 1;
% Hs(2,1) = ds(1);
% Hs(end-1,end) = ds(1);
% Hs(2:end-1,2:end-1) = Hs_v;
% 
% % Solve for curvature of the splines
% zx = Hs\bx;
% zy = Hs\by;
% 
% % Get the curvature data
% K = zx.^2 + zy.^2;
% 
% % Objective function
% f = sum(sqrt(transpose(zx)*zx + transpose(zy)*zy));