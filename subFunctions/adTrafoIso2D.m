function [ out ] = adTrafoIso2D(in, wPhi, wTV)
% the adjoint transform of dTrafoIso2D

[nx, ~] = size(in);
nx = nx/3;

out = wPhi*in(1:nx,:)...
            + wTV*([in(nx+1:2*nx,end), in(nx+1:2*nx,1:end-1)] - in(nx+1:2*nx,:))...
            + wTV*([in(3*nx,:); in(2*nx+1:3*nx-1,:)] - in(2*nx+1:3*nx,:));
        
end

