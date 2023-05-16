function out = dUpdateIso3D_w(r, k, lambda, wPhi, wTV, mu)
% for details, see documentaion provided under 
% https://github.com/LudgerS/CSreconstruction  

nx = size(r,1);

out = dTrafoIso3D(r, wPhi, wTV);

RHSx = out((nx + 1):(2*nx),:,:) - k((nx + 1):(2*nx),:,:)/mu;
RHSy = out((2*nx + 1):(3*nx),:,:) - k((2*nx + 1):(3*nx),:,:)/mu;
RHSz = out((3*nx + 1):(4*nx),:,:) - k((3*nx + 1):(4*nx),:,:)/mu;

s = sqrt(abs(RHSx).^2 + abs(RHSy).^2 + abs(RHSz).^2);

temp = max(s - lambda/mu, 0)./(s + 10^(-15));

out(1:nx,:) = shrink(out(1:nx,:) - k(1:nx,:)/mu, lambda/mu);

out((nx + 1):(2*nx),:,:) = RHSx.*temp;
out((2*nx + 1):(3*nx),:,:) = RHSy.*temp;
out((3*nx + 1):(4*nx),:,:) = RHSz.*temp;


end