function out = dUpdateIso2D_w(r, k, lambda, wPhi, wTV, mu)
% for details, see documentaion provided under 
% https://github.com/LudgerS/CSreconstruction  

nx = size(r,1);

out = dTrafoIso2D(r, wPhi, wTV);

RHSx = out((nx + 1):(2*nx),:) - k((nx + 1):(2*nx),:)/mu;
RHSy = out((2*nx + 1):(3*nx),:) - k((2*nx + 1):(3*nx),:)/mu;
s = sqrt(abs(RHSx).^2 + abs(RHSy).^2);

temp = max(s - lambda/mu, 0)./(s + 10^(-15));

out(1:nx,:) = shrink(out(1:nx,:) - k(1:nx,:)/mu, lambda/mu);

out((nx + 1):(2*nx),:) = RHSx.*temp;
out((2*nx + 1):(3*nx),:) = RHSy.*temp;


end