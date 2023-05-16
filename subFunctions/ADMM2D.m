function r = ADMM2D(kspData, iO, maxIter, mu, lambda, wPhi, wTV)
% CS reconstruction using the alternating direction method of multipliers:
% Goldstein T, O’Donoghue B, Setzer S, Baraniuk R. Fast alternating
% direction optimization methods. SIAM J Imaging Sci. 2014;7:1588–1623.
%
% Uses a fast Fourier transform–based exact inversion exploiting the 
% circulant matrix structure arising from periodic boundary conditions.
% For details consider the documentation in the repositiory and
% Chen M. On the solution of circulant linear systems. SIAM J Numer
% Anal. 1987;24:668–683.

r = bifftn(kspData);
d = dTrafoIso2D(r, wPhi ,wTV);
k = zeros(size(d));
a = 1;

rhs = @(d,k) 2*bifftn(kspData) + adTrafoIso2D(mu*d + k,wPhi,wTV);

for jj = 1:maxIter

    d = dUpdateIso2D_w(r,k,lambda,wPhi,wTV,mu);
 
    rNew = bifftn(iO.*bfftn(rhs(d,k)));
    
    kNew = k + mu*(d - dTrafoIso2D(rNew,wPhi,wTV));
        
    aNew = (1 + sqrt(1 + 4*a^2))/2;
    rNew = rNew + (a-1)/aNew*(rNew-r);
    kNew = kNew + (a-1)/aNew*(kNew-k);
    
    r = rNew;
    k = kNew;
    a = aNew;
 
end    


end

