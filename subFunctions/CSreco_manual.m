function reco = CSreco_manual(kspData, lambda, wPhi, wTV)
% performs a compressed sensing reconstruction with user specified
% regulatization strength lambda
%   kspData - k-space data with low frequencies in center and no leading
%   empty dimensions (suqeeze(kspData) == kspData).
%
%   lambda - regulatization strength. Reasonable values depend on the noise
%   level. Typical values would be between 10^-5 and 10^-1.
%
%   wPhi and wTV govern the weighting of the image l1-norm and isotropic
%   Total Variation regulatization. wPhi = 1, wTV = 0 -> only image
%   l1-norm; wPhi = 0, wTV = 1 -> only TV; wPhi = 1, wTV = 1 -> equal
%   weighting
%   For a constant c, a reconstruction with [c*lambda, wPhi, wTV] would
%   yield identical results as one with [lambda, c*wPhi, c*wTV].

%% normalization
Max = max(abs(kspData(:)));
kspData = kspData/Max;


%% change of Fourier transform convention, find undersampling mask
kspData = bfftn(cifftn(kspData));
mask = abs(kspData) > 10^(-12);

fprintf('Acceleration %.2f\n', numel(mask)/sum(mask(:)))


%% determine data dimension
dim = size(kspData);

if length(dim) == 2
    dimSwitch = '2D';
elseif length(dim) == 3
    dimSwitch = '3D';
else
    error('unsupported dimension')
end


%% parameters
mu = 1;         % ADMM step size
admmIter = 50;   


%% Precalculation
switch dimSwitch
    case '2D'
        laplaceColumn = laplaceColumn2DnonIso(dim(1), dim(2));

    case '3D'
        laplaceColumn = laplaceColumn3DnonIso(dim(1), dim(2), dim(3));

    otherwise
        error('unsupported dimension')
end

eVs = sqrt(prod(dim))*bifftn(reshape(laplaceColumn, dim));
inverseO = 2*mask + wPhi^2*1 + wTV^2*eVs; 
inverseO = inverseO.^-1;


%% Reco function
switch dimSwitch
    case '2D'
        recoFunc = @(lambda) ADMM2D(kspData, inverseO, admmIter, mu, lambda, wPhi, wTV);

    case '3D'
        recoFunc = @(lambda) ADMM3D(kspData, inverseO, admmIter, mu, lambda, wPhi, wTV);

    otherwise
        error('unsupported dimension')
end


%% Compute reconstruction
reco = recoFunc(lambda);
reco = reco*Max;

fprintf('reconstruction completed\n')
