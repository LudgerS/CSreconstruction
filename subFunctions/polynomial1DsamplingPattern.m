function mask = polynomial1DsamplingPattern(dim, usFactor, centerFraction, degree)
% polynomial undersampling mask adapted from
% Zijlstra F, Viergever MA, Seevinck PR. "Evaluation of variable density and 
% data-driven k-space undersampling for compressed sensing magnetic 
% resonance imaging". Invest Radiol. 2016;51:410–419.
%
% performs undersampling in 1 dimension, e.g. in the phase direction of a 
% 2D MR scan

nx = dim(1);

if usFactor == 1
    mask = true(dim);
else
       
nSamples = round(nx*usFactor);
nCenter = round(nSamples*centerFraction);

centerMask = false(nx, 1);
centerMask((round((nx + 1)/2 - (nCenter + 1)/2) + 1):(round((nx + 1)/2 + (nCenter + 1)/2) - 1)) = true;

weights = abs((1:nx) - (nx + 1)/2);
weights = (1 - weights/max(weights)).^degree;
weights = weights(~centerMask);

candidateLines = 1:nx;
candidateLines = candidateLines(~centerMask);

samples = datasample(candidateLines, nSamples - nCenter, 'replace', false, 'weights', weights);

mask = false(dim);
mask(repmat(centerMask, [1, dim(2)])) = true;
mask(samples, :) = true;


end

