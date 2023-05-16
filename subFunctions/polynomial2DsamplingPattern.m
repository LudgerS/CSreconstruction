function mask = polynomial2DsamplingPattern(dim, usFactor, centerFraction, degree)
% polynomial undersampling mask adapted from
% Zijlstra F, Viergever MA, Seevinck PR. "Evaluation of variable density and 
% data-driven k-space undersampling for compressed sensing magnetic 
% resonance imaging". Invest Radiol. 2016;51:410–419.
%
% performs undersampling in 2 dimensions, e.g. in the phase plane of a 3D
% MR scan

nx = dim(1);
ny = dim(2);

if usFactor == 1
    mask = true(dim);
elseif usFactor > pi/4
    error('Unsupported usFactor')
elseif centerFraction == 1
    mask = createCircle(dim, (dim + 1)/2, sqrt(centerFraction*usFactor*nx*ny/pi));
else

positions = reshape(1:nx*ny, dim);
mask = createCircle(dim, (dim + 1)/2, sqrt(centerFraction*usFactor*nx*ny/pi));

[x, y] = meshgrid(-((ny + 1)/2 - 1):(ny - (ny + 1)/2), -((nx + 1)/2-1):(nx - (nx + 1)/2));
x = x/max(x(:));
y = y/max(y(:));

weights = sqrt(x.^2 + y.^2);
weights = (max(1 - weights, 0)).^degree;
weights(mask) = 0;

samples = datasample(positions(:), round(nx*ny*usFactor) - sum(mask(:)), 'replace', false, 'weights', weights(:));
mask(samples) = true;


end

