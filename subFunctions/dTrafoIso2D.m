function [ out ] = dTrafoIso2D(image, wPhi, wTV)
% stacked and weighted image l1-norm and discrete gradient in 2D
% [l1; dx; dy]

[nx, ny] = size(image);

out = zeros(3*nx, ny);

out(1:nx,:) = wPhi*image;
out((nx+1):(2*nx),:) = wTV*([image(:,2:end), image(:,1)] - image);
out((2*nx+1):(3*nx),:) = wTV*([image(2:end,:); image(1,:)] - image);


end

