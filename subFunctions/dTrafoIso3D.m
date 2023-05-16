function out = dTrafoIso3D(image, wPhi, wTV)
% stacked and weighted image l1-norm and discrete gradient in 3D
% [l1; dx; dy; dz]

[nx, ny, nz] = size(image);

out = (1 + 1i)*ones(4*nx, ny, nz);
out(1:nx,:,:) = wPhi*image;

out((nx + 1):(2*nx),:,:) = wTV*([image(:,2:end,:), image(:,1,:)] - image);
out((2*nx + 1):(3*nx),:,:) = wTV*([image(2:end,:,:); image(1,:,:)] - image);
out((3*nx + 1):(4*nx),:,:) = wTV*(cat(3,image(:,:,2:end), image(:,:,1)) - image);



end

