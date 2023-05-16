function out = laplaceColumn2DnonIso(nx,ny)
% first column of the matrix
% dx'*dx + dy'*dy
% where dx, dy are the discrete gradient matrices in x, y direction
% with periodic boundary conditions.

out = zeros(nx*ny, 1);

out(1) = 4;
out([2, nx, nx + 1, nx*ny - nx + 1]) = -1;

end