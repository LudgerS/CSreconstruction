function res = cfftn(x)

res = 1/sqrt(numel(x))*fftshift(fftn(ifftshift(x)));

