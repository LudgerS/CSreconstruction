function res = cifftn(x)

res = sqrt(numel(x))*fftshift(ifftn(ifftshift(x)));

