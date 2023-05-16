function res = bfftn(x)

res = 1/sqrt(numel(x))*fftn(x);

