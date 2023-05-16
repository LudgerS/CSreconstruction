function [out] = shrink(in,epsilon)

out = in./abs(in+10^-15).*max(abs(in)-epsilon,0);

end

