function p = fastnormcdf(x,mu,sigma)
z = (x - mu)/sigma;
p = 0.5 * erfc(-z / sqrt(2));
end

