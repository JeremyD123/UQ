function d = trid(x)
d = length(x);
y = sum((x-1).^2) - sum(x(2:d).*x(1:d-1));
end