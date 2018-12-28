function value = hermite(n,x)
%base cases
if n == -1
    value = zeros(size(x));
    return;
end
if n == 0
    value = ones(size(x));
    return;
end

%recursion
value = x.*hermite(n-1,x)-(n-1)*hermite(n-2,x);
end