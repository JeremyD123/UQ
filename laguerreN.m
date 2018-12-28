function value = laguerreN(n,alpha,x)
    value = 0;
    for j = 0:n
        coef = gamma(n+alpha+1)/gamma(n-j+1)/gamma(alpha+j+1);
        value = value+(-1)^j/factorial(j)*coef.*(x.^j);
    end
    value = sqrt(factorial(n)*gamma(alpha+1)/gamma(n+alpha+1))*value;
end