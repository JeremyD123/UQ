function Pi = calculateStrainEnergy(lambda,mu,N)

Pi = (2*sum(lambda) + 2*sum(mu)) * 0.01^2 * 100 * 100 / N; 

fprintf('total strain energy = %.4f\n', Pi);

end