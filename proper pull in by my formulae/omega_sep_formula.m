function omega_p_separatrix = pull_in_separatrix_formula(tau_1, tau_2, k, K_vco)
%PULL_IN_SEPARATRIX_FORMULA calculates pull-in frequency by formula as a
%separatrix cycle

mu = pi*k - 1;
xi = (1 + k*tau_2*K_vco)/(2*sqrt((tau_1+tau_2)*K_vco));
eta = (k*tau_2*K_vco - mu)/(2*sqrt((tau_1+tau_2)*K_vco));
rho = sqrt(abs(xi^2 - k));
kappa = sqrt(eta^2 + k*mu);
K_vco_ht = 1 / (k*(2*tau_1 + tau_2 + 2*sqrt(tau_1*(tau_1 + tau_2))));
K_vco_fn = (2*tau_1 + tau_2 + 2*sqrt(tau_1*(tau_1 + tau_2)))/(k*tau_2^2);

% Separatrix cycle exists for K_vco > K_vco_minus
if(K_vco > K_vco_ht)
    if(K_vco < K_vco_fn)
        s_1 = ((kappa - eta)^2 + 2*xi*(kappa - eta) + k)/((kappa + eta)^2 - 2*xi*(kappa + eta) + k)*...
            exp(2*xi/rho*(atan(((xi - eta)^2 + rho^2 - kappa^2)/(2*rho*kappa)) + pi/2));
    end

    if(K_vco == K_vco_fn)
        s_1 = (kappa - eta + sqrt(k))/((kappa + eta - sqrt(k))*...
        exp(2*sqrt(k)*kappa/(kappa^2 - (eta - sqrt(k))^2)))^2;
    end

    if(K_vco > K_vco_fn)
        s_1 = ((kappa - eta + xi)^2 - rho^2)/((kappa + eta - xi)^2 - rho^2)*...
        ((((kappa + rho)^2) - (xi - eta)^2)/((kappa - rho)^2 - (xi - eta)^2))^(xi/rho);
    end

    omega_p_separatrix = (sqrt(s_1) - 1)/(sqrt(s_1) + 1)*K_vco;
end
end

