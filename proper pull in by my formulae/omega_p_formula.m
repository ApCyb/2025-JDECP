function omega_p_formula = omega_p_formula(tau_1, tau_2, k, K_vco)

mu = pi*k - 1;
K_vco_ht = 1 / (k*(2*tau_1 + tau_2 + 2*sqrt(tau_1*(tau_1 + tau_2))))

if(tau_2 == 0)
    if(K_vco <= K_vco_ht)
    omega_p_formula = K_vco;
    else
    omega_p_formula = omega_sep_formula(tau_1, tau_2, k, K_vco);
    end
else
K_vco_pt = max(K_vco_ht, mu/(k*tau_2))

if(K_vco <= K_vco_ht)
    omega_p_formula = K_vco;
else 
    if(K_vco > K_vco_ht && K_vco <= K_vco_pt)
        omega_p_formula = omega_sep_formula(tau_1, tau_2, k, K_vco);
    else
        omega_p_formula = omega_ss_formula(tau_1, tau_2, k, K_vco);
    end
end
end
