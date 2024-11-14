function omega_ss_formula = omega_ss_formula(tau_1, tau_2, k, K_vco)

accuracy = 0.00001;

mu = pi*k - 1;
K_vco_minus = 1/(k*(2*tau_1 + tau_2 + 2*sqrt(tau_1*(tau_1 + tau_2))));
K_vco_eta_0 = mu/(k*tau_2);

if(K_vco > max(K_vco_eta_0, K_vco_minus))
    xi = (k*tau_2*K_vco + 1)  /(2*sqrt(K_vco*(tau_1 + tau_2)));
    eta = (k*tau_2*K_vco - mu)/(2*sqrt(K_vco*(tau_1 + tau_2)));
    rho = sqrt(abs(xi^2 - k));
    kappa = sqrt(eta^2 + k*mu);

    r_1 = ((k*sqrt(tau_2*K_vco) + eta)^2 - kappa^2) / ((k*sqrt(tau_2*K_vco) - eta)^2 - kappa^2) *...
    ( ((k*sqrt(tau_2*K_vco) + kappa)^2 - eta^2) / ((k*sqrt(tau_2*K_vco) - kappa)^2 - eta^2))^(eta/kappa);
    r_2 = ((k*sqrt(tau_2*K_vco) + xi)^2 - rho^2) / ((k*sqrt(tau_2*K_vco) - xi)^2 - rho^2) *...
    ( ((k*sqrt(tau_2*K_vco) + rho)^2 - xi^2) / ((k*sqrt(tau_2*K_vco) - rho)^2 - xi^2))^(xi/rho) * (xi > sqrt(k)) +...
    ((k*sqrt(tau_2*K_vco) + sqrt(k)) ./ (k*sqrt(tau_2*K_vco) + sqrt(k)) *...
    exp(2*sqrt(k*tau_2*K_vco) / (k*tau_2*K_vco - 1)))^2 *(xi == sqrt(k)) +...
    (k*tau_2*K_vco + 2*xi*sqrt(tau_2*K_vco) + 1) /...
            (k*tau_2*K_vco - 2*xi*sqrt(tau_2*K_vco) + 1) *...
    exp(2*xi/rho*(atan((1 - k*tau_2*K_vco) / (2*rho*sqrt(tau_2*K_vco))) + pi/2)) * (xi < sqrt(k));
    r = max(r_1, r_2);


        omega_start = (sqrt(r) - 1)/(sqrt(r) + 1) * K_vco;
        omega_finish = omega_sep_formula(tau_1, tau_2, k, K_vco);
        if(omega_finish - omega_start < 0.00239/(K_vco_eta_0))
            omega_ss_formula = omega_finish;
        else
            omega_interval = omega_ss_intervals(tau_1, tau_2, k, K_vco, omega_start, omega_finish);
            while(omega_finish - omega_start > accuracy)
                omega_start = omega_interval(1);
                omega_finish = omega_interval(2);
                omega_interval = omega_ss_intervals(tau_1, tau_2, k, K_vco, omega_start, omega_finish);
            end
            omega_ss_formula = omega_interval(1);
        end
end
end


function omega_ss_intervals = omega_ss_intervals(tau_1, tau_2, k, K_vco, omega_start, omega_finish)

omega_step = - (omega_finish - omega_start)/5;

for omega_e_free = omega_finish:omega_step:omega_start
    res = curve1curve3intersection(tau_1, tau_2, k, K_vco, omega_e_free);

    if (curve2(tau_1, tau_2, k, K_vco, omega_e_free, res.y1, res.y0) < 0)
        if omega_e_free == omega_finish
            omega_ss_intervals = [omega_finish-abs(omega_step) omega_finish];
        else
            omega_ss_intervals = [omega_e_free omega_e_free+abs(omega_step)];
        end
        break;
    end
end
end


function curve1curve3intersection = curve1curve3intersection(tau_1, tau_2, k, K_vco, omega_e_free)
mu = pi*k - 1;
eta = (k*tau_2*K_vco - mu)/(2*sqrt(K_vco*(tau_1 + tau_2)));
kappa = sqrt(eta^2 + k*mu);

y0 = @(y1) ((1/k - omega_e_free/(k*K_vco)) * (y1 - k*tau_2*K_vco/sqrt((tau_1 + tau_2)*K_vco)*(1/k + omega_e_free/(k*K_vco)) )./...
    ( (1/k + omega_e_free/(k*K_vco)) - 1/k/sqrt((tau_1 + tau_2)*K_vco)*y1));

curve1 = @(y1) ((y1 + (kappa - eta)*(1/k + omega_e_free/(k*K_vco)))./...
            (y0(y1) - (kappa - eta)*(1/k - omega_e_free/(k*K_vco)))).^((kappa - eta)/kappa) ...
    - ...
((y0(y1) + (kappa + eta)*(1/k - omega_e_free/(k*K_vco)))./...
     (y1 - (kappa + eta)*(1/k + omega_e_free/(k*K_vco)))).^((kappa + eta)/kappa);


y1_sep = (kappa + eta)*(1/k + omega_e_free/(k*K_vco));
y1_left = y1_sep+0.000000000001;
y1_sol = fzero(curve1, [y1_left, k*sqrt(tau_2*K_vco)*(1/k + omega_e_free/(k*K_vco))]);
y0_sol = (1/k - omega_e_free/(k*K_vco)) * (y1_sol - k*tau_2*K_vco/sqrt((tau_1 + tau_2)*K_vco)*(1/k + omega_e_free/(k*K_vco)) )/...
    ( (1/k + omega_e_free/(k*K_vco)) - 1/k/sqrt((tau_1 + tau_2)*K_vco)*y1_sol);
curve1curve3intersection.y1 = y1_sol;
curve1curve3intersection.y0 = y0_sol;
end



function curve2 = curve2(tau_1, tau_2, k, K_vco, omega_e_free, y1, y0)
xi = (k*tau_2*K_vco + 1)  /(2*sqrt(K_vco*(tau_1 + tau_2)));
rho = sqrt(abs(xi^2 - k));

if(xi < sqrt(k))
        curve2 = (y1^2 - 2*xi*(1/k + omega_e_free/(k*K_vco))*y1 + k*(1/k + omega_e_free/(k*K_vco))^2)/...
                 (y0^2 + 2*xi*(1/k - omega_e_free/(k*K_vco))*y0 + k*(1/k - omega_e_free/(k*K_vco))^2) - ...
        exp(2*xi/rho*(...
        atan(((1/k - omega_e_free/(k*K_vco))*rho)/(y0 + xi*(1/k - omega_e_free/(k*K_vco)))) - ...
        atan((y1 - xi*(1/k + omega_e_free/(k*K_vco)))/((1/k + omega_e_free/(k*K_vco))*rho)) + ...
        pi/2));
end
    if(xi == sqrt(k))
        curve2 = (y1 - sqrt(k)*(1/k + omega_e_free/(k*K_vco)))/...
                 (y0 + sqrt(k)*(1/k - omega_e_free/(k*K_vco))) - ...
        exp(sqrt(k)*(1/k - omega_e_free/(k*K_vco))/(y0 + sqrt(k)*(1/k - omega_e_free/(k*K_vco))) + ...
            sqrt(k)*(1/k + omega_e_free/(k*K_vco))/(y1 - sqrt(k)*(1/k + omega_e_free/(k*K_vco))));
    end

    if(xi > sqrt(k))
        curve2 = ((y1 + (rho - xi)*(1/k + omega_e_free/(k*K_vco)))/...
                  (y0 - (rho - xi)*(1/k - omega_e_free/(k*K_vco))))^((rho - xi)/rho) - ...
                 ((y0 + (rho + xi)*(1/k - omega_e_free/(k*K_vco)))/...
                  (y1 - (rho + xi)*(1/k + omega_e_free/(k*K_vco))))^((rho + xi)/rho);
    end

end
