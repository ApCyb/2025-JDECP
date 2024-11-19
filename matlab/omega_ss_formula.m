function omega_pt = omega_ss_formula(k, K_vco, tau_1, tau_2)
%OMEGA_SS_FORMULA calculates frequency for the semistable cycle
%
% Parameters:
%   k       - slope
%   K_vco   - VCO gain
%   tau_1   - Time constant tau_1
%   tau_2   - Time constant tau_2
%
% Returns:
%   omega_pt - frequency corresponding to the birth of semistable cycle

% Step 1: Compute intermediate parameters
mu = pi*k - 1;
eta = (k*tau_2*K_vco - mu)/(2*sqrt(K_vco*(tau_1 + tau_2)));
kappa = sqrt(eta^2 + k*mu);

% Step 2: Determine the range for z1
z1_start = kappa+eta;
z1_finish = k*sqrt(tau_2*K_vco);

% Step 3: Solve for z1 using a numerical method
z1_pt = find_z1(k, K_vco, tau_1, tau_2, z1_start, z1_finish);

% Step 4: Calculate omega_pt based on the value of z1_pt
if(z1_pt == z1_start)
        % Case 1: z1_pt equals the starting point
    omega_pt = omega_sep_formula(k, K_vco, tau_1, tau_2);  % Use separatrix formula
else    % Case 2: General case where z1_pt does not equal the starting point
    z0_pt = z0(z1_pt, k, K_vco, tau_1, tau_2);    % Calculate corresponding z0
    s_pt = ((z0_pt + eta)^2 - kappa^2)/((z1_pt - eta)^2 - kappa^2)*...
((z0_pt + kappa + eta)*(z1_pt + kappa - eta)/...
    (z0_pt - kappa + eta)/(z1_pt - kappa - eta))^(eta/kappa);
    omega_pt = (sqrt(s_pt) - 1)/(sqrt(s_pt) + 1)*K_vco; % Calculate omega_pt
end
end


function found_z1 = find_z1(k, K_vco, tau_1, tau_2, z1_start, z1_finish)
%FIND_Z1 finds the value of z1 using a recursive or numerical method.
%
% Parameters:
%   k, K_vco, tau_1, tau_2 - System parameters
%   z1_start, z1_finish    - Search range for z1
%
% Returns:
%   found_z1 - Value of z1 that satisfies the main_curve condition

if(z1_finish - z1_start < 1e-12) % For convergence 
    found_z1 = z1_start;  % Return the starting value if range is negligible
else
    z1_middle = z1_start + (z1_finish-z1_start)/2; % Midpoint of the range
    if(main_curve(z1_middle, z0(z1_middle, k, K_vco, tau_1, tau_2), k, K_vco, tau_1, tau_2) < 0)
        % If main_curve is negative, recurse on the lower half
        found_z1 = find_z1(k, K_vco, tau_1, tau_2, z1_start, z1_middle);
    else
% Otherwise, find the root using a numerical solver (fzero)
        main_curve_z = @(z1) main_curve(z1, z0(z1, k, K_vco, tau_1, tau_2), k, K_vco, tau_1, tau_2);
        found_z1 = fzero(main_curve_z, [z1_middle, z1_finish]);
    end
end
end


function z0 = z0(z1, k, K_vco, tau_1, tau_2)
%Z0 calculates the value of z0 based on z1 and system parameters.
%
% Parameters:
%   z1, k, K_vco, tau_1, tau_2 - Input parameters
%
% Returns:
%   z0_val - Corresponding z0 value

mu = pi*k - 1;
xi = (k*tau_2*K_vco + 1)  /(2*sqrt(K_vco*(tau_1 + tau_2)));
eta = (k*tau_2*K_vco - mu)/(2*sqrt(K_vco*(tau_1 + tau_2)));
z0 = ((1+mu)*k*z1 - 2*(mu*xi + eta)*k)/((1+mu)*k - 2*(xi - eta)*z1);
end


function main_curve = main_curve(z1, z0, k, K_vco, tau_1, tau_2)
%MAIN_CURVE calculates the difference between the left-hand and right-hand sides of the main equation.
%
% Parameters:
%   z1, z0 - Variables for the main curve equation
%   k, K_vco, tau_1, tau_2 - System parameters
%
% Returns:
%   main_curve_val - Value of the main curve function

mu = pi*k - 1;
xi = (k*tau_2*K_vco + 1)  /(2*sqrt(K_vco*(tau_1 + tau_2)));
eta = (k*tau_2*K_vco - mu)/(2*sqrt(K_vco*(tau_1 + tau_2)));
rho = sqrt(abs(xi^2 - k));
kappa = sqrt(eta^2 + k*mu);

% Define thresholds for k and K_vco
k_fn =  2*(tau_1 + tau_2 - sqrt(tau_1*(tau_1 + tau_2))) / (pi*(2*tau_1 + tau_2 - 2*sqrt(tau_1*(tau_1 + tau_2))));
K_vco_fn = 1/(k*(2*tau_1 + tau_2 - 2*sqrt(tau_1*(tau_1 + tau_2))));

% Calculate the right-hand side of the equation
if(k < k_fn && K_vco < K_vco_fn)
    right_hand_side = (z0^2+2*xi*z0+k)/(z1^2-2*xi*z1+k)*exp(2*xi/rho*...
        (atan(rho/(z0 + xi)) - atan((z1 - xi)/rho) + pi/2));
else 
    if(k < k_fn && K_vco == K_vco_fn)
        right_hand_side = ((z0 + sqrt(k))/(z1 - sqrt(k))*...
        exp(sqrt(k)/(z0 + sqrt(k)) + sqrt(k)/(z1 - sqrt(k))))^2;
    else
        if(k < k_fn && K_vco > K_vco_fn || k >= k_fn)
            right_hand_side = (z0 + xi - rho)*(z0 + xi + rho)/(z1 - xi + rho)/(z1 - xi - rho)*...
    ((z0 + xi + rho)*(z1 + rho - xi)/(z0 + xi - rho)/(z1 - xi - rho))^(xi/rho);
        end
    end
end
% Calculate the left-hand side of the equation
left_hand_side = ((z0 + eta)^2 - kappa^2)/((z1 - eta)^2 - kappa^2)*...
    ((z0 + eta + kappa)*(z1 + kappa - eta)/(z0 + eta - kappa)/(z1 - eta - kappa))^(eta/kappa);
% Compute the main curve value as the difference
main_curve = left_hand_side - right_hand_side;
end
