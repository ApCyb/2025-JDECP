close all;
clc;
clear all;

k = 2/pi;
mu = pi*k - 1;

% a = (2000*sqrt(2) - 1)/200000/20
a = 0.005

x_array = [0.01:0.01:10 11:0.1:100 110:1:1000 1100:10:10000];
% x = (tau_1 + tau_2)K_vco

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hold-in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pull-in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xi = (a*k*x_array + 1)./(2*sqrt(x_array));
eta = (a*k*x_array - mu)./(2*sqrt(x_array));
rho = sqrt(abs(xi.^2 - k));
kappa = sqrt(eta.^2 + mu*k);



    x_ht = 1/(k*(2 - a + 2*sqrt(1 - a)));
    x_pt = mu/(a*k);
    x_nf = 1/(k*(2 - a - 2*sqrt(1 - a)));

    
    semilogx([min(x_array) x_ht], [1 1], 'black', 'LineWidth', 1);
    hold on;
grid on;
    x_sep = x_array(x_array(:) >= x_ht & x_array(:) <= x_pt);
    y_sep = omega_sep_for_diagram(a, x_sep, k);
    semilogx(x_sep, y_sep, 'black', 'LineWidth', 1);
    semilogx(x_pt, omega_sep_for_diagram(a, x_pt, k), 'x','LineWidth', 2, 'Color', 'black');

    x_ss = x_array(x_array(:) > max(x_pt, x_ht));
    y_ss = omega_ss_for_diagram(a, x_ss, k);
    semilogx(x_ss, y_ss, 'black', 'LineWidth', 1);


%%%%%%%%%%%%%%%%%%%%%%%%% Cahn %%%%%%%%%%%%%%%%%%%%%%%%%

    y_Cahn = sqrt(pi./(x_array) + 2*a);
    Cahn_est = semilogx(x_array, y_Cahn, 'green', 'LineWidth', 2, 'DisplayName','Cahn');

%%%%%%%%%%%%%%%%%%%%%%%%% Brunk Rosenkranz %%%%%%%%%%%%%%%%%%%%%%%%%
x_BR = pi/2*x_array;
D = 1/2.*(1./sqrt(x_BR) + a*sqrt(x_BR));
y_BR = tanh(1/4*log((D + sqrt(x_BR)*(1 - a) + sqrt(D.^2 + 1 - a))./...
    (D + sqrt(x_BR)*(1 - a) - sqrt(D.^2 + 1 - a))) +...
    D./(2*sqrt(1 - D.^2)).*...
    (atan((1./sqrt(x_BR) - sqrt(D.^2 + 1 - a))./sqrt(1 - D.^2)) -...
    atan((1./sqrt(x_BR) + sqrt(D.^2 + 1 - a))./sqrt(1 - D.^2)) + pi));

%   BR_est_as = semilogx([0.01 max(x_BR)], [sqrt(4/3*a-a^2/3) sqrt(4/3*a-a^2/3)], 'red', 'LineWidth', 1, 'DisplayName','BrunkRosenkranz');
% BR_est = semilogx(x_BR, y_BR, 'blue', 'LineWidth', 1, 'DisplayName','BrunkRosenkranz');

%%%%%%%%%%%%%%%%%%%%%%%%% Lindsey %%%%%%%%%%%%%%%%%%%%%%%%%

  Lindsey_est = semilogx([0.01 max(x_array)], [sqrt(4/3*a) sqrt(4/3*a)], 'b', 'LineWidth', 2, 'DisplayName','Lindsey');
 %%%%%%%%%%%%%%%%%%%%%%%%% Mengali %%%%%%%%%%%%%%%%%%%%%%%%% 
 
    x_Mengali = x_array;
    y_Mengali = 2*sqrt(a/3 + pi./(4*x_Mengali) + (pi./(x_Mengali)).^2) - 2*pi./x_Mengali;
    Mengali_est = semilogx(x_Mengali, y_Mengali, 'm', 'LineWidth', 2, 'DisplayName','Mengali');

%%%%%%%%%%%%%%%%%%%%%%%%% Endo %%%%%%%%%%%%%%%%%%%%%%%%% 
    
    y_Endo = 1/4*(pi + 4 + 2*a*x_array).*sqrt(pi./(8*x_array));
    Endo_est = semilogx(x_array, y_Endo, 'r', 'LineWidth', 3, 'DisplayName','Endo');

     %%%%%%%%%%%%%%%%%%%%%%%%% Best %%%%%%%%%%%%%%%%%%%%%%%%% 
 
    semilogx([0.01 x_pt], [sqrt(a) sqrt(a)], 'color', [1.0 0.5 0.0], 'LineWidth', 3, 'DisplayName','Best');

    x_Best = x_array(x_array(:) > x_pt);
    y_Best = sqrt(a + pi./(2*x_Best));
    Best_est = semilogx(x_Best, y_Best, 'color', [1.0 0.5 0.0], 'LineWidth', 3, 'DisplayName','Best2');


set(gca,'FontSize', 15);

xlabel('\textbf{$(\tau_1 + \tau_2)K_{\rm vco}$}','Interpreter','latex', 'fontsize', 20);  
ylabel('\textbf{$\frac{\omega_p}{K_{\rm vco}}$}','Interpreter','latex', 'fontsize', 20); 
legend([Cahn_est Lindsey_est Mengali_est Best_est Endo_est], 'Location','best');

axis([0.09 max(x_array) 0 1.01]);
xticks([0.01 0.1 1 10 100 1000 5000 20*pi*100000]);
xticklabels({'10^{-2}','10^{-1}','10^0','10^1','10^2','10^3','5\cdot10^3','20\pi\cdot10^5'});



