close all;
clc;
clear all;

k = 2/pi;
mu = pi*k - 1;


%btw, standard engineering parameters tau_1 = 0.0448, tau_2 = 0.0185 correspond to a = 0.2923

x = [0.0001:0.01:10 11:0.1:100 110:1:1000 1100:10:10000];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hold-in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
semilogx([min(x) max(x)], [1 1], 'green', 'LineWidth', 1);
hold on;
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pull-in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 a = 0.2923;

xi = (a*k*x + 1)./(2*sqrt(x));
eta = (a*k*x - mu)./(2*sqrt(x));
rho = sqrt(abs(xi.^2 - k));
kappa = sqrt(eta.^2 + mu*k);


    x_minus = 1/(k*(2 - a + 2*sqrt(1 - a)));
    x_plus = 1/(k*(2 - a - 2*sqrt(1 - a)));
    x_eta_0 = mu/(a*k);
    %%%%%% K_vco < K_vco- %%%%%%
    semilogx([min(x) x_minus], [1 1], 'black', 'LineWidth', 1);

    %%%%%% K_vco- < K_vco < K_vco_eta=0 %%%%%%
    x_sep = x(x(:) > x_minus & x(:) < x_eta_0);
    y_sep = omega_sep_for_diagram(a, x_sep, k);
    semilogx(x_sep, y_sep, 'black', 'LineWidth', 1);

    %%%% K_vco- < K_vco < K_vco_eta=0 %%%%%%
    x_ss = x(x(:) > x_eta_0);
    y_ss = omega_ss_for_diagram(a, x_ss, k);
    semilogx(x_ss, y_ss, 'black', 'LineWidth', 1);


     %%%Safonov's asymptotic formula
    fcn = @(theta) (theta^2/(sinh(theta))^2 - (1-a));
    start=eps;
    finish = 1/eps;
    theta_0 = fzero(fcn, [start, finish]);
    y_as = (sinh(theta_0)*cosh(theta_0)-theta_0)/(sinh(theta_0))^2
    semilogx([1 max(x)], [y_as y_as], 'red--', 'LineWidth', 1);

%     semilogx([x_eta_0 x_eta_0], [-10 10], 'black', 'LineWidth', 1);
    
    %lock-in
    y = omega_l_for_diagram(a, x, k);
    semilogx(x, y, 'blue', 'LineWidth', 1);

    %Lyapunov
    lyapunov_est = (1 + sqrt(a) - sqrt((1 + sqrt(a))^2 - 4*a))/(2*sqrt(a));
    semilogx([min(x) max(x)], [lyapunov_est lyapunov_est], 'red', 'LineWidth', 1);



    %lower boundary for pull-in
    r_1 = ((k*sqrt(a*x) + eta).^2 - kappa.^2) ./ ((k*sqrt(a*x) - eta).^2 - kappa.^2) .*...
    (((k*sqrt(a*x) + kappa).^2 - eta.^2) ./ ((k*sqrt(a*x) - kappa).^2 - eta.^2)).^(eta./kappa);

    r_21 = (k*a*x + 2*xi.*sqrt(a*x) + 1) ./ (k*a*x - 2*xi.*sqrt(a*x) + 1) .*...
    exp(2*xi./rho.*(atan((1 - k*a*x) ./ (2*rho.*sqrt(a*x))) + pi/2));
    r_22 = ((k*sqrt(a*x) + xi).^2 - rho.^2) ./ ((k*sqrt(a*x) - xi).^2 - rho.^2) .*...
    (((k*sqrt(a*x) + rho).^2 - xi.^2) ./ ((k*sqrt(a*x) - rho).^2 - xi.^2)).^(xi./rho);
    
    r_2 = [r_21(x(:) < x_plus), r_22(x(:) > x_plus)];
    r = max(r_1, r_2);
    y = (sqrt(r) - 1)./(sqrt(r) + 1);
%     semilogx(x, y, 'm', 'LineWidth', 1);









set(gca,'FontSize', 15)

xlabel('\textbf{$(\tau_1 + \tau_2)K_{\rm vco}$}','Interpreter','latex', 'fontsize', 20);  
ylabel('\textbf{$\frac{\omega_e^{\rm free}}{K_{\rm vco}}$}','Interpreter','latex', 'fontsize', 20); 
% legend('\omega_h', '\omega_{\rm sep}','\omega_{\rm semi}', '\omega_l', 'Location','best');
axis([0.1 max(x) 0 1.2]);
xticks([0.01 0.1 1 10 100 1000 5000]);
xticklabels({'10^{-2}','10^{-1}','10^0','10^1','10^2','10^3','5\cdot10^3'});













