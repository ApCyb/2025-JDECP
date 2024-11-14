close all;
clc;

k = 2/pi;

% Engineering parameters
tau_1 = 0.0448;
tau_2 = 0.0185;

mu = pi*k - 1;

% Critical points (formula switches)
K_vco_ht = (2*tau_1 + tau_2 - 2*sqrt(tau_1*(tau_1 + tau_2)))/(k*tau_2^2);
K_vco_fn = (2*tau_1 + tau_2 + 2*sqrt(tau_1*(tau_1 + tau_2)))/(k*tau_2^2)
K_vco_pt = mu/(k*tau_2);

% set of K_vco for modelling
K_vcos = [0.1:0.1:10 11:1:100 110:10:1000 1100:100:10000 11000:1000:100000];
n = length(K_vcos); 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hold-in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
semilogx(K_vcos, K_vcos./K_vcos, 'black', 'LineWidth', 1);
grid on;
hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pull-in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separatrix cycle exists for K_vco > K_vco_minus
K_vcos_separatrix_cycle = [K_vcos(K_vcos(:) > K_vco_ht)];

m = length(K_vcos_separatrix_cycle); 
omega_sep = zeros(1, m);
for i = 1:m
    K_vco = K_vcos_separatrix_cycle(i);
    omega_sep(i) = omega_sep_formula(tau_1, tau_2, k, K_vco);
end
semilogx(K_vcos_separatrix_cycle, omega_sep./K_vcos_separatrix_cycle, 'blue', 'LineWidth', 2);
omega_sep_formula(tau_1, tau_2, k, K_vco_fn)/K_vco_fn



K_vcos_semistable_cycle = [K_vcos(K_vcos(:) > K_vco_pt+20)];
omega_ss = zeros(1,length(K_vcos_semistable_cycle));
for j=1:length(K_vcos_semistable_cycle)
    K_vco = K_vcos_semistable_cycle(j);
    omega_ss(j) = omega_ss_formula(tau_1, tau_2, k, K_vco);
end
semilogx(K_vcos_semistable_cycle, omega_ss./K_vcos_semistable_cycle, 'red--', 'LineWidth', 1);
K_vco_fn;
omega_ss_formula(tau_1, tau_2, k, K_vco_fn)/K_vco_fn;


%%%asymptotic formula
a = tau_2/(tau_1 + tau_2)
    fcn = @(x) a*(2*x - a - x^2)/x/(x - a)- log((x^2*(1 - a)/(x-a)^2));
    start=a+eps;
    finish = sqrt(a);
    b = fzero(fcn, [start, finish])
    y_as = (-2*a*b + b^2 + a)/(2*b - b^2 - a)
    semilogx([1 max(K_vcos_semistable_cycle)], [y_as y_as], 'g', 'LineWidth', 1);

% Lyapunov fuction pull-in estimate
omega_ps_Lyapunov = zeros(1,n);
for j=1:n
    K_vco = K_vcos(j);
    omega_ps_Lyapunov(j) = K_vco*(tau_1/(2*sqrt(tau_2*(tau_1 + tau_2)) - 2*tau_2) - ...
        sqrt(tau_1^2/(2*sqrt(tau_2*(tau_1 + tau_2)) - 2*tau_2)^2-1));
end
semilogx(K_vcos, omega_ps_Lyapunov./K_vcos, 'black', 'LineWidth', 1);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lock-in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% omega_l = zeros(1, n);
% for i = 1:n
%     K_vco = K_vcos(i);
%     omega_l(i) = omega_l_formula(tau_1, tau_2, k, K_vco);
% end
% loglog(1./K_vcos, omega_l, 'green', 'LineWidth', 2);




%Plot vertical lines for K_vco_minus, K_vco_plus and K_vco_eta_0
y = 0:0.0001:1;
semilogx(K_vco_ht*y./y, y, 'black--');
semilogx(K_vco_pt*y./y, y, 'black--');



set(gca,'FontSize', 15)

xlabel('\textbf{$K_{\rm vco}$}','Interpreter','latex', 'fontsize', 20);  
ylabel('\textbf{$\frac{\omega_e^{\rm free}}{K_{\rm vco}}$}','Interpreter','latex', 'fontsize', 20); 
% legend('\omega_h', '\omega_{\rm sep}','\omega_{\rm ss}', 'Location','best');
% xticks([0.01 0.1 1 10 100 1000 10000 100000]);
% xticklabels({'10^{-2}', '10^{-1}' '10^0','10^1','10^2','10^3','10^4','10^5'});
axis([1 max(K_vcos) 0 1.2]);



