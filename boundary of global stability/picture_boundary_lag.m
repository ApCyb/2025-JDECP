close all;
clc;

k = 2/pi;

% Engineering parameters
tau_1 = 0.0448;
tau_2 = 0;

mu = pi*k - 1;

% Critical points (formula switches)
K_vco_ht = 1/(k*2*tau_1 + tau_2 + 2*sqrt(tau_1*(tau_1 + tau_2)))
 

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
    mu = pi*k - 1;
xi = (1 + k*tau_2*K_vco)/(2*sqrt((tau_1+tau_2)*K_vco));
eta = (k*tau_2*K_vco - mu)/(2*sqrt((tau_1+tau_2)*K_vco));
rho = sqrt(abs(xi^2 - k));
kappa = sqrt(eta^2 + k*mu);
    s_1 = ((kappa - eta)^2 + 2*xi*(kappa - eta) + k)/((kappa + eta)^2 - 2*xi*(kappa + eta) + k)*...
            exp(2*xi/rho*(atan(((xi - eta)^2 + rho^2 - kappa^2)/(2*rho*kappa)) + pi/2));
    omega_sep(i) = (sqrt(s_1) - 1)/(sqrt(s_1) + 1)*K_vco;
end
semilogx(K_vcos_separatrix_cycle, omega_sep./K_vcos_separatrix_cycle, 'blue', 'LineWidth', 2);
 
  

%Plot vertical lines for K_vco_minus, K_vco_plus and K_vco_eta_0
y = 0:0.0001:1;
semilogx(K_vco_ht*y./y, y, 'black--');
 


set(gca,'FontSize', 15)

xlabel('\textbf{$K_{\rm vco}$}','Interpreter','latex', 'fontsize', 20);  
ylabel('\textbf{$\frac{\omega_e^{\rm free}}{K_{\rm vco}}$}','Interpreter','latex', 'fontsize', 20); 
% legend('\omega_h', '\omega_{\rm sep}','\omega_{\rm ss}', 'Location','best');
% xticks([0.01 0.1 1 10 100 1000 10000 100000]);
% xticklabels({'10^{-2}', '10^{-1}' '10^0','10^1','10^2','10^3','10^4','10^5'});
axis([1 max(K_vcos) 0 1.2]);



