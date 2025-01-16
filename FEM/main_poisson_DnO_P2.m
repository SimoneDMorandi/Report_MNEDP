clear all
clc

%% Dati del problema

%% Dirichlet non omogeneo oppure Neumann
% nu = 1.0;
% f = @(x,y) 32*(x*(1-x) + y*(1-y));
% gD = @(x,y) x+y;
% gN = @(x,y) [1 + 16.*(-1+2.*x).*(-1+y).*y; 1 + 16.*(-1+x).*x.*(-1+2.*y)];
% u_esatta = @(x,y) (x+y) + 16 * x .* (1 - x) .* y .* (1 - y);

%% Dirichlet omogeneo
% nu = 1.0;
% f = @(x,y) 32*(x*(1-x) + y*(1-y));
% gD = @(x,y) 0;
% gN = @(x,y) [16.*(-1+2.*x).*(-1+y).*y; 16.*(-1+x).*x.*(-1+2.*y)];
% u_esatta = @(x,y) 16 * x .* (1 - x) .* y .* (1 - y);

%% Dirichlet omogeneo - Problema con dirichlet omogeneo e Neumann, da usare con marker [2 3 6 7]
% quindi qui ho diversi dati al bordo: N e D omogeneo
nu = 1.0;
f = @(x,y) -sin(pi*x).*exp(-2*y);
gD = @(x,y) 0;
gN = @(x,y) [pi*(cos(pi*x) .* exp(-2*y))./(4-pi^2); -2*(sin(pi*x) .* exp(-2*y))./(4-pi^2)];
u_esatta = @(x,y) (sin(pi*x) .* exp(-2*y))./(4-pi^2);


N_max = 3;
marker_triang = [2 3 6 7];

area_ref = 0.03;

area_min = 1e-3;
area_dec = logspace(log10(area_ref), log10(area_min), N_max);
errore_0 = zeros(N_max,1);
errore_1 = zeros(N_max,1);
h0 = zeros(N_max,1);
h1 = zeros(N_max,1);
count0 = 2;
count1 = 2;

for j=1:N_max

    area_ref = area_dec(j);
    geom = triangulator_P2(area_ref, marker_triang);
    N = geom.Nobj.N_node + geom.Nobj.N_edge - size(geom.piv.Di,1);
    
    [u_star, u0, uD] = poisson_DnO_P2(geom, nu, f, gD, gN);

    %[err, err_H1] = errore_per_P2(u_esatta, u_star, geom, gN, u0, uD);
    % 
    % if j == 1
    %     errore_0(j) = err;
    %     h0(j) = sqrt(area_ref);
    % elseif j > 1 && abs(err-errore_0(count0-1))>1e-12
    %     errore_0(count0) = err;
    %     h0(count0) = sqrt(area_ref);
    %     count0 = count0 + 1;
    % end
    % if j == 1
    %     errore_1(j) = err_H1;
    %     h1(j) = sqrt(area_ref);
    % elseif j > 1 && abs(err_H1-errore_1(count1-1))>1e-12
    %     errore_1(count1) = err_H1;
    %     h1(count1) = sqrt(area_ref);
    %     count1 = count1 + 1;
    % end
end

%% Errori

% % Fit lineare in scala log-log
% p = polyfit(log(h0(1:count0-1)), log(errore_0(1:count0-1)), 1);
% disp(['Pendenza stimata con L^2: ', num2str(p(1))]); % Stampa la pendenza stimata
% 
% % Grafico log-log
% figure;
% subplot(1,2,1);
% loglog(h0(1:count0-1), errore_0(1:count0-1), 'o', 'MarkerSize', 8, 'LineWidth', 1.5); % Pallini
% hold on;
% 
% % Aggiungi la retta di fit
% h_fit = linspace(min(h0(1:count0-1)), max(h0(1:count0-1)), N);
% errore_fit = exp(polyval(p, log(h_fit))); % Regressione in scala log-log
% loglog(h_fit, errore_fit, '--', 'LineWidth', 1.5);
% 
% % Configurazione del grafico
% xlabel('h');
% ylabel('errore_0');
% title('Convergenza in norma L^2', 'FontSize', 16);
% legend('Punti calcolati', ['Fit, pendenza = ', num2str(p(1))], 'Location', 'SouthEast');
% grid on;
% hold off;
% 
% 
% p = polyfit(log(h1(1:count1-1)), log(errore_1(1:count1-1)), 1); % Considera i primi 4 punti per il fit
% disp(['Pendenza stimata con H^1: ', num2str(p(1))]); % Stampa la pendenza stimata
% 
% % Grafico log-log
% %figure;
% subplot(1,2,2);
% loglog(h1(1:count1-1), errore_1(1:count1-1), 'o', 'MarkerSize', 8, 'LineWidth', 1.5); % Pallini
% hold on;
% 
% % Aggiungi la retta di fit
% h_fit = linspace(min(h1(1:count1-1)), max(h1(1:count1-1)), N);
% errore_fit = exp(polyval(p, log(h_fit))); % Regressione in scala log-log
% loglog(h_fit, errore_fit, '--', 'LineWidth', 1.5);
% 
% % Configurazione del grafico
% xlabel('h');
% ylabel('errore_1');
% title('Convergenza in seminorma H^1', 'FontSize', 16);
% legend('Punti calcolati', ['Fit, pendenza = ', num2str(p(1))], 'Location', 'SouthEast');
% grid on;
% hold off;


%% Plot della soluzione approssimata e della soluzione esatta
% figure
% subplot(1,2,1);
% trisurf(geom.obj.T(:,1:3), geom.obj.P(1:N,1), geom.obj.P(1:N,2), u_star);
% xlabel('X');
% ylabel('Y');
% zlabel('U');
% title('Soluzione approssimata del problema di Poisson');
% colorbar;
% 
% subplot(1,2,2);
% trisurf(geom.obj.T(:,1:3), geom.obj.P(1:N,1), geom.obj.P(1:N,2), u_esatta(geom.obj.P(1:N,1), geom.obj.P(1:N,2)));
% xlabel('X');
% ylabel('Y');
% zlabel('U');
% title('Soluzione esatta del problema di Poisson');
% colorbar;



