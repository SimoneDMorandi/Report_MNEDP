clear all
clc


%% Dirichlet non omogeneo oppure Neumann - Soluzione standard
% nu = 1.0;
% f = @(x,y) 32*(x*(1-x) + y*(1-y));
% gD = @(x,y) x+y;
% gN = @(x,y) [1 + 16.*(-1+2.*x).*(-1+y).*y; 1 + 16.*(-1+x).*x.*(-1+2.*y)];
% u_esatta = @(x,y) (x+y) + 16 * x .* (1 - x) .* y .* (1 - y);

%% Dirichlet omogeneo - Soluzione standard
% nu = 1.0;
% f = @(x,y) 32*(x*(1-x) + y*(1-y));
% gD = @(x,y) 0;
% gN = @(x,y) [16.*(-1+2.*x).*(-1+y).*y; 16.*(-1+x).*x.*(-1+2.*y)];
% u_esatta = @(x,y) 16 * x .* (1 - x) .* y .* (1 - y);

%% Dirichlet omogeneo - Paraboloide solito, va bene con tutto
% nu = 1.0;
% f = @(x,y) sin(pi*x).*sin(pi*y);
% gD = @(x,y) 0;
% gN = @(x,y) [(cos(pi*x) .* sin(pi*y))./(2*pi); (sin(pi*x) .* cos(pi*y))./(2*pi)];
% u_esatta = @(x,y) (sin(pi*x) .* sin(pi*y))./(2*pi^2);

%% Dirichlet omogeneo - Problema con dirichlet omogeneo e Neumann, da usare con marker [2 3 6 7]
% quindi qui ho diversi dati al bordo: N e D omogeneo
nu = 1.0;
f = @(x,y) -sin(pi*x).*exp(-2*y);
gD = @(x,y) 0;
gN = @(x,y) [pi*(cos(pi*x) .* exp(-2*y))./(4-pi^2); -2*(sin(pi*x) .* exp(-2*y))./(4-pi^2)];
u_esatta = @(x,y) (sin(pi*x) .* exp(-2*y))./(4-pi^2);

%%
% nu = 1.0;
% f = @(x,y)  4*pi^2 *(x.^2+y.^2).*cos(pi*(x.^2+y.^2)) + 4*pi*sin(pi*(x.^2+y.^2)) ;
% gD = @(x,y) 0;
% gN = @(x,y) [2*pi*x.*cos(pi*(x.^2+y.^2)); 2*pi*y.*cos(pi*(x.^2+y.^2))
%     ];
% u_esatta = @(x,y) log(1+(x.^2+y.^2)).*cos(5*pi*x.*y);

N = 3;
marker_triang = [2 3 6 7];
area_ref = 0.03;

area_min = 1e-3;
area_dec = logspace(log10(area_ref), log10(area_min), N);

errore_0 = zeros(N,1);
errore_1 = zeros(N,1);
h0 = zeros(N,1);
h1 = zeros(N,1);
count0 = 2;
count1 = 2;

for i=1:N
    area_ref = area_dec(i);
    geom = triangolatore(area_ref, marker_triang);
    [u_star, u0, uD] = poisson_DnO(geom, nu, f, gD, gN);
    
    [err, err_H1] = errore_DnO(u_esatta, u_star, geom, gN, u0, uD);
    if i == 1
        errore_0(i) = err;
        h0(i) = sqrt(area_ref);
    elseif i > 1 && abs(err-errore_0(count0-1))>1e-12
        errore_0(count0) = err;
        h0(count0) = sqrt(area_ref);
        count0 = count0 + 1;
    end


    %%err_H1 = (sum(errore_H1));
    if i == 1
        errore_1(i) = err_H1;
        h1(i) = sqrt(area_ref);
    elseif i > 1 && abs(err_H1-errore_1(count1-1))>1e-12
        errore_1(count1) = err_H1;
        h1(count1) = sqrt(area_ref);
        count1 = count1 + 1;
    end

end %i

%% Errori

% Fit lineare in scala log-log
p = polyfit(log(h0(1:count0-1)), log(errore_0(1:count0-1)), 1); % Considera i primi 4 punti per il fit
disp(['Pendenza stimata: ', num2str(p(1))]); % Stampa la pendenza stimata

% Grafico log-log
figure;
subplot(1,2,1);
loglog(h0(1:count0-1), errore_0(1:count0-1), 'o', 'MarkerSize', 8, 'LineWidth', 1.5); % Pallini
hold on;

% Aggiungi la retta di fit
h_fit = linspace(min(h0(1:count0-1)), max(h0(1:count0-1)), N);
errore_fit = exp(polyval(p, log(h_fit))); % Regressione in scala log-log
loglog(h_fit, errore_fit, '--', 'LineWidth', 1.5);

% Configurazione del grafico
xlabel('h');
ylabel('errore_0');
title('Convergenza in norma L^2', 'FontSize', 16);
legend('Punti calcolati', ['Fit, pendenza = ', num2str(p(1))], 'Location', 'SouthEast');
grid on;
hold off;


p = polyfit(log(h1(1:count1-1)), log(errore_1(1:count1-1)), 1); % Considera i primi 4 punti per il fit
disp(['Pendenza stimata: ', num2str(p(1))]); % Stampa la pendenza stimata

% Grafico log-log
%figure;
subplot(1,2,2);
loglog(h1(1:count1-1), errore_1(1:count1-1), 'o', 'MarkerSize', 8, 'LineWidth', 1.5); % Pallini
hold on;

% Aggiungi la retta di fit
h_fit = linspace(min(h1(1:count1-1)), max(h1(1:count1-1)), N);
errore_fit = exp(polyval(p, log(h_fit))); % Regressione in scala log-log
loglog(h_fit, errore_fit, '--', 'LineWidth', 1.5);

% Configurazione del grafico
xlabel('h');
ylabel('errore_1');
title('Convergenza in seminorma H^1', 'FontSize', 16);
legend('Punti calcolati', ['Fit, pendenza = ', num2str(p(1))], 'Location', 'SouthEast');
grid on;
hold off;


%% Plot della soluzione approssimata e della soluzione esatta
figure
subplot(1,2,1);
trisurf(geom.obj.T, geom.obj.P(:,1), geom.obj.P(:,2), u_star);
xlabel('X'); 
ylabel('Y'); 
zlabel('U'); 
title('Soluzione approssimata del problema di Poisson');
colorbar;

subplot(1,2,2);
trisurf(geom.obj.T, geom.obj.P(:,1), geom.obj.P(:,2), u_esatta(geom.obj.P(:,1), geom.obj.P(:,2)));
xlabel('X'); 
ylabel('Y'); 
zlabel('U'); 
title('Soluzione esatta del problema di Poisson');
colorbar;


