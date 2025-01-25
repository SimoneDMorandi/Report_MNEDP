clear all
close all
clc

%% Dirichlet omogeneo - Problema con dirichlet omogeneo e Neumann, da usare con marker [2 3 6 7]
% quindi qui ho diversi dati al bordo: N e D omogeneo
% nu = 1.0;
% f = @(x,y) -sin(pi*x).*exp(-2*y);
% gD = @(x,y) 0;
% gN = @(x,y) [pi*(cos(pi*x) .* exp(-2*y))./(4-pi^2); -2*(sin(pi*x) .* exp(-2*y))./(4-pi^2)];
% u_esatta = @(x,y) (sin(pi*x) .* exp(-2*y))./(4-pi^2);


%% Dati definitivi del problema considerato
nu = 1.0;

% f = @(x,y) sin(pi*x).*sin(pi*y);
% gD = @(x,y) 0;
% gN = @(x,y) [(cos(pi*x) .* sin(pi*y))./(2*pi); (sin(pi*x) .* cos(pi*y))./(2*pi)];
% u_esatta = @(x,y)  (sin(pi*x) .* sin(pi*y))./(2*pi^2);

f = @(x,y) - 4 + sin(pi*x).*sin(pi*y);
gD = @(x,y) x.^2 + y.^2;
gN = @(x,y) [2*x + (cos(pi*x) .* sin(pi*y))./(2*pi); 2*y + (sin(pi*x) .* cos(pi*y))./(2*pi)]; 
u_esatta = @(x,y) x.^2 + y.^2 + (sin(pi*x) .* sin(pi*y))./(2*pi^2);


N_max = 6;

area_ref = 0.02;

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
    geom = triangolatore_per_P2(area_ref, [1 3 5 7], [1 1 1 1]);
    N = geom.Nobj.N_node + geom.Nobj.N_edge - size(geom.piv.Di,1);
    
    [u_star, u0, uD] = poisson_DnO(geom, nu, f, gD, gN);

    [err, err_H1] = errore_per_P2(u_esatta, u_star, geom, gN, u0, uD);

    if j == 1
        errore_0(j) = err;
        h0(j) = sqrt(area_ref);
    elseif j > 1 && abs(err-errore_0(count0-1))>1e-12
        errore_0(count0) = err;
        h0(count0) = sqrt(area_ref);
        count0 = count0 + 1;
    end
    if j == 1
        errore_1(j) = err_H1;
        h1(j) = sqrt(area_ref);
    elseif j > 1 && abs(err_H1-errore_1(count1-1))>1e-12
        errore_1(count1) = err_H1;
        h1(count1) = sqrt(area_ref);
        count1 = count1 + 1;
    end
end

%% Errori

p = polyfit(log(h0(1:count0-1)), log(errore_0(1:count0-1)), 1); 
disp(['Pendenza stimata: ', num2str(p(1))]); 

figure;
subplot(1,2,1);
loglog(h0(1:count0-1), errore_0(1:count0-1), 'o', 'MarkerSize', 8, 'LineWidth', 1.5); 
hold on;

h_fit = linspace(min(h0(1:count0-1)), max(h0(1:count0-1)), N);
errore_fit = exp(polyval(p, log(h_fit))); 
loglog(h_fit, errore_fit, '--', 'LineWidth', 1.5);

xlabel('Passo di discretizzazione h');
ylabel('Errore in L^2');
title('Convergenza in norma L^2', 'FontSize', 16);
legend('Punti calcolati', ['Pendenza = ', num2str(p(1))], 'Location', 'SouthEast');
grid on;
hold off;


p = polyfit(log(h1(1:count1-1)), log(errore_1(1:count1-1)), 1); 
disp(['Pendenza stimata: ', num2str(p(1))]); 


subplot(1,2,2);
loglog(h1(1:count1-1), errore_1(1:count1-1), 'o', 'MarkerSize', 8, 'LineWidth', 1.5); 
hold on;

h_fit = linspace(min(h1(1:count1-1)), max(h1(1:count1-1)), N);
errore_fit = exp(polyval(p, log(h_fit))); 
loglog(h_fit, errore_fit, '--', 'LineWidth', 1.5);

xlabel('Passo di discretizzazione h');
ylabel('Errore in H^1');
title('Convergenza in seminorma H^1', 'FontSize', 16);
legend('Punti calcolati', ['Pendenza = ', num2str(p(1))], 'Location', 'SouthEast');
grid on;
hold off;


%% Plot della soluzione approssimata e della soluzione esatta
figure
subplot(1,2,1);
trisurf(geom.obj.T(:,1:3), geom.obj.P(:,1), geom.obj.P(:,2), u_star);
xlabel('X');
ylabel('Y');
zlabel('U');
title('Soluzione approssimata del problema');
colorbar;

subplot(1,2,2);
trisurf(geom.obj.T(:,1:3), geom.obj.P(:,1), geom.obj.P(:,2), u_esatta(geom.obj.P(:,1), geom.obj.P(:,2)));
xlabel('X');
ylabel('Y');
zlabel('U');
title('Soluzione esatta del problema');
colorbar;




