close all;
clear all;
clc

%% Dati del problema con errore spaziale
nu = 1;
u_esatta = @(x,y,t) exp(x) + exp(y) + t; 
f = @(x,y,t) 1 - exp(x) - exp(y);
u0 = @(x,y) exp(x) + exp(y);
gD = @(x,y,t) exp(x) + exp(y) + t;
gN =  @(x,y,t) [exp(x); exp(y)];

% T va da 0 a 1
n_t = 50;
delta_t = 1/n_t;

N_max = 3;
area_ref = 0.03;

area_min = 0.01;
area_dec = logspace(log10(area_ref), log10(area_min), N_max);

errore_0_cn = zeros(N_max,1);
errore_1_cn = zeros(N_max,1);
h0_cn = zeros(N_max,1);
h1_cn = zeros(N_max,1);
count0_cn = 2;
count1_cn = 2;

for i=1:N_max
    area_ref = area_dec(i);
    geom = triangolatore(area_ref, [2 4 6 8], [0 0 0 0]);
    
    [N, A, B, AD, BD] = build_A_B(geom, nu);
    [u_cn, u, uD] = build_solution_cn(geom, A, B, AD, BD, gN, gD, f, u0, n_t, delta_t, u_esatta);

    u_full_cn = zeros(size(geom.obj.P, 1), n_t+1);
    u_full_cn(geom.piv.piv > 0,:) = u_cn(:,:);

    % Errore per CRANCK NICOLSON
    [err_L2, err_H1] = errore(u_esatta, u_full_cn, geom, gN, u_cn, uD(:,end), n_t + 1, delta_t * n_t);

    if i == 1
        errore_0_cn(i) = err_L2;
        h0_cn(i) = sqrt(area_ref);
    elseif i > 1 && abs(err_L2-errore_0_cn(count0_cn-1))>1e-12
        errore_0_cn(count0_cn) = err_L2;
        h0_cn(count0_cn) = sqrt(area_ref);
        count0_cn = count0_cn + 1;
    end

    if i == 1
        errore_1_cn(i) = err_H1;
        h1_cn(i) = sqrt(area_ref);
    elseif i > 1 && abs(err_H1-errore_1_cn(count1_cn-1))>1e-12
        errore_1_cn(count1_cn) = err_H1;
        h1_cn(count1_cn) = sqrt(area_ref);
        count1_cn = count1_cn + 1;
    end

end %i

%% Grafico errore CN

p = polyfit(log(h0_cn(1:count0_cn-1)), log(errore_0_cn(1:count0_cn-1)), 1); 
disp(['Pendenza stimata: ', num2str(p(1))]); 

figure;
subplot(1,2,1);
loglog(h0_cn(1:count0_cn-1), errore_0_cn(1:count0_cn-1), 'o', 'MarkerSize', 8, 'LineWidth', 1.5); 
hold on;

h_fit = linspace(min(h0_cn(1:count0_cn-1)), max(h0_cn(1:count0_cn-1)), N);
errore_fit = exp(polyval(p, log(h_fit))); 
loglog(h_fit, errore_fit, '--', 'LineWidth', 1.5);

xlabel('Passo di discretizzazione h');
ylabel('Errore in L^2');
title('Convergenza in norma L^2 - Crank-Nicolson', 'FontSize', 16);
legend('Punti calcolati', ['Pendenza = ', num2str(p(1))], 'Location', 'SouthEast');
grid on;
hold off;

p = polyfit(log(h1_cn(1:count1_cn-1)), log(errore_1_cn(1:count1_cn-1)), 1); 
disp(['Pendenza stimata: ', num2str(p(1))]); 


subplot(1,2,2);
loglog(h1_cn(1:count1_cn-1), errore_1_cn(1:count1_cn-1), 'o', 'MarkerSize', 8, 'LineWidth', 1.5); 
hold on;

h_fit = linspace(min(h1_cn(1:count1_cn-1)), max(h1_cn(1:count1_cn-1)), N);
errore_fit = exp(polyval(p, log(h_fit))); 
loglog(h_fit, errore_fit, '--', 'LineWidth', 1.5);

xlabel('Passo di discretizzazione h');
ylabel('Errore in H^1');
title('Convergenza in seminorma H^1 - Crank-Nicolson', 'FontSize', 16);
legend('Punti calcolati', ['Pendenza = ', num2str(p(1))], 'Location', 'SouthEast');
grid on;
hold off;

%% Grafico istante finale
figure
u_esatta_full = u_esatta(geom.obj.P(:,1), geom.obj.P(:,2), delta_t*(n_t));

subplot(1, 2, 1);
trisurf(geom.obj.T(:,1:3), geom.obj.P(:,1), geom.obj.P(:,2), u_full_cn(:,end));
shading interp;
title(['Soluzione Numerica con Crank-Nicolson per t = ', num2str(delta_t*n_t)]);
xlabel('x'); ylabel('y'); zlabel('u');
colorbar;
view(3);

subplot(1, 2, 2);
trisurf(geom.obj.T(:,1:3), geom.obj.P(:,1), geom.obj.P(:,2), u_esatta_full);
shading interp;
title(['Soluzione Esatta per t = ', num2str(delta_t*n_t)]);
xlabel('x'); ylabel('y'); zlabel('u_{esatta}');
colorbar;
view(3);


%% Grafici per ricerca errori
%{
figure;
for n = 1:n_t+1

    u_esatta_full = u_esatta(geom.obj.P(:,1), geom.obj.P(:,2), delta_t*(n-1));

    % Subplot per la soluzione numerica
    subplot(1, 2, 1);
    trisurf(geom.obj.T(:,1:3), geom.obj.P(:,1), geom.obj.P(:,2), u_full_cn(:,n));
    shading interp;                         
    title(['Soluzione Numerica a t = ', num2str(delta_t*n)]);
    xlabel('x'); ylabel('y'); zlabel('u');
    colorbar;
    view(3);

    % Subplot per la soluzione esatta
    subplot(1, 2, 2);
    trisurf(geom.obj.T(:,1:3), geom.obj.P(:,1), geom.obj.P(:,2), u_esatta_full);
    shading interp;                         
    title(['Soluzione Esatta a t = ', num2str(delta_t*n)]);
    xlabel('x'); ylabel('y'); zlabel('u_{esatta}');
    colorbar;
    view(3);

    % Pausa per animazione
    pause(1);
end

figure;
for n = 1:n_t+1

    u_esatta_full = u_esatta(geom.obj.P(:,1), geom.obj.P(:,2), delta_t*(n-1));

    % Subplot per la soluzione numerica
    subplot(1, 2, 1);
    trisurf(geom.obj.T(:,1:3), geom.obj.P(:,1), geom.obj.P(:,2), u_esatta_full - u_full_cn(:,n));
    shading interp;                         
    title(['Soluzione errore cn a t = ', num2str(delta_t*n)]);
    xlabel('x'); ylabel('y'); zlabel('u');
    colorbar;
    view(3);

    % Subplot per la soluzione esatta
    subplot(1, 2, 2);
    trisurf(geom.obj.T(:,1:3), geom.obj.P(:,1), geom.obj.P(:,2), u_esatta_full - u_full_ei(:,n));
    shading interp;                         
    title(['Soluzione errore ei a t = ', num2str(delta_t*n)]);
    xlabel('x'); ylabel('y'); zlabel('u_{esatta}');
    colorbar;
    view(3);

    % Pausa per animazione
    pause(1);
end
%}