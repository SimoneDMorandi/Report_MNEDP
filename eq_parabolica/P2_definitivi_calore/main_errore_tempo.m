close all;
clear all;
clc

%% Dati del problema con errore temporale
nu = 1;
u_esatta = @(x,y,t) x.^2 + y.^2 + exp(t); % conv in spazio
f = @(x,y,t) exp(t) - 4;
u0 = @(x,y) x.^2 + y.^2 + 1;
gD = @(x,y,t) x.^2 + y.^2 + exp(t);
gN =  @(x,y,t) [2*x; 2*y];


% T va da 0 a 1
n_t = 10;
delta_t = 1/n_t;

N_max = 5;
area_ref = 0.02;

errore_0_ei = zeros(N_max,1);
errore_1_ei = zeros(N_max,1);
h0_ei = zeros(N_max,1);
h1_ei = zeros(N_max,1);
count0_ei = 2;
count1_ei = 2;

errore_0_cn = zeros(N_max,1);
errore_1_cn = zeros(N_max,1);
h0_cn = zeros(N_max,1);
h1_cn = zeros(N_max,1);
count0_cn = 2;
count1_cn = 2;


geom = triangolatore(area_ref, [2 4 6 8], [0 0 0 0]);
[N, A, B, AD, BD] = build_A_B(geom, nu);

for i=1:N_max
    n_t = 20 * i;
    delta_t = 1/n_t;
    [u_ei, u_cn, u, uD] = build_solution(geom, A, B, AD, BD, gN, gD, f, u0, n_t, delta_t, u_esatta);

    u_full_ei = zeros(size(geom.obj.P, 1), n_t+1);
    u_full_ei(geom.piv.piv > 0,:) = u_ei(:,:);
    u_full_ei(geom.piv.piv < 0,:) = uD(:,:);

    u_full_cn = zeros(size(geom.obj.P, 1), n_t+1);
    u_full_cn(geom.piv.piv > 0,:) = u_cn(:,:);
    u_full_cn(geom.piv.piv < 0,:) = uD(:,:);

    [err_L2, err_H1] = errore(u_esatta, u_full_ei, geom, gN, u_ei, uD(:,end), n_t+1, delta_t * n_t);
    if i == 1
        errore_0_ei(i) = err_L2;
        h0_ei(i) = delta_t;
    elseif i > 1 && abs(err_L2-errore_0_ei(count0_ei-1))>1e-12
        errore_0_ei(count0_ei) = err_L2;
        h0_ei(count0_ei) = delta_t;
        count0_ei = count0_ei + 1;
    end

    % Errore per CRANK NICOLSON
    [err_L2, err_H1] = errore(u_esatta, u_full_cn, geom, gN, u_cn, uD(:,end), n_t + 1, delta_t * n_t);
    if i == 1
        errore_0_cn(i) = err_L2;
        h0_cn(i) = delta_t;
    elseif i > 1 && abs(err_L2-errore_0_cn(count0_cn-1))>1e-12
        errore_0_cn(count0_cn) = err_L2;
        h0_cn(count0_cn) = delta_t;
        count0_cn = count0_cn + 1;
    end

end %i

%% Grafico convergenza errori

p = polyfit(log(h0_ei(1:count0_ei-1)), log(errore_0_ei(1:count0_ei-1)), 1); 
disp(['Pendenza stimata: ', num2str(p(1))]); 

figure;
subplot(1,2,1);
loglog(h0_ei(1:count0_ei-1), errore_0_ei(1:count0_ei-1), 'o', 'MarkerSize', 8, 'LineWidth', 1.5); 
hold on;

h_fit = linspace(min(h0_ei(1:count0_ei-1)), max(h0_ei(1:count0_ei-1)), N);
errore_fit = exp(polyval(p, log(h_fit))); 
loglog(h_fit, errore_fit, '--', 'LineWidth', 1.5);

xlabel('\Deltat');
ylabel('Errore in L^2');
title('Convergenza Eulero implicito', 'FontSize', 16);
legend('Punti calcolati', ['Pendenza = ', num2str(p(1))], 'Location', 'SouthEast');
grid on;
hold off;

p = polyfit(log(h0_cn(1:count0_cn-1)), log(errore_0_cn(1:count0_cn-1)), 1); 
disp(['Pendenza stimata: ', num2str(p(1))]); 

subplot(1,2,2);
loglog(h0_cn(1:count0_cn-1), errore_0_cn(1:count0_cn-1), 'o', 'MarkerSize', 8, 'LineWidth', 1.5); 
hold on;

h_fit = linspace(min(h0_cn(1:count0_cn-1)), max(h0_cn(1:count0_cn-1)), N);
errore_fit = exp(polyval(p, log(h_fit))); 
loglog(h_fit, errore_fit, '--', 'LineWidth', 1.5);

xlabel('\Deltat');
ylabel('Errore in L^2');
title('Convergenza Crank-Nicolson', 'FontSize', 16);
legend('Punti calcolati', ['Pendenza = ', num2str(p(1))], 'Location', 'SouthEast');
grid on;
hold off;


%% Grafico istante finale
figure
u_esatta_full = u_esatta(geom.obj.P(:,1), geom.obj.P(:,2), delta_t*(n_t));

subplot(1, 3, 1);
trisurf(geom.obj.T(:,1:3), geom.obj.P(:,1), geom.obj.P(:,2), u_full_ei(:,end));
shading interp;
title(['Soluzione Numerica con Eulero implicito per t = ', num2str(delta_t*n_t)]);
xlabel('x'); ylabel('y'); zlabel('u');
colorbar;
view(3);

subplot(1, 3, 2);
trisurf(geom.obj.T(:,1:3), geom.obj.P(:,1), geom.obj.P(:,2), u_full_cn(:,end));
shading interp;
title(['Soluzione Numerica con Crank-Nicolson per t = ', num2str(delta_t*n_t)]);
xlabel('x'); ylabel('y'); zlabel('u');
colorbar;
view(3);

subplot(1, 3, 3);
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
    pause(0.5);
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