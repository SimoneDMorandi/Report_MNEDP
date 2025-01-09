%%
clear
clc
format short e
% Equazione del calore con Dirichlet e Neumann non omogenee.

% Dati iniziali.
run("sample_square_dirichlet_short.m");
u = @(x,y) x+y +16*x.*(1-x).*y.*(1-y);
du_dx = @(x, y) 1 +16*(1-2*x).*y.*(1-y); 
du_dy = @(x, y) 1 +16*x.*(1-x).*(1-2*y);
f = @(x,y) 32*(x*(1-x) +y*(1-y));
nu = 1.0;
gd = @(x,y) x+y;
gn = @(x,y) x+y;
% function che costruisce il sistema lineare.
[A,b,ud] = assemble_poisson_D_N(geom,nu,f,gn,gd);
% Risoluzione del sistema lineare.
u_h = A\b;

% Creazione del vettore u_star per la rappresentazione grafica.
u_star= zeros(geom.Nobj.N_node,1);
for i = 1:geom.Nobj.N_node
    ii = geom.piv.piv(i);
    if ii >0
        u_star(i) = u_h(ii);
    else
        u_star(i) = ud(-ii);
    end
end

% Rappresentazione grafica della u_star.
% figure
% trisurf(geom.obj.T,geom.obj.P(:,1), geom.obj.P(:,2),u_star)
% figure
% trisurf(geom.obj.T,geom.obj.P(:,1), geom.obj.P(:,2),u(geom.obj.P(:,1), geom.obj.P(:,2)))

%% Calcolo dell'errore.
n_iter = 10;
h = linspace(0.020,0.030,n_iter);
err_0 = zeros(n_iter,1);
err_1 = zeros(n_iter,1);
global ref_area;
ref_area = 0.5;
for i= 1:n_iter
    ref_area = ref_area-(ref_area*i/10000);
    geom = triangulator(h(i));
    u_s = assemble_poisson_D_N(geom,nu,f,gn,gd);
    [E0,E1] = calculate_err(geom,u,u_star,gn,du_dx,du_dy);
    err_0(i) = E0;
    err_1(i) = E1;
end
% Rappresentazione grafica dell'errore.
figure;
loglog(h, err_0, 'o', 'MarkerFaceColor', 'r');
hold on;
p = polyfit(log(h), log(err_0), 1);
fitted_line = polyval(p, log(h));
plot(h, exp(fitted_line), 'LineWidth', 2);
figure
loglog(h, err_1, 'o', 'MarkerFaceColor', 'r');
hold on;
p = polyfit(log(h), log(err_1), 1);
fitted_line = polyval(p, log(h));
plot(h, exp(fitted_line), 'LineWidth', 2);