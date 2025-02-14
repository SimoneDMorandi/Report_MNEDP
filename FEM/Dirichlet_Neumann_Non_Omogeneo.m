clear
clc
% Heat Equation with Non-Hom. Dirichlet conditions.

% Initial Data.
nu = 1.0;
f = @(x,y) -sin(pi*x).*exp(-2*y);
gd = @(x,y) 0;
gn = @(x,y) [pi*(cos(pi*x) .* exp(-2*y))./(4-pi^2); -2*(sin(pi*x) .* exp(-2*y))./(4-pi^2)];
u= @(x,y) (sin(pi*x) .* exp(-2*y))./(4-pi^2);
marker_triang = [2 3 6 7];
area_ref = 0.03;

geom = triangulator_P2(area_ref,marker_triang)
[A,b,ud] = assemble_poisson_D_N(geom,nu,f,gd,gn);
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
figure
subplot(1,2,1)
trisurf(geom.obj.T,geom.obj.P(:,1), geom.obj.P(:,2),u(geom.obj.P(:,1), geom.obj.P(:,2)))
title('Exact solution');
subplot(1,2,2)
trisurf(geom.obj.T,geom.obj.P(:,1), geom.obj.P(:,2),u_star)
title('u_h');

%% Calcolo dell'errore.
% n_iter = 10;
% h = linspace(0.020,0.030,n_iter);
% err_0 = zeros(n_iter,1);
% err_1 = zeros(n_iter,1);
% global ref_area;
% ref_area = 0.5;
% for i= 1:n_iter
%     ref_area = ref_area-(ref_area*i/10000);
%     geom = triangulator(h(i));
%     u_s = assemble_poisson_D_N(geom,nu,f,gn,gd);
%     [E0,E1] = calculate_err(geom,u,u_star,gn,du_dx,du_dy);
%     err_0(i) = E0;
%     err_1(i) = E1;
% end
% % Rappresentazione grafica dell'errore.
% figure;
% loglog(h, err_0, 'o', 'MarkerFaceColor', 'r');
% hold on;
% p = polyfit(log(h), log(err_0), 1);
% fitted_line = polyval(p, log(h));
% plot(h, exp(fitted_line), 'LineWidth', 2);
% figure
% loglog(h, err_1, 'o', 'MarkerFaceColor', 'r');
% hold on;
% p = polyfit(log(h), log(err_1), 1);
% fitted_line = polyval(p, log(h));
% plot(h, exp(fitted_line), 'LineWidth', 2);