%% hh
clear
clc
format short e
% Heat equation with non-homog. Neumann conditions.

% Initial data.
nu = 1.0;
f = @(x,y,t) -sin(pi*x).*exp(-2*y) +t ;
gd = @(x,y,t) 0.5*(t^2);
gn = @(x,y,t) +[pi*(cos(pi*x) .* exp(-2*y))./(4-pi^2); -2*(sin(pi*x) .* exp(-2*y))./(4-pi^2)];
u = @(x,y,t) 0.5*(t^2)+ (sin(pi*x) .* exp(-2*y))./(4-pi^2);
marker_triang = [2 3 6 7];
area_ref = 0.03;
geom = triangulator_P2(area_ref,marker_triang);
N = geom.Nobj.N_node + geom.Nobj.N_edge - size(geom.piv.Di,1);
E = geom.obj.T;
P = geom.obj.P;
Piv = geom.piv.piv;
[u_star, u0, ud] = poisson_DnO_P2_t(geom, nu, f, gd, gn,1);
figure;
        % Subplot for ut
        subplot(1, 2, 1);
        trisurf(E(:,1:3), P(1:N, 1), P(1:N, 2), u(P(1:N,1),P(1:N,2),1));
        title('Exact Solution, t =0');
        xlabel('x');
        ylabel('y');
        zlabel('u(x, y, t(n))');
        colorbar;
        view(3);
        axis tight;
        % Subplot for u_h.
        subplot(1, 2, 2);
        trisurf(E(:,1:3), P(1:N, 1), P(1:N, 2), u_star);
        title('FEM Solution (u_h), t =0');
        xlabel('x');
        ylabel('y');
        zlabel('u_h(x, y, t(n))');
        colorbar;
        view(3);
        axis tight;
        drawnow;
%% Time evolution.

% Useful matrices.
[F,G,GD,B,Bd,K,Kd] = build_time(geom,nu,f,ud,gd,0);
un = u_h;
T = 10;
N_time = 1000;
delta_t = T/N_time;
for i = 1:N_time
    t(i) = i*delta_t;
end
A = B + (delta_t/2)*K;
for n = 1:N_time
    [Fn_1,udn_1] = build_time(geom,nu,f,ud,gd,t(n));
    ud_first = (1/delta_t)*(udn_1-udn);
    bt = B*un +(delta_t/2)*(-K*un + Fn -2*Bd*ud_first -Kd*udn ...
         +Fn_1 -Kd*udn_1);
    un = A\bt;
    % Plotting.
    u_star = zeros(geom.Nobj.N_node,1);
    for i = 1:geom.Nobj.N_node
        ii = Piv(i);
        if ii>0
            u_star(i) = un(ii);
        else
            u_star(i) = udn_1(-ii);
        end
    end
    if n == -10
        figure;
        % Subplot for ut
        subplot(1, 2, 1);
        trisurf(E, P(:, 1), P(:, 2), u(P(:,1),P(:,2),t(n)));
        title(['Exact Solution, t = ', num2str(t(n))]);
        xlabel('x');
        ylabel('y');
        zlabel('u(x, y, t(n))');
        colorbar;
        view(3);
        axis tight;
        % Subplot for u_h.
        subplot(1, 2, 2);
        trisurf(E,P(:, 1), P(:, 2), u_star);
        title(['FEM Solution (u_h), t = ', num2str(t(n))]);
        xlabel('x');
        ylabel('y');
        zlabel('u_h(x, y, t(n))');
        colorbar;
        view(3);
        axis tight;
        drawnow;
    end
    % Updating time variables.
    Fn = Fn_1;
    udn = udn_1;
end