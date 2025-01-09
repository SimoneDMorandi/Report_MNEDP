clear
clc
format short e
% Heat equation with non-homog. Dirichlet conditions.

% Initial data.
run("sample_square_dirichlet_short.m");
E = geom.obj.T;
P = geom.obj.P;
Piv = geom.piv.piv;
u = @(x,y,t) 0.5*(t^2) + x+y +16*x.*(1-x).*y.*(1-y);
du_dx = @(x, y) 1 +16*(1-2*x).*y.*(1-y); 
du_dy = @(x, y) 1 +16*x.*(1-x).*(1-2*y);
f = @(x,y,t) 32*(x*(1-x) +y*(1-y)) + t;
nu = 1.0;
gd = @(x,y,t) x+y +0.5*(t^2);
gn = @(x,y) 0;
[A,b,ud] = assemble_poisson_D_N_t(geom,nu,f,gn,gd);
u_h = A\b;

% Time evolution.

% Useful matrices.
[Fn,udn,B,Bd,K,Kd] = build_time(geom,nu,f,ud,gd,0);
un = u_h;
T = 10;
N_time = 200;
delta_t = T/N_time;
for i = 1:N_time
    t(i) = i*delta_t;
end
for n = 1:N_time
    A = B + (delta_t/2)*K;
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
    if n == N_time
        figure;
        % Subplot for u.
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