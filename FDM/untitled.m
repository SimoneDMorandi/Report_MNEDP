clear
clc
format long
%%
% Solving -epsilon*x'' + beta*x' = 0
% u(x=0) = 0
% u(x=L) = 1
% f(x) = sin(pi*x)
epsilon = [0.1,1,10];
beta = 10;
L = 1;
N = 10;
h = L/(N+1);
x = linspace(0, L, N+2);
f = @(x) sin(pi*x);

% Cycle on epsilon values.
for j = 1:length(epsilon)
    epsi = epsilon(j);
    Pe = (beta*h)/(2*epsi)
    % Defining the exact solution.
    c1 = -(1)/(epsi*(exp(beta/epsi)-1));
    c2 = (1)/(beta*(exp(beta/epsi)-1));
    r = beta/epsi;
    Pi = pi*epsi*(pi^2 +(r^2)/pi);
    const1 = r/(Pi) + ((2*r)/(Pi))/(exp(r)-1);
    const2 = (-(2*r)/(Pi))/(exp(r)-1);
    const3 = 1/(epsi*(pi^2+((r^2)/(pi^2))));
    const4 = -r/(Pi);
    u = @(x) const1 + const2*exp(r*x) +const3*sin(pi*x) +const4*(cos(pi*x));
    % u = @(x) (1/(-1+exp(beta/epsi)))*(-1+exp((beta*x)/epsi)); %homog.
    
    % Centered FDM solution.
    % main_diag = 2*ones(N, 1);
    % lower_diag = -(Pe+1)*ones(N-1, 1);
    % upper_diag = (Pe-1)*ones(N-1, 1);
    % A = diag(main_diag,0)+diag(lower_diag,-1)+diag(upper_diag,1);
    % b = (h^2)/(epsi)*f(x(2:N+1))';
    % u_h = [0;A\b;0];
    % figure
    % hold on
    % plot(x, u(x), 'g', 'DisplayName', 'Exact Solution');
    % plot(x,u_h,'o-','DisplayName',['Finite Difference,\epsilon=',num2str(epsi)]);
    % xlabel('x');
    % ylabel('u(x)');
    % title('Centered Finite Difference');
    % hold on
    
    % Upwind solution.
    main_diag = (2+r*h)*ones(N,1);
    lower_diag = (-1-(h*r))*ones(N-1,1);
    upper_diag = -1*ones(N-1,1);
    A = diag(main_diag, 0)+diag(lower_diag,-1)+diag(upper_diag, 1);
    b = ((h^2)/epsi)*f(x(2:N+1))';
    u_h = [0;A\b;0];
    figure
    hold on
    plot(x, u(x), 'g', 'DisplayName', 'Exact Solution');
    plot(x, u_h, 'o-', 'DisplayName', ['Finite Difference, \epsilon = ', num2str(eps)]);
    xlabel('x (Domain L)');
    ylabel('u(x)');
    title('Upwind');
    grid on;
    
    
end