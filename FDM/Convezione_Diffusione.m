% Convezione e Diffusione in una dimensione con FDM e FEM.
% Risolvo -epsilon*u(x)'' + (d/dx)beta*(x) = 0.

clear
clc
format long
% plotting == 1 per mostrare i grafici.
plotting = 0;
% Dati iniziali e inizializzazione errori.
epsilon = 0.10;
beta = 10;
L = 1;
N = 10;
f = @(x) sin(pi*x);
err_FDM = zeros(length(N),1);
err_upwind = zeros(length(N),1);
err_SUPG = zeros(length(N),1);

% Ciclo sui valori di epsilon (Pe).
for j = 1:length(epsilon)
    % Ciclo sui valori di N (errori).
    for k = 1:length(N)
    h = L/(N(k)+1);
    x = linspace(0, L, N(k)+2);
    epsi = epsilon(j);
    Pe = (beta*h)/(2*epsi)

    % Soluzione esatta.
    c1 = -(1)/(epsi*(exp(beta/epsi)-1));
    c2 = (1)/(beta*(exp(beta/epsi)-1));
    r = beta/epsi;
    Pi = pi*epsi*(pi^2 +(r^2)/pi);
    const1 = r/(Pi) + ((2*r)/(Pi))/(exp(r)-1);
    const2 = (-(2*r)/(Pi))/(exp(r)-1);
    const3 = 1/(epsi*(pi^2+((r^2)/(pi^2))));
    const4 = -r/(Pi);
    u = @(x) const1 + const2*exp(r*x) +const3*sin(pi*x) +const4*(cos(pi*x));
    % u = @(x) (1/(-1+exp(beta/epsi)))*(-1+exp((beta*x)/epsi)); % Soluzione omogenea.
    
    % Soluzione con FDM.
    main_diag = 2*ones(N(k), 1);
    lower_diag = -(Pe+1)*ones(N(k)-1, 1);
    upper_diag = (Pe-1)*ones(N(k)-1, 1);
    A = diag(main_diag,0)+diag(lower_diag,-1)+diag(upper_diag,1);
    b = (h^2)/(epsi)*f(x(2:N(k)+1))';
    u_h = [0;A\b;0];
    err_FDM(k) = sqrt(h * sum((u(x) - u_h').^2)) 

    % Plot
    if plotting == 1
    figure
    hold on
    plot(x, u(x), 'g', 'DisplayName', 'Exact Solution');
    plot(x,u_h,'o-','DisplayName',['Finite Difference,Pe=',num2str(Pe)]);
    xlabel('x');
    ylabel('u(x)');
    title('Centered Finite Difference');
    grid on
    legend show
    hold off
    end

    % Soluzione Upwind.
    main_diag = (2+r*h)*ones(N(k),1);
    lower_diag = (-1-(h*r))*ones(N(k)-1,1);
    upper_diag = -1*ones(N(k)-1,1);
    A = diag(main_diag, 0)+diag(lower_diag,-1)+diag(upper_diag, 1);
    b = ((h^2)/epsi)*f(x(2:N(k)+1))';
    u_h = [0;A\b;0];
    err_upwind(k) = sqrt(h * sum((u(x) - u_h').^2))
    if plotting == 1
    figure
    hold on
    plot(x, u(x), 'g', 'DisplayName', 'Exact Solution');
    plot(x, u_h, 'o-', 'DisplayName', ['Finite Difference, Pe = ', num2str(Pe)]);
    xlabel('x');
    ylabel('u(x)');
    title('Upwind');
    legend show
    grid on
    hold off
    end

    % Soluzione SUPG.
    u_h = assemble_SUPG(epsi,beta,f,h,N(k),Pe);
    u_h = [0;u_h;0];
    err_SUPG(k) = sqrt(h * sum((u(x) - u_h').^2))
    if plotting == 1
    figure
    hold on
    plot(x, u(x), 'g', 'DisplayName', 'Exact Solution');
    plot(x, u_h, 'o-', 'DisplayName', ['SUPG, Pe = ', num2str(Pe)]);
    xlabel('x');
    ylabel('u(x)');
    title('SUPG');
    legend show
    grid on
    hold off
    end
    end
end
