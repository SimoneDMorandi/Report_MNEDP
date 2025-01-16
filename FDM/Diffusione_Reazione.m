clear
clc
format short e
% Solving -epsilon*u'' + sigma*u = f
%          u(0) = 0, u(L) = 0;
%          f = (sigma-espilon)*exp(x);
plotting_f = 0;
epsilon = 0.001;
L = 1;
sigma = 500;
N = 50;
err_FDM = zeros(length(N),length(epsilon));
err_FEM = zeros(length(N),length(epsilon));
err_M = zeros(length(N),length(epsilon));
for j = 1:length(epsilon)
epsi = epsilon(j);
alpha = sqrt(sigma/epsi);
f = @(x) (sigma-epsi)*exp(x);
%u_hom = @(x) (exp(-alpha*x)-exp(-alpha*x))/(exp(-alpha)-exp(alpha));
u = @(x) exp(x)+ ((-exp(-alpha)-exp(1))/(exp(alpha)-exp(-alpha)))*exp(alpha*x) +((-exp(alpha)-exp(1))/(exp(alpha)-exp(-alpha)))*exp(-alpha*x);
    h = L/(N+1);
    Se = (sigma*(h^2))/(epsi)
    x = linspace(0,L,N+2);
    % FDM solution.
    main_diag = -(2+Se)*ones(N,1);               
    lower_diag = ones(N-1,1);  
    upper_diag = ones(N-1,1);
    A = diag(main_diag,0)+diag(lower_diag,-1)+diag(upper_diag,1);
    b = -((h^2)/epsi)*f(x(2:end-1))';
    u_h = A\b;
    err_FDM = sqrt(h * sum((u(x(2:end-1)) - u_h').^2)) 
    if plotting_f == 1
        figure
        hold on;
        plot(x(2:end-1), u(x(2:end-1)), 'g', 'DisplayName', 'Exact Solution');
        plot(x(2:end-1), u_h, 'o-', 'DisplayName', ['Finite Difference, Se = ', num2str(Se)]);
        xlabel('x');
        ylabel('u(x)');
        title('Centered Finite Difference');
        grid on;
        legend show;
    end

    % FEM solution.
    Sef = ((sigma*(h^2))/(6*epsi))
    main_diag = (2+4*Sef)*ones(N,1);               
    lower_diag = (Sef-1)*ones(N-1,1);  
    upper_diag = (Sef-1)*ones(N-1,1);
    A = diag(main_diag,0)+diag(lower_diag,-1)+diag(upper_diag,1);
    b = ((h^2)/epsi)*f(x(2:end-1))';
    u_h = A\b;
    err_FEM = sqrt(h * sum((u(x(2:end-1)) - u_h').^2)) 
    if plotting_f == 1
        figure
        hold on;
        plot(x(2:end-1), u(x(2:end-1)), 'g', 'DisplayName', 'Exact Solution');
        plot(x(2:end-1), u_h, 'o-', 'DisplayName', ['Finite Difference, Sef = ', num2str(Sef)]);
        xlabel('x');
        ylabel('u(x)');
        title('Finite Elements');
        grid on;
        legend show;
    end
    
    % Mass lumping.
    r_0 = (2/3)*ones(N,1);
    r_1 = (1/6)*ones(N-1,1);
    R = (sigma*h)*(diag(r_0)+diag(r_1,1)+diag(r_1,-1));
    M = zeros(N,N); 
    for i = 1:N
        M(i,i) = sum(R(i,:));
    end
    main_diag = (6*ones(N-2,1))';
    main_diag = Sef*[5,main_diag,5] + 2*ones(N,1)';             
    lower_diag = -ones(N-1,1);  
    upper_diag = -ones(N-1,1);
    A = diag(main_diag,0)+diag(lower_diag,-1)+diag(upper_diag,1);
    u_h = A\b;
    err_Mass = sqrt(h * sum((u(x(2:end-1)) - u_h').^2)) 
    if plotting_f == 1
        figure
        hold on;
        plot(x(2:end-1), u(x(2:end-1)), 'g', 'DisplayName', 'Exact Solution');
        plot(x(2:end-1), u_h, 'o-', 'DisplayName', ['Mass Lumping Solution, Sef = ', num2str(Sef)]);
        xlabel('x');
        ylabel('u(x)');
        title('Mass Lumping');
        grid on;
        legend show;
    end
end
