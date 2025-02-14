%% 
clear
clc

% Risoluzione sul quadrato [0,1]x[0,1] del problema laplacian(u) = f(x,y) con Dirichlet omogeneo.

% Soluzione con FDM.
L = 0.5;
N= 50;
h = (2*L)/(N+1);
mu = 1;
d = ones(N,1);
c = -ones(N-1,1);
f = @(x,y) x+y;
x = (1:N)'*h;
y = (1:N)'*h;
[X,Y] = meshgrid(x,y);
b_2 = zeros(N*N,1);
D = (mu/(h^2))*(diag(4*d,0)+diag(c,1)+diag(c,-1));
C = (mu/(h^2))*diag(-ones(N,1),0);
A = kron(eye(N), D);

% Creazione dei blocchi sottodiagonali.
sub_block = diag(ones(N-1, 1), -1);
super_block = diag(ones(N-1, 1), 1);
sub_matrix = kron(sub_block, C);
super_matrix = kron(super_block, C);
A = A + sub_matrix + super_matrix;
for m = 1:N
    for l = 1:N
        j = (m-1)*N + l;
        b_2(j) = f(x(l),y(m));
    end
end

% Escludo i vertici del quadrato non presente nel dominio a L.
exclude = (X>=L) & (Y>=L);
keep = ~exclude(:);
A = A(keep,keep);
b = b_2(keep);
u = A\b;

% Ripristino il dominio iniziale con i valori annullati.
u_sol = zeros(N*N,1);
u_sol(keep) = u;
u_sol = reshape(u_sol,[N,N]);

% Plot di u_sol.
figure
subplot(1,2,1)
surf(X,Y,u_sol);
xlabel('x');
ylabel('y');
zlabel('u(x)');
title('L-shaped Domain Solution with FDM');

% Ricerca del primo autovalore e autofunzione.
% Soluzione del problema A*u = lambda*u, con la stessa A di prima.
[V,D] = eigs(A,1,'smallestabs');
lambda = D(1,1);
u_eig = V(:,1);

% Plot.
u_sol_eig = zeros(N*N,1);
u_sol_eig(keep) = u_eig;
u_sol_eig = reshape(u_sol_eig,[N,N]);
subplot(1,2,2)
s = surf(X,Y,-u_sol_eig);
xlabel('x');
ylabel('y');
zlabel('u(x)');
title(sprintf('First Eigenfunction, \\lambda = %.4f', lambda));


%% FEM
 
% Soluzione con FEM, utilizzando i P2.
clear
clc

% Dati iniziali.
area_ref = 1e-3;
marker_triang = [1, 3, 5, 7];
f = @(x,y) x+y;
gD = @(x,y) 0;
gN = @(x,y) 0;
nu = 1.0;
% Risolvo su [0,1]x[0,1]
geom = triangulator_P_2_eig(area_ref, marker_triang);
N = geom.Nobj.N_node + geom.Nobj.N_edge - size(geom.piv.Di,1);    
[~,~, uD,A,b,M] = poisson_DnO_P2_eig(geom, nu, f, gD, gN); % Calcolo anche la matrice di massa M per dopo.

% Annullo i dof sul quadrante in alto a destra e calcolo la nuova soluzione.
    x_threshold = 0.5;
    y_threshold = 0.5;
    top_right_nodes = find(geom.obj.P(:,1) > x_threshold + eps & geom.obj.P(:,2) > y_threshold + eps);
    top_right_dofs = geom.piv.piv(top_right_nodes); 
    top_right_dofs = top_right_dofs(top_right_dofs > eps);
    for dof = top_right_dofs'
        A(dof, :) = 0;               
        A(:, dof) = 0;              
        A(dof, dof) = 1;             
        b(dof) = 0;
        M(dof, :) = 0;               
        M(:, dof) = 0;              
        M(dof, dof) = 0;             
    end
    u0 = gmres(A, b,[],1e-12,size(A,1));
    u_star = zeros(N,1);
    for i = 1:N
        ii = geom.piv.piv(i);
        if ii > 0
            u_star(i) = u0(ii);
        else
            u_star(i) = uD(-ii);
        end
    end

% Plot.
figure
subplot(1,2,1);
s = trisurf(geom.obj.T(:,1:3), geom.obj.P(1:N,1), geom.obj.P(1:N,2), u_star);
xlabel('X');
ylabel('Y');
zlabel('u(x,y)');
title('L-shaped domain solution with FEM');
s.EdgeColor = 'Black';

% Calcolo autovettore e autofunzione.
% Risoluzione del problema A*u = lambda*M*u, con M matrice di massa e Dirichlet omogeneo.
[V,Lambda] = eigs(A,M,1,'smallestabs');
lambda = Lambda(1,1);
u_eig =V(:,1);
u_eig_star = zeros(N,1);
for i = 1:N
    ii = geom.piv.piv(i);
    if ii > 0
        u_eig_star(i) = u_eig(ii);
    end
end

% Plot.
subplot(1,2,2)
s = trisurf(geom.obj.T(:,1:3), geom.obj.P(1:N,1), geom.obj.P(1:N,2), u_eig_star, 'EdgeColor', 'none');
xlabel('X');
ylabel('Y');
zlabel('u(x,y)');
s.EdgeColor = 'Black';
title(sprintf('First Eigenfunction, \\lambda = %.4f', lambda));
