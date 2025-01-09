clear
clc

% Risoluzione su 2Lx2L.
L = 1;
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

% Escludo i vertici del quadrato non presente nel dominio omega.
exclude = (X>=L) & (Y>=L);
keep = ~exclude(:);
A = A(keep,keep);
b = b_2(keep);
u = A\b;

% Ripristino il dominio iniziale con i valori annullati.
u_sol = zeros(N*N,1);
u_sol(keep) = u;
u_sol = reshape(u_sol,[N,N]);

% Rappresentazione grafica di u.
figure
surf(X,Y,u_sol);


% Ricerca del primo autovalore e autofunzione.
[V,D] = eigs(A,1,'smallestabs');
lambda = D(1,1);
u_eig = V(:,1);

% Rappresentazione grafica
u_sol_eig = zeros(N*N,1);
u_sol_eig(keep) = u_eig;
u_sol_eig = reshape(u_sol_eig,[N,N]);
figure
s = surf(X,Y,-u_sol_eig);
grid off
axis off
s.EdgeColor = 'none';
s.EdgeColor = 'black';