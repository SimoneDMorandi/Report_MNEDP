function [u_h] = assemble_SUPG(epsi,beta,f,h,N,Pe)
% Assembla la soluzione numerica per il problema SUPG in una dimensione con i P1.

% Definisco termine convettivo, diffusivo e termine noto.
x = linspace(0,1,N+2)';
d_0 =  2*ones(N,1);
d_1 = -ones(N-1,1);
D = (epsi/h)*(diag(d_0)+diag(d_1,1)+diag(d_1,-1));
c_1 = ones(N-1,1);
C = (beta/2)*(diag(c_1,1)+diag(-c_1,-1));
b = h*(f(x(2:N+1)));

% Costruisco la correzione con tau.
if Pe <= 1
    tau = (1/3)*((h^2)/(4*epsi));
else
    tau = h/(2*beta);
end
t_0 = (2*beta^2)*ones(N,1);
t_1 = -(beta^2)*ones(N-1,1);

% Aggiungo la correzione alla matrice del termine diffusivo e al termine noto.
D = D + tau*(beta^2 /h)*(diag(d_1,1)+diag(d_0,0) +diag(d_1,-1)); 
b = b -(tau*beta)*(f(x(2:N+1)));
A = D+C;

% Risoluzione del sistema lineare.
u_h = A\b;