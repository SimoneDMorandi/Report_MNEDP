function [u_h] = assemble_SUPG(epsi,beta,f,h,N,Pe)
x = linspace(0,1,N-1)';
d_0 =  2*ones(N-1,1);
d_1 = -ones(N-2,1);
D = (epsi/h)*(diag(d_0)+diag(d_1,1)+diag(d_1,-1));

c_1 = ones(N-2,1);
C = (beta/2)*(diag(c_1,1)+diag(-c_1,-1));
A = D+C+R;
b = h*(f(x));
if Pe <= 1
    tau = (1/3)*((h^2)/(4*epsi));
else
    tau = h/(2*beta);
end
T = tau*(beta^2)*D
A = A + T
% Risoluzione del sistema lineare.
u_h = A\b