clear
clc
%% Equazione del calore con Dirichlet non omogeneo

% Dati iniziali.
% Il triangolatore fornisce le matrici per costruire il sistema lineare.
run("sample_square_dirichlet_short.m");
u = @(x,y) x+y +16*x.*(1-x).*y.*(1-y);
f = @(x,y) 32*((y-y^2)+(x-x^2));
nu = 1.0;
gd =  @(x,y) x+y;

% function che costruisce il sistema lineare.
[A,b,ud] = assemble_poisson_No_D(geom,nu,f,gd);
% Risoluzione del sistema lineare.
u_h = A\b;


% Creazione del vettore per rappresentazione grafica.
u_star= zeros(geom.Nobj.N_node,1);
for i = 1:geom.Nobj.N_node
    ii = geom.piv.piv(i);
    if ii > 0
        u_star(i) = u_h(ii);
    else
        u_star(i) = ud(-ii);
    end
end
% Rappresentazione Grafica.
figure
trisurf(geom.obj.T,geom.obj.P(:,1),geom.obj.P(:,2),u_star)
figure
trisurf(geom.obj.T,geom.obj.P(:,1), geom.obj.P(:,2),u(geom.obj.P(:,1), geom.obj.P(:,2)))