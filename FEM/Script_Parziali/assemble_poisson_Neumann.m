function [A,b] = assemble_poisson_Neumann(geom,nu,f,gn)
E = geom.obj.T;
P = geom.obj.P;
Piv = geom.piv.piv;
N = geom.Nobj.N_node-size(geom.piv.Ne,1)-size(geom.piv.Di,1);
A = zeros(N,N);
b = zeros(N,1);
run("sample_square_dirichlet_short.m");
for e = 1:length(E)
    d_x1 = P(E(e,3),1) - P(E(e,2),1);
    d_x2 = P(E(e,1),1) - P(E(e,3),1);
    d_x3 = P(E(e,2),1) - P(E(e,1),1);
    d_y1 = P(E(e,2),2) - P(E(e,3),2);
    d_y2 = P(E(e,3),2) - P(E(e,1),2);
    d_y3 = P(E(e,1),2) - P(E(e,2),2);
    d_x = [d_x1,d_x2,d_x3];
    d_y = [d_y1,d_y2,d_y3];
    Area = geom.support.TInfo(e).Area; 
    x_bar = (geom.obj.P(E(e,1),1) + geom.obj.P(E(e,2),1) + geom.obj.P(E(e,3),1)) / 3;
    y_bar = (geom.obj.P(E(e,1),2) + geom.obj.P(E(e,2),2) + geom.obj.P(E(e,3),2)) / 3;
    for j = 1:3
        jj = Piv(E(e,j));
        if jj > 0
            for k = 1:3
                kk = Piv(E(e,k));
                if kk > 0
                    djk = ((nu')/(4*Area))*(d_y(k)*d_y(j) +d_x(k)*d_x(j));
                    %cjk = (1/6)*(beta(1)*d_y(k) + beta(2)*d_x(k));
                    %rjk = ((sigma*Area)/12)*(1+j==k);
                    A(jj,kk) = A(jj,kk) +djk; %+cjk +rjk
                end
            end
            b(jj) = b(jj) + (1/3)*(Area*f(x_bar,y_bar));
        end
    end
end
% Ciclo che aggiunge i contributi di Dirichlet (bn)
for e = 1: size(geom.piv.Ne)
    l = geom.piv.Ne(e);
    ni = geom.obj.E(l,1);
    nf = geom.obj.E(l,2);
    xi = Piv(ni,1);
    xf = Piv(nf,1);
    yi = Piv(ni,2);
    yf = Piv(nf,2);
    ji = Piv(ni);
    jf = Piv(nf);
    E = sqrt((xf-xi)^2 + (yf-yi)^2);
    if ji > 0
        b(ji) = b(ji) + (1/6)*E*(2*gd(xi,yi)+gd(xf,yf));
    end
    if jf > 0
        b(jf) = b(jf) +(1/6)*E*(gd(xi,yi)+2*gd(xf,yf));
    end       
end