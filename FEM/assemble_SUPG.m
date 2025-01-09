function [A,b] = assemble_SUPG(geom,eps,beta,f)
E = geom.obj.T;
P = geom.obj.P;
Piv = geom.piv.piv;
N = geom.Nobj.N_node -size(geom.piv.Di,1);
A = zeros(N,N);
b = zeros(N,1);
beta = [beta,beta];
for e = 1:length(E)
    Area = geom.support.TInfo(e).Area;
    x1 = P(E(e,3),1) - P(E(e,2),1);
    x2 = P(E(e,1),1) - P(E(e,3),1);
    x3 = P(E(e,2),1) - P(E(e,1),1);
    y1 = P(E(e,2),2) - P(E(e,3),2);
    y2 = P(E(e,3),2) - P(E(e,1),2);
    y3 = P(E(e,1),2) - P(E(e,2),2);
    B_t = (1/(2*Area))*[y2-y3,y3-y1;x3-x2,x1-x3]
    phi_l_1 = @(x,y) x + beta*(x);
    phi_l_2 = @(x,y) y+beta*(y);
    phi_l_3 = @(x,y) 1-x-y + beta*(-x-y);
    xbar = (P(E(e,1),1) + P(E(e,2),1) + P(E(e,3),1))/3;
    ybar = (P(E(e,1),2) + P(E(e,2),2) + P(E(e,3),2))/3;
    for j = 1:3
        jj = Piv(E(e,j));
        if jj > 0
            for k = 1:3
                kk = Piv(E(e,k));
                if kk > 0
                    A(jj,kk) =sum( (beta*B_t*[1;0] + beta*B_t*[0;1] + beta*B_t*[-1;-1])* ...
                    (beta*B_t.*(phi_l_1(xbar,ybar)) +beta*B_t.*(phi_l_2(xbar,ybar))+ ...
                     beta*B_t.*(phi_l_3(xbar,ybar)))*2*Area);
                end
            b(jj) = b(jj); %(1/3)*(Area*f(xbar,ybar));
            end
        end
    end
end