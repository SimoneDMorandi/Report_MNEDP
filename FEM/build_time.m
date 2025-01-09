function [F,ud,B,Bd,K,Kd] = build_time(geom,nu,f,ud,gd,t)
E = geom.obj.T;
P = geom.obj.P;
Piv = geom.piv.piv;
N = geom.Nobj.N_node -size(geom.piv.Di,1);
Nd = size(geom.piv.Di,1);
B = zeros(N,N);
K = zeros(N,N);
Bd = zeros(N,Nd);
Kd = zeros(N,Nd);
F = zeros(N,1);
for e = 1:length(E)
    d_x1 = P(E(e,3),1)-P(E(e,2),1);
    d_x2 = P(E(e,1),1)-P(E(e,3),1);
    d_x3 = P(E(e,2),1)-P(E(e,1),1);
    d_y1 = P(E(e,2),2)-P(E(e,3),2);
    d_y2 = P(E(e,3),2)-P(E(e,1),2);
    d_y3 = P(E(e,1),2)-P(E(e,2),2);
    d_x = [d_x1,d_x2,d_x3];
    d_y = [d_y1,d_y2,d_y3];
    Area = geom.support.TInfo(e).Area;
    x_bar = (P(E(e,1),1) + P(E(e,2),1) + P(E(e,3),1))/3;
    y_bar = (P(E(e,1),2) + P(E(e,2),2) + P(E(e,3),2))/3;
    for j = 1:3
        jj = Piv(E(e,j));
        if jj > 0
            for k = 1:3
                kk = Piv(E(e,k));
                if kk > 0
                    Kjk = ((nu')/(4*Area))*(d_y(k)*d_y(j) +d_x(k)*d_x(j));
                    bjk = (Area/12)*(1+ j==k);  
                    B(jj,kk) = B(jj,kk) +bjk;
                    K(jj,kk) = K(jj,kk) + Kjk;
                else
                    Kd(jj,-kk) = Kd(jj,-kk) +((nu')/(4*Area))*(d_y(k)*d_y(j)+d_x(k)*d_x(j));
                    Bd(jj,-kk) = B(jj,-kk) + (Area/12)*(1+j==k);
                end
            end
            F(jj) = F(jj) + (1/3)*(Area*f(x_bar,y_bar,t));
        end
    end
end
for i = 1:size(geom.piv.Di,1)
    n = geom.piv.Di(i,1);
    ud(i) = gd(P(n,1),P(n,2),t);
end
end