function [E0,E1] = calculate_err(geom,u,u_s,gn,du_dx,du_dy)
E = geom.obj.T;
P = geom.obj.P;
run("nodes_weights.m")

E0 = zeros(size(u_s,1),1);
E1 = zeros(size(u_s,1),1);

% Errore in L2.
for e = 1:geom.Nobj.N_ele
    Area = geom.support.TInfo(e).Area;
    v1 = E(e,1);
    v2 = E(e,2);
    v3 = E(e,3);

    x1 = P(v1,1); y1 = P(v1,2);
    x2 = P(v2,1); y2 = P(v2,2);
    x3 = P(v3,1); y3 = P(v3,2);

    err_2_T = 0;
    semi_norm_H1_T = 0;

    for q = 1:length(omega)
        % Coord baricentriche.
        lambda1 = xhat(q);
        lambda2 = yhat(q);
        lambda3 = 1- lambda1 -lambda2;

        % Coord cartesiane.
        xq = lambda1*x1 + lambda2*x2 +lambda3*x3;
        yq = lambda1*y1 + lambda2*y2 + lambda3*y3;
        u_h = u_s(v1)*lambda1 + u_s(v2)*lambda2 + u_s(v3)*lambda3;

        % Aggiungo all'errore L2.
        integrale = (u(xq,yq)-u_h)^2;
        err_2_T = err_2_T + omega(q)*integrale;
        
        % Aggiungo alla seminorma in H1
        B_t = (1/(2*Area))*[y2-y3,y3-y1; x3-x2,x1-x3];
        grad_u = [du_dx(xq, yq); du_dy(xq, yq)];
        grad_uh = u_s(v1)*lambda1*B_t*[1;0] + u_s(v2)*lambda2*B_t*[0;1] +...
            + u_s(v3)*lambda3*B_t*[-1;-1];
        integrale_H1 = (grad_u-grad_uh)'*(grad_u-grad_uh);
        semi_norm_H1_T = semi_norm_H1_T + omega(q)*integrale_H1;
    end
    E0(e) = 2*Area*err_2_T;
    E1(e) = 2*Area*semi_norm_H1_T;
end
E0 = sqrt(sum(E0));
E1 = sqrt((sum(E1)));
end