function [err, err_H1] = errore_per_P2(u_esatta, u_star, geom, gN, u0, uD)

    errore_L2 = zeros(size(u_star,1),1);
    errore_H1 = zeros(size(u_star,1),1);
    
    % Punti di quadratura e pesi
    [zita, csi, eta, omega] = int_nodes_weights(5, 'triangle');
    
    % Calcolo l'errore in norma L^2
    for e = 1:geom.Nobj.N_ele
        area = geom.support.TInfo(e).Area;
        
        v1 = geom.obj.T(e,1);
        v2 = geom.obj.T(e,2);
        v3 = geom.obj.T(e,3);
    
        x1 = geom.obj.P(v1,1); y1 = geom.obj.P(v1,2);
        x2 = geom.obj.P(v2,1); y2 = geom.obj.P(v2,2);
        x3 = geom.obj.P(v3,1); y3 = geom.obj.P(v3,2);
    
        x_quad = zita * x3 + csi * x1 + eta * x2;
        y_quad = zita * y3 + csi * y1 + eta * y2;
    
        % Approssimazione dell'errore quadrato su T
        error_squared_T = 0;
        seminorma_H1_squared_T = 0;
    
        for q = 1:length(omega)
    
            grad_phi_local = [
                4*csi(q) - 1, 0;
                0, 4*eta(q) - 1;
                4*(csi(q)+eta(q)) - 3, 4*(csi(q)+eta(q)) - 3;
                4*eta(q), 4*csi(q);
                -4*eta(q), 4*(1 - csi(q) - 2*eta(q));
                4*(1 - 2*csi(q) - eta(q)), -4*csi(q)
                ];
            phi_local = [2*csi(q)*(csi(q) - 1/2);
                2*eta(q)*(eta(q) - 1/2);
                2*(1 - csi(q) - eta(q))*((1 - csi(q) - eta(q)) - 1/2);
                4*eta(q)*csi(q);
                4*eta(q)*(1 - csi(q) - eta(q));
                4*(1 - csi(q) - eta(q))*csi(q)];
            
            u_h = 0;
            for i = 1 : 6
                if geom.piv.piv(geom.obj.T(e,i)) > 0
                    u_h = u_h + u0(geom.piv.piv(geom.obj.T(e,i))) * phi_local(i);
                elseif geom.piv.piv(geom.obj.T(e,i)) < 0
                    u_h = u_h + uD(-geom.piv.piv(geom.obj.T(e,i))) * phi_local(i);
                end
            end
    
            % Errore in norma L2 al quadrato nel punto di quadratura
            integrand = (u_esatta(x_quad(q),y_quad(q)) - u_h)^2;
            error_squared_T = error_squared_T + omega(q) * integrand;
    
            % Errore in seminorma H1 nel punto di quadratura
            B_T = 1/(2*area)*[(y2-y3) (y3-y1); (x3-x2) (x1-x3)];
            grad_uh = 0;
            for i = 1:6
                if geom.piv.piv(geom.obj.T(e,i)) > 0
                    grad_uh = grad_uh + u0(geom.piv.piv(geom.obj.T(e,i))) * B_T * grad_phi_local(i,:)';
                else
                    grad_uh = grad_uh + uD(-geom.piv.piv(geom.obj.T(e,i))) * B_T * grad_phi_local(i,:)';
                end
            end
    
            integrand_H1 = (gN(x_quad(q),y_quad(q)) - grad_uh)' * (gN(x_quad(q),y_quad(q)) - grad_uh);
            seminorma_H1_squared_T = seminorma_H1_squared_T + omega(q) * integrand_H1;
        end %q

        errore_L2(e) = 2 * area * error_squared_T;
        errore_H1(e) = 2 * area * seminorma_H1_squared_T;

    end %e
    
    err = sum(errore_L2);
    err = sqrt(err);
    
    
    err_H1 = (sum(errore_H1));
    err_H1 = sqrt(err_H1);
  
end


