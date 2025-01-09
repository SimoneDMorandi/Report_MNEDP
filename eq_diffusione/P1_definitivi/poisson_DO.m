function [u_star, u0] = poisson_DO(geom, nu, f)

    % Numero di nodi effettivamente da considerare
    N = geom.Nobj.N_node - size(geom.piv.Di,1) - size(geom.piv.Ne,1);
    A = zeros(N, N);
    b = zeros(N, 1);

    for e = 1 : geom.Nobj.N_ele

       area = geom.support.TInfo(e).Area;

       % Calcolo dei vertici
       v1 = geom.obj.T(e,1);
       v2 = geom.obj.T(e,2);
       v3 = geom.obj.T(e,3);

       % Baricentro del triangolo
       x_bar = (geom.obj.P(v1,1) + geom.obj.P(v2,1) + geom.obj.P(v3,1)) / 3;
       y_bar = (geom.obj.P(v1,2) + geom.obj.P(v2,2) + geom.obj.P(v3,2)) / 3;

       % Calcolo delle delta 
       delta_x = zeros(3,1);
       delta_y = zeros(3,1);
       delta_x(1) = geom.obj.P(v3,1) - geom.obj.P(v2,1);
       delta_x(2) = geom.obj.P(v1,1) - geom.obj.P(v3,1);
       delta_x(3) = geom.obj.P(v2,1) - geom.obj.P(v1,1);
       delta_y(1) = geom.obj.P(v2,2) - geom.obj.P(v3,2);
       delta_y(2) = geom.obj.P(v3,2) - geom.obj.P(v1,2);
       delta_y(3) = geom.obj.P(v1,2) - geom.obj.P(v2,2);

       for j = 1 : 3

        jj = geom.piv.piv(geom.obj.T(e,j));

        if (jj > 0)
            for k = 1 : 3
                kk = geom.piv.piv(geom.obj.T(e,k));
                if (kk > 0)
                A(jj,kk) = A(jj,kk) + nu * (delta_y(k) * delta_y(j) + delta_x(k) * delta_x(j))/(4*area);
                end % if su kk 
                
            end % for su k

            % Valuto f nel baricentro del triangolo
            b(jj) = b(jj) + area * f(x_bar, y_bar) / 3;

        end % if su jj

      end %for su j

    end % for su e

    u0 = gmres(A,b,[],1e-12,size(A,1));

    % Impongo le condizioni di Dirichlet al bordo nella posizione corretta
    counter = 1;
    u_star = zeros(geom.Nobj.N_node,1);
    for i = 1:geom.Nobj.N_node
        if geom.piv.piv(i) > 0
            u_star(i) = u0(counter);
            counter = counter + 1;
        end
    end

end