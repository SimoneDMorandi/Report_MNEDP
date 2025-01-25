 function [u_star, u0, uD] = poisson_DnO(geom, nu, f, gD, gN)
    
    % Numero di nodi effettivamente da considerare
    N = geom.Nobj.N_node - size(geom.piv.Di,1);
    A = zeros(N, N);
    AD = zeros(N, size(geom.piv.Di,1));
    b = zeros(N, 1);
    uD = zeros(size(geom.piv.Di,1),1);

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
                else
                    AD(jj, -kk) = AD(jj, -kk) +  nu * (delta_y(k) * delta_y(j) + delta_x(k) * delta_x(j))/(4*area);
                end % if su kk 
                
            end % for su k

            % Valuto la forzante nel baricentro del triangolo
            b(jj) = b(jj) + area * f(x_bar, y_bar) / 3;

        end % if su jj

      end %for su j

    end % for su e

    for i = 1:size(geom.piv.Di,1)
        n = geom.piv.Di(i,1);
        uD(i) = gD(geom.obj.P(n,1),geom.obj.P(n,2));
    end

    % Impongo le condizioni di Neumann
    for e = 1 : size(geom.piv.Ne,1)

        % Id lato con condizione di Neumann
        id_edge = geom.piv.Ne(e,1);
        marker = geom.piv.Ne(e,2);

        % Vertici che delimitano il lato
        v1 = geom.obj.E(id_edge,1);
        v2 = geom.obj.E(id_edge,2);

        % Dof e coordinate dei vertici
        dof_v1 = geom.piv.piv(v1);
        dof_v2 = geom.piv.piv(v2);

        % Coordinate dei vertici
        coord_v1 = geom.obj.P(v1,:);
        coord_v2 = geom.obj.P(v2,:);

        L = norm(coord_v2 - coord_v1, 2);
        n = zeros(2,1);
        
        % Impongo il valore del versore normale
        switch marker
            case 2
                n = [0;-1];

            case 4
                n = [1;0];

            case 6
                n = [0;1];

            case 8
                n = [-1;0];
        end % switch

        if dof_v1 > 0
            b(dof_v1) = b(dof_v1) + L/2* gN(coord_v1(1), coord_v1(2))'*n;
        end
        if dof_v2 > 0
            b(dof_v2) = b(dof_v2) + L/2* gN(coord_v2(1), coord_v2(2))'*n;
        end
            
    end % for e

    u0 = gmres(A, b - AD*uD,[],1e-12,size(A,1));

    u_star = zeros(size(geom.obj.P,1),1);

    % Impongo le condizioni di Dirichlet non omogeneo al bordo
    for i = 1:geom.Nobj.N_node
        ii = geom.piv.piv(i);
        if ii > 0
            u_star(i) = u0(ii);
        else
            u_star(i) = uD(-ii);
        end
    end

end