function [u_cn, u, uD] = build_solution_cn(geom, A, B, AD, BD, gN, gD, f, u0, n_t, delta_t, u_esatta)

% Punti di quadratura e pesi
[zita, csi, eta, omega] = int_nodes_weights(5, 'triangle');
nq = length(omega);

N = geom.Nobj.N_node + geom.Nobj.N_edge - size(geom.piv.Di,1);
u_cn = zeros(N,n_t+1);
u = zeros(N,n_t+1);
uD = zeros(size(geom.piv.Di,1),n_t);

% Avanzamento in tempo
u_cn(:,1) = u0(geom.obj.P(geom.piv.piv > 0, 1), geom.obj.P(geom.piv.piv > 0, 2));
u(:,1) = u_esatta(geom.obj.P(geom.piv.piv > 0, 1), geom.obj.P(geom.piv.piv > 0, 2), 0);
for n = 1 : n_t

    uD_0 = zeros(size(geom.piv.Di,1),1);
    uD_1 = zeros(size(geom.piv.Di,1),1);

    % Calcolo la forzante al tempo tn e al tempo tn+1
    F_0 = zeros(N, 1);
    F_1 = zeros(N, 1);
    for e = 1 : geom.Nobj.N_ele
        area = geom.support.TInfo(e).Area;

        % Calcolo dei vertici
        v1 = geom.obj.T(e,1);
        v2 = geom.obj.T(e,2);
        v3 = geom.obj.T(e,3);

        x1 = geom.obj.P(v1,1); y1 = geom.obj.P(v1,2);
        x2 = geom.obj.P(v2,1); y2 = geom.obj.P(v2,2);
        x3 = geom.obj.P(v3,1); y3 = geom.obj.P(v3,2);

        % Trasformazione baricentrica -> globale
        x_quad = zita * x3 + csi * x1 + eta * x2;
        y_quad = zita * y3 + csi * y1 + eta * y2;

        F_local0 = zeros(6,1);
        F_local1 = zeros(6,1);

        for q = 1 : nq
            phi_local = [2*csi(q)*(csi(q) - 1/2);
            2*eta(q)*(eta(q) - 1/2);
            2*(1 - csi(q) - eta(q))*((1 - csi(q) - eta(q)) - 1/2);
            4*eta(q)*csi(q);
            4*eta(q)*(1 - csi(q) - eta(q));
            4*(1 - csi(q) - eta(q))*csi(q)];
            for l=1:6
                F_local0(l) = F_local0(l) + omega(q) * phi_local(l) * f(x_quad(q),y_quad(q),delta_t*(n-1)) * 2 * area;
                F_local1(l) = F_local1(l) + omega(q) * phi_local(l) * f(x_quad(q),y_quad(q),delta_t*(n)) * 2 * area;
            end % l
        end % q

        for l = 1 : 6
            jj = geom.piv.piv(geom.obj.T(e,l));
            if (jj > 0)
                F_0(jj) = F_0(jj) + F_local0(l);
                F_1(jj) = F_1(jj) + F_local1(l);
            end % if su jj
        end % for su l
    end % for su e (elementi della triangolazione)


    % Impongo le condizioni di Neumann
    for e = 1:size(geom.piv.Ne, 1)
        % Id lato con condizione di Neumann
        id_edge = geom.piv.Ne(e, 1);
        marker = geom.piv.Ne(e, 2);

        % Vertici che delimitano il lato
        v1 = geom.obj.E(id_edge, 1);
        v2 = geom.obj.E(id_edge, 2);
        vm = geom.obj.E(id_edge, end);

        % DOF e coordinate dei vertici
        dof_v1 = geom.piv.piv(v1);
        dof_v2 = geom.piv.piv(v2);
        dof_vm = geom.piv.piv(vm);
        coord_v1 = geom.obj.P(v1, :);
        coord_v2 = geom.obj.P(v2, :);
        coord_vm = geom.obj.P(vm, :);

        % Lunghezza del lato
        L = norm(coord_v2 - coord_v1, 2);

        % Normale al lato
        nn = zeros(2, 1);
        switch marker
            case 2, nn = [0; -1];
            case 4, nn = [1; 0];
            case 6, nn = [0; 1];
            case 8, nn = [-1; 0];
        end

        % Punti di quadratura e pesi per il lato [0,1]
        xi_q = [0.1127016654, 0.5, 0.8872983346]; % Punti di quadratura
        omega_lato = [5/18, 8/18, 5/18]; % Pesi di quadratura

        % Contributo al vettore b
        contributo_0 = zeros(3, 1); % Per phi_1, phi_2, phi_3
        contributo_1 = zeros(3, 1);

        for k = 1:3

            xq = (1 - xi_q(k)) * coord_v1 + xi_q(k) * coord_v2;

            phi_local = [
                2 * (1 - xi_q(k)) * ((1 - xi_q(k)) - 0.5);
                4 * xi_q(k) * (1 - xi_q(k));
                2 * xi_q(k) * (xi_q(k) - 0.5)
                ];

            for j_index = 1:3
                contributo_0(j_index) = contributo_0(j_index) + omega_lato(k) * gN(xq(1), xq(2), delta_t*(n-1))' * nn * phi_local(j_index);
                contributo_1(j_index) = contributo_1(j_index) + omega_lato(k) * gN(xq(1), xq(2), delta_t*(n))' * nn * phi_local(j_index);
            end
        end

        contributo_0 = L * contributo_0;
        contributo_1 = L * contributo_1;

        if dof_v1 > 0
            F_0(dof_v1) = F_0(dof_v1) + contributo_0(1);
            F_1(dof_v1) = F_1(dof_v1) + contributo_1(1);
        end
        if dof_vm > 0
            F_0(dof_vm) = F_0(dof_vm) + contributo_0(2);
            F_1(dof_vm) = F_1(dof_vm) + contributo_1(2);
        end
        if dof_v2 > 0
            F_0(dof_v2) = F_0(dof_v2) + contributo_0(3);
            F_1(dof_v2) = F_1(dof_v2) + contributo_1(3);
        end
    end % e


   

    for i = 1:size(geom.piv.Di,1)
        uD_0(i) =gD(geom.obj.P(geom.piv.Di(i,1),1),geom.obj.P(geom.piv.Di(i,1),2),delta_t * (n-1));
        uD_1(i) =gD(geom.obj.P(geom.piv.Di(i,1),1),geom.obj.P(geom.piv.Di(i,1),2),delta_t * (n));
    end

    % Crank-Nicolson
    u_cn(:,n+1) = gmres((B + delta_t/2 * A),(B * u_cn(:,n) + ...
        delta_t/2*(-A*u_cn(:,n) + F_0 - AD * uD_0 + F_1 - AD*uD_1)) ...
        - BD * (uD_1 - uD_0), [], 1e-12, size(B + delta_t/2 * A,1));
    u(:,n+1) = u_esatta(geom.obj.P(geom.piv.piv > 0, 1),geom.obj.P(geom.piv.piv > 0, 2),delta_t*(n));
    uD(:,n) = uD_0;

    if (n == n_t), uD(:,n+1) = uD_1; 

end %n


end
