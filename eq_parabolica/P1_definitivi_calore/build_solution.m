function [u_ei, u_cn, u, uD] = build_solution(geom, A, B, AD, BD, gN, gD, f, u0, n_t, delta_t, u_esatta)

% Punti di quadratura e pesi
[zita, csi, eta, omega] = int_nodes_weights(5, 'triangle');
nq = length(omega);

N = geom.Nobj.N_node - size(geom.piv.Di,1); % dof
u_ei = zeros(N,n_t+1);
u_cn = zeros(N,n_t+1);
u = zeros(N,n_t+1);
uD = zeros(size(geom.piv.Di,1),n_t);

% Avanzamento in tempo
u_ei(:,1) = u0(geom.obj.P(geom.piv.piv > 0, 1), geom.obj.P(geom.piv.piv > 0, 2));
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

        F_local0 = zeros(3,1);
        F_local1 = zeros(3,1);

        for q = 1 : nq
            phi_local = [csi(q), eta(q), 1 - csi(q) - eta(q)];
            for l=1:3
                F_local0(l) = F_local0(l) + omega(q) * phi_local(l) * f(x_quad(q),y_quad(q),delta_t*(n-1)) * 2 * area;
                F_local1(l) = F_local1(l) + omega(q) * phi_local(l) * f(x_quad(q),y_quad(q),delta_t*(n)) * 2 * area;
            end % l
        end % q

        for l = 1 : 3
            jj = geom.piv.piv(geom.obj.T(e,l));
            if (jj > 0)
                F_0(jj) = F_0(jj) + F_local0(l);
                F_1(jj) = F_1(jj) + F_local1(l);
            end % if su jj
        end % for su l
    end % for su e (elementi della triangolazione)

    % Condizioni di Neumann
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
        nn = zeros(2,1);

        % Impongo il valore del versore normale
        switch marker
            case 2
                nn = [0;-1];

            case 4
                nn = [1;0];

            case 6
                nn = [0;1];

            case 8
                nn = [-1;0];
        end % switch

        if dof_v1 > 0
            F_0(dof_v1) = F_0(dof_v1) + L/2* gN(coord_v1(1), coord_v1(2),delta_t*(n-1))'*nn;
            F_1(dof_v1) = F_1(dof_v1) + L/2* gN(coord_v1(1), coord_v1(2),delta_t*(n))'*nn;
        end
        if dof_v2 > 0
            F_0(dof_v2) = F_0(dof_v2) + L/2* gN(coord_v2(1), coord_v2(2),delta_t*(n-1))'*nn;
            F_1(dof_v2) = F_1(dof_v2) + L/2* gN(coord_v2(1), coord_v2(2),delta_t*(n))'*nn;
        end

    end % for e

    for i = 1:size(geom.piv.Di,1)
        uD_0(i) =gD(geom.obj.P(geom.piv.Di(i,1),1),geom.obj.P(geom.piv.Di(i,1),2),delta_t * (n-1));
        uD_1(i) =gD(geom.obj.P(geom.piv.Di(i,1),1),geom.obj.P(geom.piv.Di(i,1),2),delta_t * (n));
    end

    % Eulero implicito
    u_ei(:,n+1) = gmres(B + delta_t * A, B * u_ei(:,n) + ...
        delta_t * (F_1 - AD * uD_1) - BD * uD_1 + BD * uD_0, [], 1e-12, size(B + delta_t * A,1));

    % Crank-Nicolson
    % u_cn(:,n+1) = gmres((B + delta_t/2 * A),(B * u_cn(:,n) + ...
    %     delta_t/2*(-A*u_cn(:,n) + F_0 - BD * (uD_1 - uD_0)/delta_t  ...
    %     -AD * uD_0 + F_1 - BD * (uD_1 - uD_0)/delta_t +...
    %     - AD*uD_1)), [], 1e-12, size(B + delta_t * A,1));
    u_cn(:,n+1) = gmres((B + delta_t/2 * A),(B * u_cn(:,n) + ...
        delta_t/2*(-A*u_cn(:,n) + F_0 - AD * uD_0 + F_1 - AD*uD_1)) ...
        - BD * (uD_1 - uD_0), [], 1e-12, size(B + delta_t/2 * A,1));
    u(:,n+1) = u_esatta(geom.obj.P(geom.piv.piv > 0, 1),geom.obj.P(geom.piv.piv > 0, 2),delta_t*(n));
    uD(:,n) = uD_0;
    if (n == n_t), uD(:,n+1) = uD_1;
end %n


end
