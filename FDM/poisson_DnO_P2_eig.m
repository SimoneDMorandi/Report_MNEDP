function [u_star, u0, uD,A,b,M] = poisson_DnO_P2_eig(geom, nu, f, gD, gN)

% Risolve l'equazione di diffusione con anche la matrice di massa a parte.
% Punti di quadratura e pesi
[zita, csi, eta, omega] = int_nodes_weights(5, 'triangle');
nq = length(omega);

N = geom.Nobj.N_node + geom.Nobj.N_edge - size(geom.piv.Di,1);
M = zeros(N,N);
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

    x1 = geom.obj.P(v1,1); y1 = geom.obj.P(v1,2);
    x2 = geom.obj.P(v2,1); y2 = geom.obj.P(v2,2);
    x3 = geom.obj.P(v3,1); y3 = geom.obj.P(v3,2);

    % Trasformazione baricentrica -> globale
    x_quad = zita * x3 + csi * x1 + eta * x2;
    y_quad = zita * y3 + csi * y1 + eta * y2;

    % Calcolo delle delta
    delta_x = zeros(3,1);
    delta_y = zeros(3,1);
    delta_x(1) = x3 - x2;
    delta_x(2) = x1 - x3;
    delta_x(3) = x2 - x1;
    delta_y(1) = y2 - y3;
    delta_y(2) = y3 - y1;
    delta_y(3) = y1 - y2;

    B_1 = 1/(2*area) * [delta_y(1) delta_x(1); delta_y(2) delta_x(2)]; % B^-1

    % Matrice e termine noto locale
    A_local = zeros(6, 6);
    b_local = zeros(6,1);
    M_local = zeros(6,6); % Aggiungo la matrice di massa.
    for q = 1:nq
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
        for l = 1:6
            b_local(l) = b_local(l) + omega(q) * phi_local(l) * f(x_quad(q),y_quad(q)) * (2 * area);
            for k = 1:6
                A_local(l, k) = A_local(l, k) + ...
                    omega(q) * grad_phi_local(k,:) * B_1 * B_1' * grad_phi_local(l,:)' * (2 * area);
                M_local(l,k) = M_local(l, k) + ...
                omega(q) * phi_local(l) * phi_local(k) * (2 * area); % Aggiungo la matrice di massa.
            end
        end
    end

    % Assemblaggio della matrice globale
    for l = 1 : 6

        jj = geom.piv.piv(geom.obj.T(e,l));

        if (jj > 0)
            for k = 1 : 6
                kk = geom.piv.piv(geom.obj.T(e,k));
                if (kk > 0)
                    A(jj,kk) = A(jj,kk) + nu * A_local(l,k);
                    M(jj, kk) = M(jj, kk) + M_local(l, k); % Aggiungo la matrice di massa.
                else
                    AD(jj, -kk) = AD(jj, -kk) +  nu * A_local(l,k);
                end % if su kk

            end % for su k

            b(jj) = b(jj) + b_local(l);

        end % if su jj

    end %for su l

end % for su e

for i = 1:size(geom.piv.Di,1)
    n = geom.piv.Di(i,1);
    uD(i) = gD(geom.obj.P(n,1),geom.obj.P(n,2));
end

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
    n = zeros(2, 1);
    switch marker
        case 2, n = [0; -1];
        case 4, n = [1; 0];
        case 6, n = [0; 1];
        case 8, n = [-1; 0];
    end

    % Punti di quadratura e pesi per il lato [0,1]
    xi_q = [0.1127016654, 0.5, 0.8872983346]; % Punti di quadratura
    omega_lato = [5/18, 8/18, 5/18]; % Pesi di quadratura

    % Contributo al vettore b
    contributo = zeros(3, 1); % Per phi_1, phi_2, phi_3

    for k = 1:3

        xq = (1 - xi_q(k)) * coord_v1 + xi_q(k) * coord_v2;

        phi_local = [
            2 * (1 - xi_q(k)) * ((1 - xi_q(k)) - 0.5);
            4 * xi_q(k) * (1 - xi_q(k));
            2 * xi_q(k) * (xi_q(k) - 0.5)
        ];

        for j_index = 1:3
            contributo(j_index) = contributo(j_index) + omega_lato(k) * gN(xq(1), xq(2))' * n * phi_local(j_index);
        end
    end

    contributo = L * contributo;

    % Aggiorna il vettore dei carichi
    if dof_v1 > 0, b(dof_v1) = b(dof_v1) + contributo(1); end
    if dof_vm > 0, b(dof_vm) = b(dof_vm) + contributo(2); end
    if dof_v2 > 0, b(dof_v2) = b(dof_v2) + contributo(3); end
end % e

u0 = gmres(A, b - AD*uD,[],1e-12,size(A,1));

u_star = zeros(N,1);

% Impongo le condizioni di Dirichlet non omogeneo al bordo
for i = 1:N
    ii = geom.piv.piv(i);
    if ii > 0
        u_star(i) = u0(ii);
    else
        u_star(i) = uD(-ii);
    end
end

end