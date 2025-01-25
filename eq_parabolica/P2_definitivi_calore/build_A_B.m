function [N, A, B, AD, BD] = build_A_B(geom, nu)

% Punti di quadratura e pesi
[zita, csi, eta, omega] = int_nodes_weights(5, 'triangle');
nq = length(omega);

N = geom.Nobj.N_node + geom.Nobj.N_edge - size(geom.piv.Di,1);
A = zeros(N, N);
B = zeros(N, N);
AD = zeros(N, size(geom.piv.Di,1));
BD = zeros(N, size(geom.piv.Di,1));

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

    % Assemblo la matrice di rigidezza e la matrice di massa locali
    A_local = zeros(6,6);
    B_local = zeros(6,6);
    for q = 1 : nq

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

        for l=1:6
            for k = 1:6
                A_local(l,k) = A_local(l,k) +  ...
                    omega(q) * grad_phi_local(k,:) * B_1 * B_1' * grad_phi_local(l,:)' * (2 * area);
                B_local(l,k) = B_local(l,k) + omega(q) * phi_local(k) * phi_local(l)' * 2* area;
            end
        end
    end

    % Assemblo le matrici globali
    for j = 1 : 6

        jj = geom.piv.piv(geom.obj.T(e,j));

        if (jj > 0)
            for k = 1 : 6
                kk = geom.piv.piv(geom.obj.T(e,k));
                if (kk > 0)
                    A(jj,kk) = A(jj,kk) + nu * A_local(j,k);
                    B(jj,kk) = B(jj,kk) + B_local(j,k);
                else
                    AD(jj, -kk) = AD(jj, -kk) +  nu * A_local(j,k);
                    BD(jj, -kk) = BD(jj, -kk) + B_local(j,k);
                end % if su kk

            end % for su k

        end % if su jj

    end %for su j

end % for su e

end
