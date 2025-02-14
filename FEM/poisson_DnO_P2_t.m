function [u_star, u0, uD] = poisson_DnO_P2_t(geom, nu, f, gd, gn,t)

% Quad. nodes and weights.
[zita, csi, eta, omega] = int_nodes_weights(5, 'triangle');
nq = length(omega);
N = geom.Nobj.N_node + geom.Nobj.N_edge - size(geom.piv.Di,1);
A = zeros(N, N);
AD = zeros(N, size(geom.piv.Di,1));
b = zeros(N, 1);
uD = zeros(size(geom.piv.Di,1),1);
for e = 1 : geom.Nobj.N_ele

    area = geom.support.TInfo(e).Area;

    v1 = geom.obj.T(e,1);
    v2 = geom.obj.T(e,2);
    v3 = geom.obj.T(e,3);

    x1 = geom.obj.P(v1,1); y1 = geom.obj.P(v1,2);
    x2 = geom.obj.P(v2,1); y2 = geom.obj.P(v2,2);
    x3 = geom.obj.P(v3,1); y3 = geom.obj.P(v3,2);

    x_quad = zita * x3 + csi * x1 + eta * x2;
    y_quad = zita * y3 + csi * y1 + eta * y2;

    delta_x = zeros(3,1);
    delta_y = zeros(3,1);
    delta_x(1) = x3 - x2;
    delta_x(2) = x1 - x3;
    delta_x(3) = x2 - x1;
    delta_y(1) = y2 - y3;
    delta_y(2) = y3 - y1;
    delta_y(3) = y1 - y2;

    B_1 = 1/(2*area) * [delta_y(1) delta_x(1); delta_y(2) delta_x(2)]; 

    A_local = zeros(6, 6);
    b_local = zeros(6,1);
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
            b_local(l) = b_local(l) + omega(q) * phi_local(l) * f(x_quad(q),y_quad(q),t) * (2 * area);
            for k = 1:6
                A_local(l, k) = A_local(l, k) + ...
                    omega(q) * grad_phi_local(k,:) * B_1 * B_1' * grad_phi_local(l,:)' * (2 * area);
            end
        end
    end

    for l = 1 : 6

        jj = geom.piv.piv(geom.obj.T(e,l));

        if (jj > 0)
            for k = 1 : 6
                kk = geom.piv.piv(geom.obj.T(e,k));
                if (kk > 0)
                    A(jj,kk) = A(jj,kk) + nu * A_local(l,k);
                else
                    AD(jj, -kk) = AD(jj, -kk) +  nu * A_local(l,k);
                end

            end

            b(jj) = b(jj) + b_local(l);

        end 

    end

end

for i = 1:size(geom.piv.Di,1)
    n = geom.piv.Di(i,1);
    uD(i) = gd(geom.obj.P(n,1),geom.obj.P(n,2),t);
end

for e = 1:size(geom.piv.Ne, 1)
    id_edge = geom.piv.Ne(e, 1);
    marker = geom.piv.Ne(e, 2);

    v1 = geom.obj.E(id_edge, 1);
    v2 = geom.obj.E(id_edge, 2);
    vm = geom.obj.E(id_edge, end);

    dof_v1 = geom.piv.piv(v1);
    dof_v2 = geom.piv.piv(v2);
    dof_vm = geom.piv.piv(vm);
    coord_v1 = geom.obj.P(v1, :);
    coord_v2 = geom.obj.P(v2, :);
    coord_vm = geom.obj.P(vm, :);

    L = norm(coord_v2 - coord_v1, 2);

    n = zeros(2, 1);
    switch marker
        case 2, n = [0; -1];
        case 4, n = [1; 0];
        case 6, n = [0; 1];
        case 8, n = [-1; 0];
    end

    xi_q = [0.1127016654, 0.5, 0.8872983346]; 
    omega_lato = [5/18, 8/18, 5/18];

    contributo = zeros(3, 1);
    for k = 1:3

        xq = (1 - xi_q(k)) * coord_v1 + xi_q(k) * coord_v2;

        phi_local = [
            2 * (1 - xi_q(k)) * ((1 - xi_q(k)) - 0.5);
            4 * xi_q(k) * (1 - xi_q(k));
            2 * xi_q(k) * (xi_q(k) - 0.5)
        ];

        for j_index = 1:3
            contributo(j_index) = contributo(j_index) + omega_lato(k) * gn(xq(1), xq(2),t)' * n * phi_local(j_index);
        end
    end

    contributo = L * contributo;

    if dof_v1 > 0, b(dof_v1) = b(dof_v1) + contributo(1); end
    if dof_vm > 0, b(dof_vm) = b(dof_vm) + contributo(2); end
    if dof_v2 > 0, b(dof_v2) = b(dof_v2) + contributo(3); end
end

u0 = gmres(A, b - AD*uD,[],1e-12,size(A,1));

u_star = zeros(N,1);

for i = 1:N
    ii = geom.piv.piv(i);
    if ii > 0
        u_star(i) = u0(ii);
    else
        u_star(i) = uD(-ii);
    end
end

end