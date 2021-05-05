function [y_new] = calc_true_cost_linear(reference, noise, controller, b, lx, G, c_true)
    y_e = 0;
    e = [];
    uh = 0;
    u_l = 0;
    y_all=[];
    z = tf('z');
    G1 = G*z;

    N = size(reference, 1);
    %n_a1 = size(c_est, 1)-1;
    n_a2 = size(c_true, 1)-1;
    for k=2:N
        e = [reference(k-1) - y_e; e];
        u_l = calculate_u(controller, b, e, u_l, lx);
        p_n = u_l;
        f = [];
        for l = 1:n_a2+1
            f = [f p_n.^(l-1)];
        end
        uh = [uh; f*c_true];
        y =  lsim(G1, uh);
        y_e = y(end, 1)+noise(k);
        y_all = [y_all; y_e];
    end


    y_new = lsim(G1, uh)+noise;
   

end

