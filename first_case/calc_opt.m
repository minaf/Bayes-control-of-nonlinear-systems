function [y_new] = calc_opt(reference, noise, controller, b, lx, G)
    y_e = 0;
    e = [];
    uh = 0;
    u_l = 0;
    y_all=[];
    z = tf('z');
    G1 = G*z;

    N = size(reference, 1);
    
    for k=2:N
        e = [reference(k-1) - y_e; e];
        u_l = calculate_u(controller, b, e, u_l, lx);
        uh = [uh; u_l];
        y =  lsim(G1, uh);
        y_e = y(end, 1)+noise(k);
        y_all = [y_all; y_e];
    end
    y_new = lsim(G1, uh)+noise;
   

end

