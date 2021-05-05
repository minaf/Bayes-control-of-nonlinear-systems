function [y_new] = calc_pw(reference, noise, controller, b, lx, G, u_grid, theta)
    y_e = 0;
    e = [];
    uh = 0;
    u_l = 0;
    y_all=[];
    z = tf('z');
    G1 = G*z;
    m = size(u_grid, 2);
    x_grid = zeros(1, m);
    for k=1:m
        x_grid(1, k) = calculate_output_pw(u_grid(1,k), u_grid, theta);
    end

    N = size(reference, 1);
    
    for k=2:N
        e = [reference(k-1)-y_e; e];
        u_l = calculate_u(controller, b, e, u_l, lx);
        p_n = cancel_nonlinearity_pw(u_l, u_grid, x_grid, theta);
        
        if p_n>1
            p_n= p_n-1;
        elseif p_n<-1
            p_n = p_n+1;
        else
            p_n=0;
        end
        
        uh = [uh; p_n];
        y =  lsim(G1, uh);
        y_e = y(end, 1)+noise(k,1);
        y_all = [y_all; y_e];
    end


    y_new = lsim(G1, uh)+noise;
   

end

