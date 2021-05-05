function [y_new] = calc_linear(reference, noise, controller, b, lx, G)
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
        p_n = u_l;
        if p_n>1
            p_n= p_n-1;
        elseif p_n<-1
            p_n = p_n+1;
        else
            p_n=0;
        end
        
        uh = [uh; p_n];
        y =  lsim(G1, uh);
        y_e = y(end, 1)+noise(k);
        y_all = [y_all; y_e];
    end


    y_new = lsim(G1, uh)+noise;
   

end

