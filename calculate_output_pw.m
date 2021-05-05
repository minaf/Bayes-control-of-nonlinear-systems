function [output] = calculate_output_pw(input,u_grid, theta)
m = size(theta, 1);
f1 = 1;
    for l=2:m
        tmp = 0;
        if input<=u_grid(l) && input>= u_grid(l-1)
            tmp = input-u_grid(l-1);
        end
  
        if input>u_grid(l)
            tmp = u_grid(l)-u_grid(l-1);
        end
        f1 = [f1 tmp];
    end
    output = f1*theta;
end

