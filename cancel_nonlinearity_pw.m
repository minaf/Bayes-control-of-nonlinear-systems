function [output] = cancel_nonlinearity_pw(inp, input_grid, output_grid, theta)
inp
output_grid
m = size(input_grid, 2);

if inp<min(output_grid(1, 1), output_grid(1, end) )
    inp = min(output_grid(1, 1), output_grid(1, end) );
end

if inp>max(output_grid(1, 1), output_grid(1, end) )
    inp = max(output_grid(1, 1), output_grid(1, end) );
end


tmp = 1;
for k = 2:m
    if inp<=max(output_grid(1, k), output_grid(1, k-1)) && inp>=min(output_grid(1, k), output_grid(1, k-1))
        tmp = k;
        break;
    end
end


f1 =1;
for l=2:tmp-1
    f1 = [f1 input_grid(l)-input_grid(l-1)];
end

n = f1*theta(1:tmp-1, 1);

input_grid(1, tmp-1)

output = (-n+inp)/theta(tmp, 1)+input_grid(1, tmp-1);
