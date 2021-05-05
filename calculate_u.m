function [y_d] = calculate_u(controller, b, e, u_l, lx)
    
tmp = size(e, 1);
if tmp<lx
    e = [e; zeros(lx-tmp, 1)];
end
y_d = controller'*e(1:lx, 1)-b*u_l;

end

