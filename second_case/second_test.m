N = 500;
sigma = 10;
u = randn(N, 1);
ts = -1;
c_true = [0 0 0 1];
n_a = size(c_true, 2)-1;

c_true = [0 2.8 -4.5 5.7];
F = [];
for k = 1:N
    f = [];
    for l = 1:n_a+1
        f = [f u(k).^(l-1)];
    end
    F = [F; f];
end

%v = u./(1+u.^2);
v = F*c_true';

G = tf([0, 1, 0.6],[1, -1, 0.8], ts);

y0 = lsim(G, v);

alpha = exp(-0.5);
Mod = [0 0 0 (1-alpha)^2; 1 -2*alpha alpha^2 0];
Mz = tf(Mod(1, :),Mod(2, :),ts);

b = -1;
lx = 3;


mc = 5;
tolerance = 1e-3;
n = 100;

m1 = 10; %number of grid points for a piecewise linear function
u_min = min(u);
u_max = max(u);
u_grid = linspace(u_min, u_max, m1);

F = [];
for k=1:N
    f1 = 1;
    for l=2:m1
        tmp = 0;
        if u(k)<=u_grid(l) && u(k)>= u_grid(l-1)
            tmp = u(k)-u_grid(l-1);
        end
         
        if u(k)>u_grid(l)
            tmp = u_grid(l)-u_grid(l-1);
        end
        f1 = [f1 tmp];
    end
    F = [F; f1];
end

%% ideal controller

fu = @(x) fun_ideal(x, G, Mz, ts, b, lx);
[x, fval]= fmincon(fu, zeros(lx, 1), [], []);
controller_ideal = x;

%%

 MC = 10; %number of different g samples
 M = 200; %length of horizon
 cost_all = zeros(mc, 4);
%%
for kl = 1:mc
%%
kl
o = sigma*randn(N, 1);
y = y0+o;

[P, m, theta] = hammerstein_kernel_ident(F, u, y, n, tolerance);

%%
controller = bc_new_fmin(P, m, Mz, lx, b, ts, M, MC);

cost= calculate_cost(controller, Mz, G, b, ts);

%%
F1 = [];
n_a = 5;
for k = 1:N
    f = [];
    for l = 1:n_a+1
        f = [f u(k).^(l-1)];
    end
    F1 = [F1; f];
end

[P, m, theta1] = hammerstein_kernel_ident(F1, u, y, n, tolerance);

%%
controller1 = bc_new_fmin(P, m, Mz, lx, b, ts, M, MC);
cost1= calculate_cost(controller1, Mz, G, b, ts);
%%
z = [y u];
orders = [0 n 0];
option = arxRegulOptions('RegularizationKernel','TC');
arxOpt = arxOptions;

[Lambda,Reg] = arxRegul(z,orders,option);
arxOpt.Regularization.Lambda = Lambda;
arxOpt.Regularization.R = Reg;
model = arx(z,orders,arxOpt);
med = getpvec(model, option);


Phi = toeplitz(u(n:end), u(n:-1:1));
y_all = y(n:end);
Nn = size(y_all, 1);


lambda_est = (y_all-Phi*med)'*(y_all-Phi*med)/Nn;
K = (eye(n)/(Lambda*Reg))*lambda_est;

R = K - K*Phi'/(Phi*K*Phi'+lambda_est*eye(size(Phi, 1)))*Phi*K;

R = (R'+R)/2;


controller2 = bc_new_fmin(P, med, Mz, lx, b, ts, M, MC);

cost2= calculate_cost(controller2, Mz, G, b, ts);
%%

r = [0; 50*randn(M-1, 1)];

yd = lsim(Mz, r); %reference
o1 = sigma*randn(M, 1); %another noise realization

y_new = calc_pw(r, o1, controller, b, lx, G, c_true', u_grid, theta);
y_c = calc_pol(r, o1, controller1, b, lx, G, theta1, c_true');
y_opt = calc_pol(r, o1, controller_ideal, b, lx, G, c_true', c_true');
y_lin = calc_linear(r, o1, controller2, b, lx, G, c_true');

cost_all(kl, 1) = norm(y_new-yd, 2)^2/M; % cost of piecewise linear app
cost_all(kl, 2) = norm(y_c-yd, 2)^2/M; % cost of polynomial app
cost_all(kl, 3) = norm(y_opt-yd, 2)^2/M; %cost of optimal
cost_all(kl, 4) = norm(yd-y_lin, 2)^2/M; %cost of linear model (ignoring nonlinearity)
end

