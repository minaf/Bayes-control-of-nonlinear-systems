N = 500;
u = randn(N, 1);
ts = -1;
n_a = 5;

%c_true = [0 2.8 -4.5 5.7];
F = [];
v = zeros(N, 1);
for k = 1:N
    if u(k)>1
        v(k) = u(k)-1;
    end
    
    if u(k)<-1
        v(k) = u(k)+1;
    end
    
end

%v = u./(1+u.^2);
%v = F*c_true';

G = tf([0, 0, 0, 0.28261, 0.50666],[1,-1.41833,1.58939,-1.31608,0.88642], ts);
%G = tf([0, 1], [1 -0.5], -1);
y0 = lsim(G, v);


controller_ideal = [0.2045; -0.2715; 0.2931; -0.2396; 0.1643; 0.0084];
b = -1;
lx = size(controller_ideal, 1);
K_i = tf(controller_ideal', [1 b zeros(1, lx-2)], -1);

Mz = minreal(G*K_i/(1+G*K_i));

b = -1;
lx = 6;


mc = 50;
tolerance = 1e-3;
n = 100;
m1 = 10;
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
F1 = [];
for k = 1:N
    f = [];
    for l = 1:n_a+1
        f = [f u(k).^(l-1)];
    end
    F1 = [F1; f];
end

sigma = 0.5;
sigma_all = zeros(2, mc);
%% ideal controller
% 
% fu = @(x) fun_ideal(x, G, Mz, ts, b, lx);
% [x, fval]= fmincon(fu, zeros(lx, 1), [], []);
% controller_ideal = x;

%%

%sigma = 0.5;
 MC = 10;
 M = 200;
 cost_all = zeros(mc, 5);
%%
for kl = 1:mc
%%
kl
o = sigma*randn(N, 1);
y = y0+o;

% [P1, m1, theta_t] = hammerstein_pw(F, u, y, n, tolerance);
% 
%%
[P, m, theta, sig] = hammerstein_kernel_ident(F, u, y, n, tolerance);
sigma_all(kl, 1) = sig;
%%
controller = bc_new_fmin(P, m, Mz, lx, b, ts, M, MC);


%%
% [P1, m1, theta_t] = hammerstein_pw(F, u, y, n, tolerance);
% 
% %%
% tic
% controller_t = bc_new_fmin(P1, m1, Mz, lx, b, ts, M, MC);
% toc


%%

[P1, m1, theta1, sig] = hammerstein_kernel_ident(F1, u, y, n, tolerance);
sigma_all(kl, 2) = sig;
%[P1, m1, theta1_t] = hammerstein_kernel_ident(F1, u, y, n, tolerance);
% 
%%
tic
controller1 = bc_new_fmin(P1, m1, Mz, lx, b, ts, M, MC);
toc

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


%%
%K = tf(controller', [1 b zeros(1, lx-2)], -1);
r = [0; randn(M-1, 1)];
%r = [0; ones(20, 1); 2*ones(20, 1);ones(20, 1);3*ones(20, 1); zeros(19, 1)];
yd = lsim(Mz, r);

o1 = sigma*randn(M, 1);

y_new = calc_pw(r, o1, controller, b, lx, G, u_grid, theta);
y_c = calc_pol(r, o1, controller1, b, lx, G, theta1);
y_opt = calc_opt(r, o1, controller_ideal, b, lx, G);
y_lin = calc_linear(r, o1, controller2, b, lx, G);

cost_all(kl, 1) = norm(y_new-yd, 2)^2/N; %cost of piecewise linear app
cost_all(kl, 2) = norm(y_c-yd, 2)^2/N; %cost of polynomial app
cost_all(kl, 3) = norm(y_opt-yd, 2)^2/N; %cost of optimal
cost_all(kl, 4) = norm(yd-y_lin, 2)^2/N; %cost of linear model (ignoring nonlinearity)
end

