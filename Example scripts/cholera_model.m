clear;

T = 250;
dt = .05;%e-2;
N_quad = T/dt;

 % uncertain params
bH = 7.5;
kL = 1e6;
X = 168/5;
Z = 70;
g = 7/5;

%% generate realisations

input_means = [bH, kL, X, Z, g];
input_ranges = 0.2.*input_means;

N_p = length(input_means); N = 100;
% 
% [U, N] = general.generate_legendre_samples(grid_h, N_p);
U = general.generate_legendre_samples(N_p, N); 
I = general.generate_inputs(input_ranges, input_means, U);

%% generate N realisations of f, and then centre f
[y, t] = f_cholera(I, 0:dt:T, 2, 45);
% [t2, y2] = f_cholera(I, 0:dt:100, 2, 45);


figure(1);
loglog(t, y);
ylim([5e-2 1e4]);
xlim([5e-2 2e2]);

% figure(1);
% loglog(t, y(:, 1:3));
% ylim([5e-2 1e4]);
% xlim([5e-2 2e2]);
% 
% figure(2);
% loglog(t, y(:, 4:5));
% xlim([5e-2 2e2]);

% %% set time grid % time discr may not always be unique
% N_t = length(time_grid);    % #points in time discretisation
% N_h = N_t - 1;               % #intervals in time discretisation
% 
G_alg2 = gen_alg2(y, t, U, 4, true, 'total');

semilogx(G_alg2')
% G_alg2 = gen_alg2(y2, t2, U, 4, true, 'main')
% 
% G_alg3 = gen_alg3(fmk, time_grid, U, 1e-5, 4, 'main');

% G_alg3 = gen_alg3(y, t, U, 1e-3, 3, 'total')
% figure(2);
% semilogx(0:dt:T, G_alg2);

