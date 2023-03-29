clear;
import PCE_methods.*; import general.*;

%% set time grid
T = 10;
dt = .01;
N_quad = int64(T/dt) + 1; %length(0:dt:T);

%% generate realisations
input_ranges = [0.25 1.25 0.5];
input_means = [0.5 25/8 -1];

N_p = 3; grid_h = 0; N = 1000;

% [U, N] = general.generate_legendre_samples(grid_h, N_p);
U = general.generate_legendre_samples(N_p, N); 
I = general.generate_model_inputs(input_ranges, input_means, U);

%% generate N realisations of f, and then centre f
t = 0:dt:T;
fmk = f_mechosc(I, t);
mean_process = sum(fmk)/size(fmk,1);
fc = fmk - mean_process; 

%% create projection matrix and calculate the PCE coeffs
N_ord = 4;
[coefficients, basis_index, proj_matrix] = PCE_methods.PCE(fc, U, N_ord);
% [basis_index, N_PC] = PCE_methods.generate_basis_index(N_ord, N_p);
% proj_matrix = PCE_methods.generate_projection_matrix(basis_index, U, N_p, N_PC, N_PC_quad);
% coefficients = PCE_methods.calculate_PCE_coefficients(proj_matrix, fc);


%% check approximations
y = PCE_methods.generate_approximation(proj_matrix, coefficients);
yp = PCE_methods.generate_approximation(proj_matrix, coefficients, sum(basis_index, 2) == 0);

figure(1);
foo = f_mechosc(input_means, t);
plot(t, y + mean_process); hold on;
plot(t, foo); plot(t, yp + mean_process); plot(t, min(fmk));
plot(t, max(fmk)); hold off;

%% calculate sensitivity indices
G_tot_td = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, true, 'total', dt);
G_tot_pw = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, false, 'total', dt);

G_main_td = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, true, 'main', dt);
G_main_pw = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, false, 'main', dt);

figure(2);
plot(linspace(0,T, length(G_tot_pw)), G_tot_pw); hold on;
plot(linspace(0,T, length(G_tot_td)), G_tot_td); hold off;

figure(3);
plot(linspace(0,T, length(G_main_pw)), G_main_pw); hold on;
plot(linspace(0,T, length(G_main_td)), G_main_td); hold off;