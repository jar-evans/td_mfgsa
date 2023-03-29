clear;
import PCE_methods.*; import KL_methods.*; import general.*;

%% set time grid
T = 10;
dt = .1;
N_quad = int64(T/dt) + 1; %length(0:dt:T);

%% turn plots on/off
u_PLOTS = true; 
e_PLOTS = true;
v_PLOTS = true;
approx_PLOTS = false;

%% generate realisations
input_ranges = [0.25 1.25 0.5];
input_means = [0.5 25/8 -1];

N_p = 3; N = 100;

% [U, N] = general.generate_legendre_samples(grid_h, N_p);
U = general.generate_legendre_samples(N_p, N); 
I = general.generate_model_inputs(input_ranges, input_means, U);

%% generate N realisations of f, and then centre f
fmk = f_mechosc(I, 0:dt:T);
mean_process = sum(fmk)/N;
fc = fmk - mean_process;

%% form the covariance matrix (discretised covariance function) and compute dicretised KL modes
N_KL = 5;
[u, e] = KL_methods.eigen_decomp(fc, N_KL);

KL_methods.plot_variance(u, e, fc, v_PLOTS);
KL_methods.plot_eigpairs(u, e, u_PLOTS, e_PLOTS);

fi = KL_methods.generate_KL_modes(fc, e);

approx = KL_methods.generate_approximation(fi, e, approx_PLOTS, mean_process);

N_ord = 4;
[coefficients, basis_index, ~] = KL_methods.PCE_KL_modes(fi, U, N_ord);

G_main = KL_methods.calculate_sensitivity_indices(basis_index, coefficients, 'main')
G_tot = KL_methods.calculate_sensitivity_indices(basis_index, coefficients, 'total')




















