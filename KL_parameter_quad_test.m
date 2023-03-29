clear;
import PCE_methods.*; import KL_methods.*; import general.*;

dt = 0.1; T = 10; time_grid = 0:dt:T; tol = 1e-8;
N_p = 3; ord = 5; N_t = length(time_grid);

%% generate realisations
input_ranges = [0.25 1.25 0.5];
input_means = [0.5 25/8 -1];


[U, w] = spquad(N_p, ord);
I = general.generate_model_inputs(input_ranges, input_means, U);
[u_sp, fc, e] = get_ev(I, w/(2^N_p), time_grid, 8);
fi = KL_methods.generate_KL_modes(fc, e);
[coefficients, basis_index, ~] = KL_methods.PCE_KL_modes(fi, U, 4);
G_tot_sp = KL_methods.calculate_sensitivity_indices(basis_index, coefficients, 'total')
th_den = sum(u_sp);


[U, w] = mv_lgwt(N_p, ord);
I = general.generate_model_inputs(input_ranges, input_means, U);
u_lg = get_ev(I, w, time_grid, 8);

N = 125; w = ones(1,N)/(N-1); U = general.generate_legendre_samples(N_p, N); 
I = general.generate_model_inputs(input_ranges, input_means, U);
u_mc1 = get_ev(I, w, time_grid, 8);

N = 441; w = ones(1,N)/(N-1); U = general.generate_legendre_samples(N_p, N); 
I = general.generate_model_inputs(input_ranges, input_means, U);
[u_mc2, fc, e] = get_ev(I, w, time_grid, 8);
fi = KL_methods.generate_KL_modes(fc, e);
[coefficients, basis_index, ~] = KL_methods.PCE_KL_modes(fi, U, 4);
G_tot_441 = KL_methods.calculate_sensitivity_indices(basis_index, coefficients, 'total')
th_den = sum(u_mc2);

N = 1000; w = ones(1,N)/(N-1); U = general.generate_legendre_samples(N_p, N); 
I = general.generate_model_inputs(input_ranges, input_means, U);
u_mc3 = get_ev(I, w, time_grid, 8);

N = 10000; w = ones(1,N)/(N-1); U = general.generate_legendre_samples(N_p, N); 
I = general.generate_model_inputs(input_ranges, input_means, U);
[u_mc4, fc, e] = get_ev(I, w, time_grid, 8);
fi = KL_methods.generate_KL_modes(fc, e);
[coefficients, basis_index, ~] = KL_methods.PCE_KL_modes(fi, U, 4);
G_tot_10k = KL_methods.calculate_sensitivity_indices(basis_index, coefficients, 'total')
th_den = sum(u_mc4);



figure();
plot(cumsum(u_mc1)/sum(u_mc1)); hold on;
plot(cumsum(u_mc2)/sum(u_mc2));
plot(cumsum(u_mc3)/sum(u_mc3));
plot(cumsum(u_mc4)/sum(u_mc4));
plot(cumsum(u_sp)/sum(u_sp)); 
plot(cumsum(u_lg)/sum(u_lg)); legend('MC 125','MC 441', 'MC 1000', 'MC 10000', 'SP', 'LG'); hold off;


figure();
semilogy(u_mc1/u_mc1(1)); hold on;
semilogy(u_mc2/u_mc2(1));
semilogy(u_mc3/u_mc3(1));
semilogy(u_mc4/u_mc4(1));
semilogy(u_sp/u_sp(1));
semilogy(u_lg/u_lg(1)); legend('MC 125','MC 441', 'MC 1000', 'MC 10000', 'SP', 'LG'); hold off;





function [u, fc, e] = get_ev(I, w, time_grid, N_KL)

    %% generate N realisations of f, and then centre f
    fmk = f_mechosc(I, time_grid); [N, N_t] = size(fmk);
    mean_process = sum(fmk)/N; fc = fmk - mean_process;     
    
    Klm = zeros(N_t);
    
    for k = 1:N               %optimise
        for l = 1:N_t
            for m = 1:N_t
                Klm(l, m) = Klm(l, m) + w(k)*fc(k,l)*fc(k,m);
            end
        end
    end
    
%     Klm = Klm/(N-1); % R^N_quad x N_quad
%      Klm = Klm/0.1;

    [e, U] = eig(Klm);  % [vectors, values] (for w = 1)
    [u, ~] = sort(diag(U), 'descend'); %[u, vec_order]
%     ind = 1:length(u);
    
    %N_KL = N_quad; % non truncated form
%     norm_u = u/max(u);
%     ind = ind(norm_u > tol);
%     vec_order = vec_order(norm_u > tol);
    
    u = u(1:N_KL); % u(ind);
    e = e(:, 1:N_KL); % e(:,vec_order);
    N_KL = length(u); % truncated form
end