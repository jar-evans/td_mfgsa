function U = generate_legendre_samples(N_p, N, grid_h)

    if exist('grid_h', 'var')

        grids = cell(1, N_p);    
        
        [grids{1:N_p}] = ndgrid(-1:grid_h:1);

        N = numel(cell2mat(grids(1)));
        U = zeros(N, N_p);

        for grid = 1:N_p
            U(1:N, grid) = reshape(cell2mat(grids(grid)), N, 1);
        end

    else
        U = 2*rand(N, N_p) - 1;
    end

    
end

% [X1,X2] = ndgrid(-1:1:1);
% N = numel(X1);
% U_a = reshape(X1, N, 1);
% U_b = reshape(X2, N, 1);
% 
% Utest = horzcat(U_a, U_b)
% 
% [X1,X2,X3] = ndgrid(-1:0.12:1);
% N = numel(X1);
% U_a = reshape(X1, N, 1);
% U_b = reshape(X2, N, 1);
% U_l = reshape(X3, N, 1);
% 
% % U_a = 2*rand(N, 1) - 1;
% % U_b = 2*rand(N, 1) - 1;
% % U_l = 2*rand(N, 1) - 1;
% 
% U = horzcat(U_a, U_b, U_l);