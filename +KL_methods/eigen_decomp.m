function [u, e] = eigen_decomp(fc, N_KL, w)

    
    [N, N_t] = size(fc);
    Klm = zeros(N_t);

    if nargin < 3
        w = ones(1,N)/(N-1);
    end
    
                  %optimise
      for l = 1:N_t
          for m = 1:N_t
              for k = 1:N
                  Klm(l, m) = Klm(l, m) + w(k)*fc(k,l)*fc(k,m);
%                   if l ~= m
%                     Klm(m, l) = Klm(m, l) + w(k)*fc(k,l)*fc(k,m);
%                   end
              end
          end
      end
    
    
%     Klm = Klm/(N-1); % R^N_quad x N_quad
    
    %% solve the eigenpair problem (W^1/2 K W^1/2 ui = lambda ui)
    %W = diag(w);
    
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