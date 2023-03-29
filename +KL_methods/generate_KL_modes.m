function fi = generate_KL_modes(fc, e)

    N = size(fc, 1);
    N_KL = size(e, 2);
    fi = zeros(N_KL, N);

    %w = %;%/0.1;
    
    for i = 1:N_KL
        fi(i, :) = fc*e(:,i); %w*(fc*e(:,i));
    end

end

