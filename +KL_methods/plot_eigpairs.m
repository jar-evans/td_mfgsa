function plot_eigpairs(u, e, u_PLOTS, e_PLOTS)
    if (u_PLOTS)
        figure();
        plot(cumsum(u)/sum(u));
        title('Truncation level quality')
        
        figure();
        norm_u = u/u(1);
        semilogy(norm_u);
        title('Normed eigenvalues')
    end
    
    if (e_PLOTS)
        figure();
        plot(e); 
        title('Eigenvectors')
    end

end

