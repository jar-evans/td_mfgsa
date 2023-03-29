classdef KLPC_SURROGATE < TD_SURROGATE
    %PC_SURROGATE Summary of this class goes here
    %   Detailed explanation goes here

    properties
        eigenvectors
        eigenvalues
        modes
        coefficients
        basis_index
    end

    methods
        function obj = KLPC_SURROGATE(params, N_KL, N_ord, w)
            %PC_SURROGATE Construct an instance of this class
            %   Detailed explanation goes here

            obj@TD_SURROGATE(params);  
            [obj.eigenvalues, obj.eigenvectors] = KL_methods.eigen_decomp(params.fc, N_KL, w);
            obj.modes = KL_methods.generate_KL_modes(params.fc, obj.eigenvectors);
            [obj.coefficients, obj.basis_index, ~] = KL_methods.PCE_KL_modes(obj.modes, params.U, N_ord);
        

        end

        function SI = calculate_sensitivity_indices(obj, type)
            SI = KL_methods.calculate_sensitivity_indices(obj.basis_index, obj.coefficients, type);
        end

        function plot_eigpairs(obj)
            KL_methods.plot_eigpairs(obj.eigenvalues, obj.eigenvectors, true, true)
        end


        function plot_variance(obj)
            KL_methods.plot_variance(obj.eigenvalues, obj.eigenvectors, obj.fc, true)
        end
    end
end
