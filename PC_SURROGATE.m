classdef PC_SURROGATE < TD_SURROGATE
    %PC_SURROGATE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        coefficients
        basis_index
    end
    
    methods
        function obj = PC_SURROGATE(params, N_ord)
            %PC_SURROGATE Construct an instance of this class
            %   Detailed explanation goes here

            obj@TD_SURROGATE(params);            
            [obj.coefficients, obj.basis_index, ~] = PCE_methods.PCE(params.fc, params.U, N_ord);

        end

        function w = generate_realisations(obj, I)
            w = PCE_methods.generate_realisations(obj.mean_process, obj.coefficients, obj.basis_index, I);
        end
        
        function SI = calculate_sensitivity_indices(obj, generalised, type)
            SI = PCE_methods.calculate_sensitivity_indices(obj.basis_index, obj.coefficients, generalised, type, obj.time_grid);
        end

        function v = pointwise_variance(obj)
            v = sum(obj.coefficients, 1);
            v = v.^2;
        end

    end
end

