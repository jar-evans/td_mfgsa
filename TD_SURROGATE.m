classdef TD_SURROGATE
    %TD_SURROGATE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        f
        fmk
        fc
        U
        I
        mean_process
        time_grid{mustBeNumeric}
        input_range{mustBeNumeric}
        input_mean{mustBeNumeric}
    end
    
    methods
        function obj = TD_SURROGATE(f, time_grid, input_range, input_mean)
            %TD_SURROGATE Construct an instance of this class
            %   Detailed explanation goes here

            if nargin == 4
                obj.f = f;
                obj.time_grid = time_grid;
                obj.input_mean = input_mean;
                obj.input_range = input_range;

            elseif nargin == 1
                p = properties(f);
                for i=1:length(p)
                    obj.(p{i}) = f.(p{i});
                end
            end
        end       


        function [U, I, w] = sample_inputs(obj, N)
            N_p = length(obj.input_mean);    

            if ~N
                [U, w] = spquad(N_p, 4);
                warning('sparse quad. in param. space');
            else
                U = general.generate_legendre_samples(N_p, N);
                w = ones(1, N)/(N-1);
            end

            I = general.generate_model_inputs(obj.input_range, obj.input_mean, U);            
        end

        function obj = sample_model(obj, N)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            [obj.U, obj.I] = obj.sample_inputs(N);
            obj.fmk = obj.f(obj.I);
        end

        function obj = centre_process(obj)
            obj.mean_process = sum(obj.fmk)/size(obj.fmk,1);
            obj.fc = obj.fmk - obj.mean_process;
            
        end

%         function generate_MC_SI(obj)
%             PC_obj = PC_SURROGATE(obj);
%         end

        function PC_obj = generate_PC_surrogate(obj, N, N_ord)

            obj = obj.sample_model(N);
            obj = obj.centre_process();

            if ~exist('N_ord', 'var')
                N_ord = 4;
            end

            PC_obj = PC_SURROGATE(obj, N_ord);
        end

        function KLPC_obj = generate_KLPC_surrogate(obj, N, N_KL, N_ord)


            [obj.U, obj.I, w] = obj.sample_inputs(N);
            obj.fmk = obj.f(obj.I);


            obj = obj.centre_process();

            if ~exist('N_ord', 'var')
                N_ord = 4;
            end

            if ~exist('N_KL', 'var')
                N_KL = 8;
            end

            KLPC_obj = KLPC_SURROGATE(obj, N_KL, N_ord, w);
        end
    end
        
end

