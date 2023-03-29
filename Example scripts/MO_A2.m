clear;

PB_CHOICE = 2; %1: Mech. oscillator, 2: SIR
ALG_CHOICE = 3; %1: MC, 2: PCE, 3: KL-PCE


switch PB_CHOICE
    case 1

        %% set basic params for surrogate
        time_grid = 0:0.1:10;
        input_range = [0.25 1.25 0.5];
        input_mean = [0.5 25/8 -1];
        f = @(I) f_mechosc(I, time_grid);
        
        N = 5000; N_ord = 4;
        
        %% establish surrogate
        base = TD_SURROGATE(f, time_grid, input_range, input_mean);

    case 2

        %% set basic params for surrogate
        time_grid = 0:0.1:10;
        input_range = [0.25 1.25 0.5];
        input_mean = [0.5 25/8 -1];
        f = @(I) f_mechosc(I, time_grid);
        
        N = 5000; N_ord = 4;
        
        %% establish surrogate
        base = TD_SURROGATE(f, time_grid, input_range, input_mean);
end

switch ALG_CHOICE
    case 1
    case 2
        pc_surrogate = base.generate_PC_surrogate(N, N_ord);
        
        %% generate new realisations of the process from the surrogate
        % (not necessary just fun to see)

%         [U_new, I_new] = pc_surrogate.sample_inputs(100);
%         w = pc_surrogate.generate_realisations(U_new); 
% 
%         figure();
%         plot(pc_surrogate.time_grid, w); hold on;
%         plot(pc_surrogate.time_grid, pc_surrogate.mean_process, 'PC'); hold off;
        

        %% calculate sensitivity indices
        G_tot_td = pc_surrogate.calculate_sensitivity_indices(true, 'total');
        G_tot_pw = pc_surrogate.calculate_sensitivity_indices(false, 'total');
        
        G_main_td = pc_surrogate.calculate_sensitivity_indices(true, 'main');
        G_main_pw = pc_surrogate.calculate_sensitivity_indices(false, 'main');
        
        figure(); clf
        plot(linspace(0,10,length(G_tot_pw)), G_tot_pw);
        bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex'; xlabel('Time (s)');
        title('Pointwise total Si mech. oscillator')

        figure(); clf
        plot(linspace(0,10,length(G_tot_td)), G_tot_td);
        bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex'; xlabel('Time (s)');
        title('Generalised total Si mech. oscillator')

%         figure(2); clf
% h = boxplot([pc_st(:,end,1,1), pc_st(:,end,1,2), pc_st(:,end,1,3),...
%              pc_st(:,end,2,1), pc_st(:,end,2,2), pc_st(:,end,2,3), ...
%              pc_st(:,end,3,1), pc_st(:,end,3,2), pc_st(:,end,3,3)],...
%     'Colors',[blue; red; something; blue; red; something; blue; red; something],'Whisker',10,...
%     'labels',{'$s_t^1$','$s_t^2$','$s_t^3$', '$s_t^1$','$s_t^2$','$s_t^3$', '$s_t^1$','$s_t^2$','$s_t^3$'}); hold on;
% set(h,{'linew'},{2}); grid on
% legend(flipud(findall(gca,'Tag','Box')), {'p = 100', 'p = 1000', 'p = 10000'},...
%     'Location','NorthEast','interpreter','latex'); legend boxoff
% bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
% title('Total SI spread for mech. oscillator - 1FHF - PCE','interpreter','latex')
        
        figure(); clf
        plot(linspace(0,10,length(G_tot_pw)), G_tot_pw); hold on;
        plot(linspace(0,10,length(G_tot_td)), G_tot_td); hold on; hold off;
        bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';  xlabel('Time (s)');
        title('Gen. and pointwise total Si mech. oscillator')

        figure(); clf
        plot(linspace(0,10,length(G_main_pw)), G_main_pw); hold on;
        plot(linspace(0,10,length(G_main_td)), G_main_td); hold on; hold off;
        bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';  xlabel('Time (s)');
        title('Gen. and pointwise main Si mech. oscillator')
        
        figure();
        semilogy(pc_surrogate.time_grid, pc_surrogate.pointwise_variance())
    case 3

        klpc_surrogate = base.generate_KLPC_surrogate(N, 1e-6, 4);
% 
        klpc_surrogate.plot_eigpairs();
        klpc_surrogate.plot_variance();

        %% calculate sensitivity indices
        G_tot_td = klpc_surrogate.calculate_sensitivity_indices('total')
%         klpc_surrogate.eigenvalues

%         
%         G_main_td = klpc_surrogate.calculate_sensitivity_indices('main');

%         Ts = 0.1:0.1:10;
%         tot = zeros(length(Ts), 3);
% 
%         i = 1;
%         for T = Ts
%             time_grid = 0:0.01:T;
%             f = @(I) f_mechosc(I, time_grid);
%             base = TD_SURROGATE(f, time_grid, input_range, input_mean);
%             klpc_surrogate = base.generate_KLPC_surrogate(N, 1e-1, 4);
% 
%             tot(i, :) = klpc_surrogate.calculate_sensitivity_indices('total');
%             i = i + 1;
%         end
% 
% 
%         figure();
%         plot(Ts, tot);


end