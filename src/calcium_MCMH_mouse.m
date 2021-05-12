

% Created on mrt 19 16:03
% Author: Melisa
% Version for CN43

data_dir = '/scratch/melisa/pFC_Criticality/data/activity/';
model_dir = '/scratch/melisa/pFC_Criticality/data/model/';
HC_dir = '/scratch/melisa/pFC_Criticality/data/heat_capacity/';
figures_dir = '/scratch/melisa/pFC_Criticality/figures/';

mice = ['32363'; '32364'; '32365'; '56165'];
sessions = {['1','2'];['1','2'];['2'];['1','2','4']};
trials = ['1','6','11','16'];

ncells = 25;
random_samples = 30;
n_group_cells = 20;
combos = combntns(1:random_samples,ncells);
size_combos = size(combos);
sample = randsample(size_combos(1),random_samples);

mouse = '56165'
session = '1'
        for day = [1,2,3,4]
            
            trial_id = trials(day)

            file_name = strcat('mouse_',mouse,'_session_',session,'_trial_',trial_id,...
            '_v1.4.20.3.0.1.1.0.mat');

            model_name = strcat('model_mouse_',mouse,'_session_',session,'_trial_',trial_id,...
            '_v1.4.20.3.0.1.1.0_ncells_',int2str(ncells),'.mat');
            model_dir = strcat(model_dir,model_name);

            heat_capacity_name = strcat('heat_capacity_mouse_',mouse,'_session_',session,...
                '_trial_',trial_id,'_v1.4.20.3.0.1.1.0_ncells_',int2str(ncells),'.mat');
            heat_capacity_dir = strcat(HC_dir, heat_capacity_name);

            load(strcat(data_dir,file_name));

   
            for trial = 1:5 % iteration over trials
    
                x = calcium_binary.rest_binary{trial};
                x_size = size(x);
                nneuros = x_size(1);
    
                for i=1:n_group_cells

                    model = maxent.createModel(ncells,'ising');
                    model = maxent.trainModel(model,x(combos(sample(i),:),:),'threshold',1);
                    model_ising{trial,i} = model;
                    
                    C_N = MCMH(ncells,model, 20000);
                    heat_capacity{trial,i} = C_N;
                    
                end
            end
            save(model_dir,'model_ising')
            save(heat_capacity_dir,'heat_capacity')
            
            figure(1)
            for trial = 1:5    
                subplot(2,3,trial)
                for i=1:n_group_cells
                    plot(0.5:0.1:3,heat_capacity{trial,i})
                    hold on
                    xlabel('T') 
                    ylabel('C_N') 
                    ylim([1 12])
                end
            end
            subplot(2,3,2)
            title(strcat('Heat Capacity Mouse:', int2str(mouse), ' session:', int2str(session), ' day:', int2str(day)))
            figure_name = strcat(figures_dir, 'heat_capacity_mouse_',mouse,'_session_',session,'_trial_',trial,'_v1.4.20.3.0.1.1.0_ncells_,',int2str(ncells),'.png');
            saveas(1,figure_name)
            
        end
