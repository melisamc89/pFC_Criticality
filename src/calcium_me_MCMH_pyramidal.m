

% Created on mrt 19 16:03
% Author: Melisa
% Version for CN43

data_dir = '/scratch/melisa/pFC_Criticality/data/activity/';
MODEL_dir = '/scratch/melisa/pFC_Criticality/data/model/';
HC_dir = '/scratch/melisa/pFC_Criticality/data/heat_capacity/';
figures_dir = '/scratch/melisa/pFC_Criticality/figures/';

mice = ['32363'; '32364'; '32365'; '56165'; '56166'];
sessions = {['1','2'];['1','2'];'2';['1','2','4'];['1','3']};
trials = ['1','6','11','16'];

ncells = 25;
random_samples = 30;
n_group_cells = 10;
combos = combntns(1:random_samples,ncells);
size_combos = size(combos);
sample = randsample(size_combos(1),random_samples);

prompt = 'mouse id : '
mouse_it = input(prompt)

mouse = mice(mouse_it,:)
nsessions = size(sessions{mouse_it});
    for session_it = 1:nsessions(2)
        
        session = sessions{mouse_it}(session_it)
        
        for day = [1,2,3,4]
            
            trial_id = trials(day);

            file_name = strcat('mouse_',mouse,'_session_',session,'_trial_',trial_id,...
            '_v1.4.20.3.0.1.1.0.mat');

            model_name = strcat('model_pyramidal_mouse_',mouse,'_session_',session,'_trial_',trial_id,...
            '_v1.4.20.3.0.1.1.0_ncells_',int2str(ncells),'.mat');
            
            model_file_dir = strcat(MODEL_dir,model_name);

            heat_capacity_name = strcat('hc_pyramidal_mouse_',mouse,'_session_',session,...
                '_trial_',trial_id,'_v1.4.20.3.0.1.1.0_ncells_',int2str(ncells),'.mat');
            
            heat_capacity_dir = strcat(HC_dir, heat_capacity_name);

            load(strcat(data_dir,file_name));

            number_of_cells_flag = 0;
            for trial = 1:5 % iteration over trials
    
                data = calcium_binary.rest_trace{trial};
                data_mean  = mean(data,2);
                firing_rate = data_mean * 10;
                index = find(firing_rate >1 & firing_rate<12);

                %[sorted_data index] = sort(data_mean, 'descend');
                
                data_binary = calcium_binary.rest_binary{trial};
                %data_limit = size(data,1)*50/100;
                data_limit = size(index,1)
                
                if data_limit > random_samples
                    number_of_cells_flag = 1;
                    x = data_binary(index,:);
                    x_size = size(x);
                    nneuros = x_size(1);

                    for i=1:n_group_cells

                        model = maxent.createModel(ncells,'ising');
                        model = maxent.trainModel(model,x(combos(sample(i),:),:),'threshold',1);
                        model_ising{trial,i} = model;

                        C_N = MCMH(ncells,model, 100000);
                        heat_capacity{trial,i} = C_N;

                    end
                end
            end
            
            if number_of_cells_flag == 1
                save(model_file_dir,'model_ising')
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
                title(strcat('Heat Capacity Pyramidal cells. Mouse:', mouse, ' session:', session, ' day:', int2str(day)))
                figure_name = strcat(figures_dir, 'hc_pyramidal_mouse_',mouse,'_session_',session,'_trial_',trial_id,'_v1.4.20.3.0.1.1.0_ncells_',int2str(ncells),'.png');
                saveas(1,figure_name)
                close(1)
            end
        end
    end