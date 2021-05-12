

% Created on mrt 19 16:03
% Author: Melisa
% Version for CN43

data_dir = '/scratch/melisa/pFC_Criticality/data/activity/';
MODEL_dir = '/scratch/melisa/pFC_Criticality/data/model/';
HC_dir = '/scratch/melisa/pFC_Criticality/data/size_effect/';
figures_dir = '/scratch/melisa/pFC_Criticality/figures/size_effect/';

mice = ['32363'; '32364'; '32365'; '56165'; '56166'];
sessions = {['1','2'];['1','2'];'2';['1','2','4'];['1','3']};
trials = ['1','6','11','16'];

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

            model_name = strcat('model_mouse_',mouse,'_session_',session,'_day_',int2str(day),...
            '_v1.4.20.3.0.1.1.0','.mat');
            
            model_file_dir = strcat(MODEL_dir,model_name);

            heat_capacity_name = strcat('heat_capacity_mouse_',mouse,'_session_',session,...
                '_day_',int2str(day),'_v1.4.20.3.0.1.1.0','.mat');
            
            heat_capacity_dir = strcat(HC_dir, heat_capacity_name);

            load(strcat(data_dir,file_name));

            number_of_cells_flag = 0;
            init_time = 1
            for trial = 1:5 % iteration over trials
                data = calcium_binary.rest_binary{trial};
                data_size = size(data);
                all_data(:,init_time:init_time+data_size(2)-1) = data;
                init_time = init_time+data_size(2);
            end
            data_mean  = mean(all_data,2);
            firing_rate = data_mean *10;    
            index = find(firing_rate >1 & firing_rate<7);

            data_limit = size(index,1);
                
            if data_limit > 100               
                number_of_cells_flag = 1;
                x = all_data(index,:);
                x_size = size(x);
                nneuros = x_size(1);

                counter = 1
                for i=[20,40,60,80]
                    
                    for k =1:5
                        selected_cells = randsample(data_limit,i);
                        model = maxent.createModel(i,'ising');
                        model = maxent.trainModel(model,x(selected_cells,:),'threshold',3);
                        model_ising{trial,counter} = model;

                        C_N = MCMH(i,model, 500000);
                        heat_capacity{counter,k} = C_N;
                    end
                    counter = counter +1 
                end
            end
            
            if number_of_cells_flag == 1
                save(model_file_dir,'model_ising')
                save(heat_capacity_dir,'heat_capacity')
            end
        end
    end