

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
        
        session = sessions{mouse_it}(session_it);

        for day = [1,2,3,4]
            
            trial_id = trials(day);

            file_name = strcat('mouse_',mouse,'_session_',session,'_trial_',trial_id,...
            '_v1.4.20.3.0.1.1.0.mat');

            model_name = strcat('model_mouse_',mouse,'_session_',session,'_trial_',trial_id,...
            '_v1.4.20.3.0.1.1.0','.mat');
            
            model_file_dir = strcat(MODEL_dir,model_name);

            heat_capacity_name = strcat('heat_capacity_mouse_',mouse,'_session_',session,...
                '_trial_',trial_id,'_v1.4.20.3.0.1.1.0','.mat');
            
            heat_capacity_dir = strcat(HC_dir, heat_capacity_name);

            load(strcat(data_dir,file_name));

            number_of_cells_flag = 0;
            for trial = 1:5 % iteration over trials
   
                data = calcium_binary.rest_trace{trial};
                data_mean  = mean(data,2);
                firing_rate{day,trial} = data_mean *10;
            end
        end
  
        
        figure(1)
        for day = [1,2,3,4]
            for trial = 1:5    
                subplot(4,5,(day-1)*5+trial)
                hist = histogram(firing_rate{day,trial})
                plot(hist.BinEdges(2:end),hist.Values)
                xlabel('MeanFiringRate [Hz]') 
                ylabel('#') 
                title(strcat('Day:',int2str(day),', trial:',int2str(trial)))
            end
        end
    end