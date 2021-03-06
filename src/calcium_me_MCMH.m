

% Created on mrt 19 16:03
% Author: Melisa

data_dir = '/home/melisa/Documents/criticality/data/calcium_activity_converted/'
result_dir = '/home/melisa/Documents/criticality/results/'
figures_dir = '/home/melisa/Documents/criticality/figures/'

mouse = '56165'
session = '1'
trial = '6'

file_name = strcat('mouse_',mouse,'_session_',session,'_trial_',trial,'_v1.4.20.3.0.1.1.0.mat')

load(strcat(data_dir,file_name));

ncells = 15
random_samples = 20

file_results_name = strcat('heat_capacity_mouse_',mouse,'_session_',session,'_trial_',trial,'_v1.4.20.3.0.1.1.0_ncells_',int2str(ncells),'.mat');

combos = combntns(1:random_samples,ncells);
size_combos = size(combos);
sample = randsample(size_combos(1),random_samples);

   
for trial = 1:5 % iteration over trials
    
    x = calcium_binary.rest_binary{trial};
    x_size = size(x);
    nneuros = x_size(1);
    
    for i=1:20

        model = maxent.createModel(ncells,'ising');
        model = maxent.trainModel(model,x(combos(sample(i),:),:),'threshold',1);

        limited_empirical_distribution = maxent.getEmpiricalModel(x(combos(mod(sample(i)*i,nneuros)+1,:),:),'min_count',2);

        model_logprobs = maxent.getLogProbability(model,limited_empirical_distribution.words);
        npatterns = size(model_logprobs); 

        for j=1:npatterns(2)
            entropy_states(j) = sum(model_logprobs(-model_logprobs<-model_logprobs(j)));
        end
        entropy_states = log(entropy_states);
        

        C_N = MCMH(ncells,model, 100000);
        heat_capacity{day,i} = C_N;

    end
end

results = strcat(result_dir,file_results_name);
save(results,'heat_capacity')


figure(1)
for day = 1:4    
subplot(2,3,day)
for i=1:20
plot(0.5:0.1:3,heat_capacity{day,i})
hold on
xlabel('T') 
ylabel('C_N') 
ylim([1 12])
end
end
subplot(2,3,2)
title(strcat('Heat Capacity Mouse:', mouse, ' session:', session))


figure_name = strcat(figures_dir, 'heat_capacity_mouse_',mouse,'_session_',session,'_trial_',trial,'_v1.4.20.3.0.1.1.0_ncells_,',ncells_s,'.png');
saveas(1,figure_name)



