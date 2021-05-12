
% Created on mrt 18 13:32
% Auth
% Created on mrt 19 16:03
% Author: Melisa

data_dir = '/home/melisa/Documents/criticality/data/calcium_activity_converted/'
file_name = 'mouse_32363_session_2_trial_11_v1.4.20.3.0.1.1.0.mat'


load(strcat(data_dir,file_name))

ncells = 15
random_samples = 25

combos = combntns(1:random_samples,ncells);
size_combos = size(combos);
sample = randsample(size_combos(1),random_samples);

   
for day = 1:5
    
x = calcium_binary.rest_binary{day};
x_size = size(x);
nneuros = x_size(1);
subplot(2,3,day)

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

scatter(-model_logprobs/ncells, entropy_states/ncells)
hold on;
clear entropy_states model_logprobs

end
end