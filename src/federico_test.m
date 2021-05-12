

% Created on mrt 19 16:03
% Author: Melisa

data_dir = '/home/melisa/Documents/criticality/data/calcium_activity_converted/'
result_dir = '/home/melisa/Documents/criticality/results/'
figures_dir = '/home/melisa/Documents/criticality/figures/'

mouse = '56165'
session = '1'
trial = '6'

file_name = strcat('mouse_',mouse,'_session_',session,'_trial_',trial,'_v1.4.20.3.0.1.1.0_10hz.mat')

load(strcat(data_dir,file_name));

ncells = 20
ncells_s = '15'
random_samples = 25

file_results_name = strcat('heat_capacity_mouse_',mouse,'_session_',session,'_trial_',trial,'_v1.4.20.3.0.1.1.0_ncells_',ncells_s,'.mat');

combos = combntns(1:random_samples,ncells);
size_combos = size(combos);
sample = randsample(size_combos(1),random_samples);

   clear h_i_ref;
   clear w_ij_ref;
for day = [1 2 3 4 5]
    
x = horzcat(calcium_binary.rest_binary{day:day+0});
x_size = size(x);
nneuros = x_size(1);


for i=1:1

model = maxent.createModel(ncells,'ising');
model = maxent.trainModel(model,x(combos(sample(i),:),:),'threshold',1);

limited_empirical_distribution = maxent.getEmpiricalModel(x(combos(mod(sample(i)*i,nneuros)+1,:),:),'min_count',2);




h_i_ref(day,:) = model.factors(1:ncells);
w_ij_ref(day,:) = model.factors(ncells+1:end);


end
end

%%
for d1 = [1 2 3 4 5]
for dd = [1 2 3 4 5]
h_i_1 = h_i_ref(d1,:);
w_ij_1 = w_ij_ref(d1,:);

h_i_2 = h_i_ref(dd,:);
w_ij_2 = w_ij_ref(dd,:);


figure(1)
subplot(5,5,(d1-1)*5+dd)
scatter(w_ij_1,w_ij_2,50,'b','filled')
axis equal
refline(1,0)
lsline
figure(2)
subplot(5,5,(d1-1)*5+dd)
scatter(h_i_1,h_i_2,50,'b','filled')
axis equal
refline(1,0)
lsline
end
end