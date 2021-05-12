
clear all
dir_path = '/scratch/melisa/pFC_Criticality/data/activity/';
MODEL_dir = '/scratch/melisa/pFC_Criticality/data/model/';
HC_dir = '/scratch/melisa/pFC_Criticality/data/heat_capacity/';
figures_dir = '/scratch/melisa/pFC_Criticality/figures/';

data_path = 'mouse_56165_session_1_trial_1_v1.4.20.3.0.1.1.0.mat';

load(strcat(dir_path,data_path))

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
index = find(firing_rate >4 & firing_rate<7);
if length(index) > 50
    spikes = all_data(index,:);
    spikes50 = spikes(1:50,:);
end

%%
% randomly divide it into a training set and a test set (so we can verify how well we trained)
[ncells,nsamples] = size(spikes50);
idx_train = randperm(nsamples,ceil(nsamples/1));
idx_test = setdiff(1:nsamples,idx_train);
samples_train = spikes50(:,idx_train);
%samples_test = spikes15(:,idx_test);

% create a k-pairwise model (pairwise maxent with synchrony constraints)
model = maxent.createModel(ncells,'ising');

% train the model to a threshold of one standard deviation from the error of computing the marginals.
% because the distribution is relatively small (15 dimensions) we can explicitly represent all 2^15 states 
% in memory and train the model in an exhaustive fashion.
model = maxent.trainModel(model,samples_train,'threshold',5);

n_cells = 50;

% Get parameters of the model
h_i = model.factors(1:n_cells);
w_ij = squareform(model.factors(n_cells+1:end))*0.5;

%%
Pool_N = 20;
p_cells = 10;

Pools = zeros(Pool_N,p_cells);

for pp = 1:Pool_N
    y = datasample(1:n_cells,p_cells,'Replace',false);
    Pools(pp,:) = y;
end
T = 1;
n_sweeps = 20000;
Generations_N = 100;
Fit_Epoch = zeros(Pool_N,Generations_N);
for generations = 1:Generations_N
    C_N = zeros(1,Pool_N);
    for gg = 1:Pool_N
        Pool_take = Pools(gg,:);
        h_i_g = h_i(Pool_take);
        w_ij_g = w_ij(Pool_take,Pool_take);
    
        patt = randi(2,p_cells,1)-1; % Initialize
        pp = patt; 
        E_state = sum(pp'.*h_i_g) + 1/2*(pp'*w_ij_g)*pp;
        P_stae = exp(-E_state/T);
        E_MH = [];
        kk=0;
        for t = 1:n_sweeps*p_cells
            pp = patt; 
            flip = randi(p_cells,1); % Pick a site to flip
            pp(flip) = -pp(flip)+1; % Flip (notice the 0,1 convention)
            E_new = sum(pp'.*h_i_g) + 1/2*(pp'*w_ij_g)*pp;
            P_new = exp(-E_new/T);
            take = P_new/P_stae; % METROPOLIS-HASTINGS
            if(rand(1)<take)
                patt = pp;
                P_stae = P_new;
                E_state = E_new;
            end
            if(mod(t,p_cells)==0 && t>p_cells*10) %% Sample energy every N steps and remove burn-in period
                kk = kk+1;
                E_MH(kk) = E_state;
            end
        end
        C_N(gg) = var(E_MH)/T^2;
        Fit_Epoch(gg,generations) = C_N(gg);
    end
    %EVOLVE POOLS 
    %Mix
    new_pool = zeros(size(Pools));
    for n_pool = 1:Pool_N
        % Pick 2 parents with probability proportional to their fit
        C_N_Prob = [0, cumsum(C_N/sum(C_N))];
        p1 = discretize(rand(1),C_N_Prob);
        p2 = discretize(rand(1),C_N_Prob);
        C1 = Pools(p1,:);
        C2 = Pools(p2,:);
        Uni = cell(2,1);
        Uni{1} = setdiff(C1,C2);
        Uni{2} = setdiff(C2,C1);     
        if numel(Uni{1})>0
        % Keep Common neurons 
            Common = intersect(C1,C2);
        % Complete the pool by selecting from the 2 parents proportionally to their
        % relative fit
            Add = [];
            for pick = 1:p_cells-numel(Common)
                nn = discretize(rand(1), [0, C_N(p1)/(C_N(p1)+C_N(p2)),1]);
                np = randi(numel(Uni{nn}),1);
                Add = cat(1,Add(:),Uni{nn}(np));
                Uni{nn}(np)=[];
            end
            Common = cat(2,Common,Add');
        else
            Common = C1;    
        end
        % Random Noise (substitute cells in the population with probability p)
        Take_From = setdiff(1:n_cells,Common);
        for ss = 1:p_cells
            if(rand(1)<0.02)
                Common(ss) = Take_From(randi(numel(Take_From),1));
            end
        end
        new_pool(n_pool,:) = Common;
    end
 
    Pools = new_pool;
    c(generations,:) = C_N;
end
