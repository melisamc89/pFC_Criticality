%Created on mrt 19
%Author Federico Stella
%IMPORTANT, MODELS ARE ESTIMATED WITH THE sigma=[0,1] CONVENTION
%THis function estimates heat capactiy of an ising model using MCMH
%inputs : n_cells -> number of cells
%         model   -> ising model from ME toolboox
%         n_sweeps -> number of monte carlo iterations
% output : C_N : heat capacity in the range [0.5,0.75,0.8,0.9,1,1.1,1.2,1.25,1.50,1.75,2,2.25,2.50] 'temperature'

function C_N = MCMH(n_cells, model, n_sweeps)
    t_i = 0;
    C_N = [];
    for T = [0.5,0.75,0.8,0.9,1,1.1,1.2,1.25,1.50,1.75,2,2.25,2.50]
        t_i = t_i +1; 

        % Get parameters of the model
        h_i = model.factors(1:n_cells);
        w_ij = squareform(model.factors(n_cells+1:end));

        patt = randi(2,n_cells,1)-1; % Initialize
        pp = patt; 

        E_state = sum(pp'.*h_i) + 1/2*(pp'*w_ij)*pp;
        P_stae = exp(-E_state/T);

        E_MH = [];
        kk=0;
        for t = 1:n_sweeps*n_cells
            pp = patt; 
            flip = randi(n_cells,1); % Pick a site to flip
            pp(flip) = -pp(flip)+1; % Flip (notice the 0,1 convention)

            E_new = sum(pp'.*h_i) + 1/2*(pp'*w_ij)*pp;
            P_new = exp(-E_new/T);
            take = P_new/P_stae; % METROPOLIS-HASTINGS

            if(rand(1)<take)
                patt = pp;
                P_stae = P_new;
                E_state = E_new;

            end


            if(mod(t,n_cells)==0 && t>n_cells*10) %% Sample energy every N steps and remove burn-in period
                kk = kk+1;
                E_MH(kk) = E_state;
            end

        end

    C_N(t_i) = var(E_MH)/T^2;

    end

end