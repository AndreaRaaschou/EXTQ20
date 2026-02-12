function main_ecoevo(sigma_a_ind)

% This is the main script for the eco-evolutionary dynamics analysis

close all

rand('seed',fix(sum(1000*clock)));
randn('seed',fix(sum(2000*clock)+1));

sigma_a_vec=[0.1 0.3 0.5 0.7 1];
sigma_a=sigma_a_vec(sigma_a_ind);

%% Set up the parameters and initial conditions of the ecological model

% Habitat variables
K0=[10000];
U=[0];
sigma_K=[1];

% Consumer variables
N0=[1];
V=[-0];
prey_sp_id=[0];
prey_sp_id_counter=max(prey_sp_id);

m=length(V);
intrin_growth_comp=1;
r_A=ones(m,1)*intrin_growth_comp;
mut_N=1*10e-3;
sigma_mut_N = 0.02;

% Predator variables (NOTE: In this simulation we omit predators, i.e. we set predator abudnance to zero)
pred_sp_id=[0];
pred_sp_id_counter=max(pred_sp_id);

P0=[0];
Z=[0];

n=length(Z);
intrin_death_pred=-0.2;
mu_A=ones(n,1)*intrin_death_pred;
mut_P=1*10e-2;
sigma_mut_P = 0.02;
bmax=0.0001;
sigma_b=0.1;
cP=0.3;

% Implementation parameters
t_end = 1000000;
dt=0.5;
t_evo=0;
mut_range=-3:0.01:3;

%% Variables for data output (keep: used later)
prey_dist_data{1,1}=1;
prey_dist_data{1,2}=[V; N0; prey_sp_id];

pred_dist_data{1,1}=1;
pred_dist_data{1,2}=[Z; P0; pred_sp_id];

prey_fitland_data{1,1}=1;
prey_fitland_data{1,2}=0;

pred_fitland_data{1,1}=1;
pred_fitland_data{1,2}=0;

%% Plot the resource landscape
count=0;
range=-3:0.01:3;
for i=range
    count=count+1;
    K_A(count) = K0*exp(-(i-U).^2/2/sigma_K^2);
end

figure(1)
subplot(4,1,1)
plot(range,K_A);
title('resource distribution')

%% Set up the vectors for the ODE system (NOTE: We omit predators in this simulations but parameters are still sett to keep the GLV intact)
abun_vec=[N0; P0];
growth_death_vec=[r_A; mu_A];
mut_vec=[ones(m,1)*mut_N; ones(n,1)*mut_P];
trait_vec=[V'; Z'];

%% Save the parameters in a struct that will be piped into functions later
parameters.N0=N0;
parameters.P0=P0;
parameters.V=V;
parameters.Z=Z;
parameters.m=m;
parameters.n=n;
parameters.r_A=r_A;
parameters.mu_A=mu_A;
parameters.K0=K0;
parameters.U=U;
parameters.sigma_K=sigma_K;
parameters.sigma_a=sigma_a;
parameters.bmax=bmax;
parameters.sigma_b=sigma_b;
parameters.cP=cP;
parameters.t_end=t_end;
parameters.abun_vec=abun_vec;
parameters.growth_death_vec=growth_death_vec;
parameters.intrin_growth_comp=intrin_growth_comp;
parameters.intrin_death_pred=intrin_death_pred;
parameters.mut_N=mut_N;
parameters.mut_P=mut_P;
parameters.dt=dt;

%% Compute equilibrium community before evolution
mutant_trait=1;
mut_type_flag=1;

simsteps=t_end/dt;
[t,y,A,mut_fit, K_A, alpha, a] = mutfit_and_popequi_func(parameters, mutant_trait, mut_type_flag, simsteps, 0);

% Update abundances from equilibrium
abun_vec=y(end,:)';


% Update the system according to the latest simulated equilibrium
[N0,P0,m,n,V,Z,r_A,mu_A, abun_vec, growth_death_vec, mut_vec,trait_vec, prey_sp_id, pred_sp_id] = ...
    uppdate_sys_func(abun_vec,m,n,V,Z,intrin_growth_comp,intrin_death_pred,mut_N,mut_P, prey_sp_id, pred_sp_id);

% Update parameter struct
parameters.N0=N0;
parameters.P0=P0;
parameters.V=V;
parameters.Z=Z;
parameters.m=m;
parameters.n=n;
parameters.r_A=r_A;
parameters.mu_A=mu_A;
parameters.abun_vec=abun_vec;
parameters.growth_death_vec=growth_death_vec;

%% Evolutionary time loop
while t_evo < 5000

    % Compute fitness landscapes (N_fitland used later for saving/plotting)
    simsteps=3;
    count=1;
    for j=mut_range
        mutant_trait=j;

        mut_type_flag=1;
        [t,y,A,mut_fit, K_A, alpha, a] = mutfit_and_popequi_func(parameters, mutant_trait, mut_type_flag, simsteps, 1);
        N_fitland(count)=max(mut_fit,0);

        mut_type_flag=2;
        [t,y,A,mut_fit, K_A, alpha, a] = mutfit_and_popequi_func(parameters, mutant_trait, mut_type_flag, simsteps, 1);
        P_fitland(count)=max(mut_fit,0);

        count=count+1;
    end

    if t_evo==0
        prey_fitland_data{1,1}=t_evo;
        prey_fitland_data{1,2}=N_fitland;

        pred_fitland_data{1,1}=t_evo;
        pred_fitland_data{1,2}=P_fitland;
    end

    % Mutation selection rates
    w=abun_vec.*mut_vec;
    w_tot=sum(w);

    mut_sp=find(rand<cumsum(w/w_tot),1,'first');
    mut_fit=[];

   % if growth_death_vec(mut_sp)>0 % consumer mutation

        mutant_trait=trait_vec(mut_sp) + sigma_mut_N*randn;
        mut_type_flag=1;

        t_evo=t_evo+1;

        simsteps=3;
        [t,y,A,mut_fit, K_A, alpha, a] = mutfit_and_popequi_func(parameters, mutant_trait, mut_type_flag, simsteps, 1);

        if mut_fit>0

            disp('consumer mutated with positive fitness')

            V_tmp=V; V_tmp(mut_sp)=mutant_trait;
            Z_tmp=Z;

            parameters.V=V_tmp;

            simsteps=t_end/dt;
            [t,y,A,mut_fit, K_A, alpha, a] = mutfit_and_popequi_func(parameters, mutant_trait, mut_type_flag, simsteps, 0);

            abun_vec_tmp=y(end,:)';

            [mm nn]=size(y);
            for i=1:nn
                if length(find(y(:,i)<=0))>1
                    abun_vec(i)=0;
                end
            end

            [N0_tmp,P0_tmp,m_tmp,n_tmp,V_tmp,Z_tmp,r_A_tmp,mu_A_tmp, abun_vec_tmp, growth_death_vec_tmp, mut_vec_tmp,trait_vec_tmp, prey_sp_id_tmp, pred_sp_id_tmp] = ...
                uppdate_sys_func(abun_vec_tmp,m,n,V_tmp,Z_tmp,intrin_growth_comp,intrin_death_pred,mut_N,mut_P, prey_sp_id, pred_sp_id);

            parameters.V=V_tmp;
            parameters.Z=Z_tmp;
            parameters.m=m_tmp;
            parameters.n=n_tmp;
            parameters.r_A=r_A_tmp;
            parameters.mu_A=mu_A_tmp;
            parameters.abun_vec=abun_vec_tmp;
            parameters.growth_death_vec=growth_death_vec_tmp;
            parameters.N0=N0_tmp;
            parameters.P0=P0_tmp;

            % Mutual invasibility check
            simsteps=3;
            [t,y,A,mutual_inv, K_A, alpha, a] = mutfit_and_popequi_func(parameters, trait_vec(mut_sp), mut_type_flag, simsteps, 1);

            if mutual_inv > 0

                disp('consumer mutated with mutual invasibility')

                % Include the mutant and keep the mutating species
                N0_tmp=[1; N0]; 
                V_tmp=[mutant_trait V]; 
                m_tmp=length(V_tmp); 
                r_A_tmp=ones(m_tmp,1)*intrin_growth_comp;

                n_tmp=length(Z); 
                prey_sp_id_tmp=[prey_sp_id(mut_sp) prey_sp_id]; 
                pred_sp_id_tmp=pred_sp_id;

                abun_vec_tmp=[N0_tmp; P0];
                growth_death_vec_tmp=[r_A_tmp; mu_A];
                mut_vec_tmp=[ones(m_tmp,1)*mut_N; ones(n,1)*mut_P];
                trait_vec_tmp=[V_tmp'; Z'];

                parameters.V=V_tmp;
                parameters.Z=Z;
                parameters.m=m_tmp;
                parameters.n=n_tmp;
                parameters.r_A=r_A_tmp;
                parameters.mu_A=mu_A;
                parameters.abun_vec=abun_vec_tmp;
                parameters.growth_death_vec=growth_death_vec_tmp;
                parameters.N0=N0_tmp;
                parameters.P0=P0;

                simsteps=(t_end/dt)*100;
                [t,y,A,mutual_inv, K_A, alpha, a] = mutfit_and_popequi_func(parameters, trait_vec(mut_sp), mut_type_flag, simsteps, 0);

                abun_vec_tmp=y(end,:)';

                [mm nn]=size(y);
                for i=1:nn
                    if length(find(y(:,i)<=0))>1
                        abun_vec(i)=0;
                    end
                end
                Z_tmp=Z;

                [N0_tmp,P0_tmp,m_tmp,n_tmp,V_tmp,Z_tmp,r_A_tmp,mu_A_tmp, abun_vec_tmp, growth_death_vec_tmp, mut_vec_tmp,trait_vec_tmp,prey_sp_id_tmp, pred_sp_id_tmp] = ...
                    uppdate_sys_func(abun_vec_tmp,m_tmp,n_tmp,V_tmp,Z_tmp,intrin_growth_comp,intrin_death_pred,mut_N,mut_P,prey_sp_id_tmp, pred_sp_id_tmp);

            end

            % Update the system
            N0=N0_tmp; P0=P0_tmp;
            V=V_tmp; Z=Z_tmp;
            abun_vec=abun_vec_tmp; trait_vec=trait_vec_tmp;
            m=m_tmp; n=n_tmp;
            r_A=r_A_tmp; mu_A=mu_A_tmp;
            growth_death_vec=growth_death_vec_tmp;
            mut_vec=mut_vec_tmp;
            prey_sp_id=prey_sp_id_tmp; pred_sp_id=pred_sp_id_tmp;


        %end
    end

    % Save data 
    prey_dist_data{end+1,1}=t_evo;
    prey_dist_data{end,2}=[V; N0'; prey_sp_id];

%     com_info_data{end+1,1}=t_evo;
%     com_info_data{end,2}=K_A;
%     com_info_data{end,3}=alpha;
%     com_info_data{end,4}=a;

    prey_fitland_data{end+1,1}=t_evo;
    prey_fitland_data{end,2}=N_fitland;

    % Update parameters 
    parameters.V=V;
    parameters.Z=Z;
    parameters.m=m;
    parameters.n=n;
    parameters.r_A=r_A;
    parameters.mu_A=mu_A;
    parameters.abun_vec=abun_vec;
    parameters.growth_death_vec=growth_death_vec;
    parameters.N0=N0;
    parameters.P0=P0;

    % Plot every 50 steps
    if mod(t_evo,50)==1
        figure(1)
        subplot(4,1,2)
        cla
        hold on
        plot(mut_range,N_fitland,'color',[0.4 0.4 0.4])
        ylabel('Cons. Fitness')
        ylim([0 0.01])

        subplot(4,1,3)
        hold on
        plot(parameters.V,t_evo,'.','color',[0.4 0.4 0.4])
        xlim([-3 3])
        ylabel('Evo. time')
        xlabel('Trait value')

        subplot(4,1,4)
        cla
        hold on
        for i=1:length(V)
            bar(V(i),N0(i),'facecolor',[0.4 0.4 0.4],'barwidth',0.1)
        end
        for i=1:length(Z)
            bar(Z(i),P0(i),'facecolor','r','barwidth',0.1)
        end
        xlim([-3 3])
        ylabel('N*')
        xlabel('Trait')
        title('Trait distribution')
    end

end


end
