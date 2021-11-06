% NOTE for users:
%  This script produces Supp Figure 3, the first row of which also
%  corresponds to the data plots in main text Figure 5C,E,F. It also
%  produces a figure showing the plot in main text Figure 5D.

clear all;
addpath('functions'); %keep helper functions here

do_calc = 1;
do_save = 0; savefilename = ['Fig5_and_suppFigure3_general_' char(datetime('now','Format','yyyy-MM-dd_HH_mm_ss')) '.mat'];
do_plot = 1;

if do_calc

    sim = struct();

    % -- basic properties of objects, rewards, and timing
    sim.obj.n = 9;

    sim.rew.max = 5;
    sim.rew.all = 0:sim.rew.max;
    sim.rew.n = numel(sim.rew.all);

    sim.t.explore = 0.5; % time to explore one object
    sim.t.consume = 1.25; % time to consume one reward
    sim.t.travel = 8; % time to travel to next trial
    sim.t.max_trial_duration = sim.t.explore.*sim.obj.n + sim.t.consume.*sim.obj.n + sim.t.travel;

    % exponentially decaying reward sizes (with time spent in envir)
    % decay factor = 1 - (fraction of reward that disappears per unit time since start of trial)
    sim.t.reward_size_decay_factor = 0.725;

    
    % -- models of foraging behavior (for which we'll be calculating the optimal policy)
    sim.model = struct();
    sim.model.name = {'nexp','e and nexp','r and nexp','optimal'};
    sim.model.n = numel(sim.model.name);
    
    
    % -- environments
    sim.env.name = {'E1','E2','E3','E4','E5'};
    sim.env.n = numel(sim.env.name);

    % probability of entering each environment at the start of the trial
    sim.env.penter = ones(1,sim.env.n) ./ sim.env.n; % uniform probability
    assert(all(sim.env.penter >= 0) && sum(sim.env.penter) == 1,'environment-entering probabilities must form a proper probability distribution');
    
    % -- environment parameter sets
    sim.env.param = struct();
    
    % environments can differ in:
    % - max baiting probability (for best object in environment)
    % - exponentially decaying baiting probability (with object #)
    % - reward distributions for baited objects (p(rew size | baited))
    sim.env.param.set = {};
    
    min_prew_weight = 0.0005;
    
    % the different sets of environments
    scale_factor = 1;
    sim.env.param.set{end+1} = struct('paramset_name','pbait max', ...
        'pbaitmax',[.1 .3 .5 .7 .9], ...
        'tau',ones(1,sim.env.n).*0.7.*scale_factor, ...
        'prew_given_bait',{arrayfun(@(ei) (sim.rew.all > 0)./sum(sim.rew.all > 0),1:sim.env.n,'uniform',0)});
    
    sim.env.param.set{end+1} = struct('paramset_name','pbait decay rate', ...
        'pbaitmax',ones(1,sim.env.n).*0.4, ...
        'tau', [0.900 0.618 0.424 0.291 0.200].*scale_factor, ...
        'prew_given_bait',{arrayfun(@(ei) (sim.rew.all > 0)./sum(sim.rew.all > 0),1:sim.env.n,'uniform',0)});

    sim.env.param.set{end+1} = struct('paramset_name','rewsize distrib', ...
        'pbaitmax',ones(1,sim.env.n).*0.4, ...
        'tau',ones(1,sim.env.n).*0.7.*scale_factor, ...
        'prew_given_bait',{arrayfun(@(z) normalize([0 max(min_prew_weight,linspace((1-z),z,sim.rew.n-1).^2)],1),linspace(0,1,sim.env.n),'uniform',0)});

    sim.env.param.set{end+1} = struct('paramset_name','rewsize distrib, fixed prew', ...
        'pbaitmax',ones(1,sim.env.n).*0.2, ...
        'tau',zeros(1,sim.env.n), ...
        'prew_given_bait',{arrayfun(@(z) normalize( ...
              ([0 max(min_prew_weight,linspace((1-z),z,sim.rew.n-1).^2)])*(scale_factor) ...
            + [0 ones(1,sim.rew.n-1)]*mean(max(min_prew_weight,linspace((1-z),z,sim.rew.n-1).^2))*(1-scale_factor) ...
            ,1),linspace(0,1,sim.env.n),'uniform',0)});
    
    % since we are changing all 3 types of parameters here, we only change 
    % each one by 1/3 as much as before
    all3_scale_factor = (1/3);
	sim.env.param.set{end+1} = struct('paramset_name','ALL THREE', ...
        'pbaitmax',linspace(.5-.4*all3_scale_factor,.5+.4*all3_scale_factor,5), ...
        'tau',exp(linspace(-0.8574 + (1.5041)*0.5*all3_scale_factor,-0.8574 - (1.5041)*0.5*all3_scale_factor,sim.env.n)).*scale_factor, ...
        'prew_given_bait',{arrayfun(@(z) normalize( ...
              ([0 max(min_prew_weight,linspace((1-z),z,sim.rew.n-1).^2)])*(all3_scale_factor) ...
            + [0 ones(1,sim.rew.n-1)]*mean(max(min_prew_weight,linspace((1-z),z,sim.rew.n-1).^2))*(1-all3_scale_factor) ...
            ,1),linspace(0,1,sim.env.n),'uniform',0)});

    % control condition
    sim.env.param.set{end+1} = struct('paramset_name','FIXED PREW DISTRIB AND REW SIZE', ...
        'pbaitmax',ones(1,sim.env.n).*0.2, ...
        'tau',zeros(1,sim.env.n), ...
        'prew_given_bait',{arrayfun(@(z) normalize([0 max(min_prew_weight,[0 0 1 0 0])],1),linspace(0,1,sim.env.n),'uniform',0)});

    sim.env.param.n = numel(sim.env.param.set);
    
    for epsi = 1:sim.env.param.n
        % resulting baiting probability distributions, in the format pbait_for_obj{ei}(oi)
        sim.env.param.set{epsi}.pbait_for_obj = arrayfun(@(ei) sim.env.param.set{epsi}.pbaitmax(ei) .* 2.^(-sim.env.param.set{epsi}.tau(ei).*((1:sim.obj.n)-1)),1:sim.env.n,'uniform',0);

        % resulting reward probability distributions, in the format: prewsize(rewid,objid)
        sim.env.param.set{epsi}.prew_for_obj = arrayfun(@(ei) sim.env.param.set{epsi}.prew_given_bait{ei}' * sim.env.param.set{epsi}.pbait_for_obj{ei},1:sim.env.n,'uniform',0);
        sim.env.param.set{epsi}.prew_for_obj = cellfun(@(prew) [1-sum(prew(2:end,:),1) ; prew(2:end,:)],sim.env.param.set{epsi}.prew_for_obj,'uniform',0);
    end
    
    % results of setting these environmental params
    % (i.e. resulting settings of rstate and state)
    sim.env.param.result = cell(1,sim.env.param.n);

    % estimate results for specific reward sequences in the task
    median_rew = median(sim.rew.all(sim.rew.all > 0));
    sim.task.rewseq.name = {'Norew','Rew'};
    sim.task.rewseq.color = {'b','r'};
    sim.task.rewseq.linestyle = {'-',':'};
    sim.task.rewseq.r = {zeros(1,sim.obj.n),[median_rew zeros(1,sim.obj.n-1)]};
    sim.task.rewseq.n = numel(sim.task.rewseq.name);

    
    % --- MAIN LOOP
    
    % do simulation for each set of environment parameters
    for epsi = 1:sim.env.param.n
        
        % set environment parameters for this simulation
        efnames = fieldnames(sim.env.param.set{epsi});
        for efni = 1:numel(efnames)
            efn = efnames{efni};
            sim.env.(efn) = sim.env.param.set{epsi}.(efn);
        end

        % -- rstates: sub-states defined by the explored reward sequence
        sim.rstate = struct();
        sim.rstate.nexp = 0;
        sim.rstate.r = nans(1,sim.obj.n);

        for oi = 1:sim.obj.n
            for ri = 1:sim.rew.n
                okprev = sim.rstate.nexp == oi-1 ...
                    & (isnan(sim.rstate.r(:,1)) | sim.rstate.r(:,1) <= sim.rew.all(ri));
                rnew = [sim.rew.all(ri)*ones(sum(okprev),1) sim.rstate.r(okprev,1:(end-1))];

                sim.rstate.nexp = [sim.rstate.nexp ; oi.*ones(size(rnew,1),1)];
                sim.rstate.r = [sim.rstate.r ; rnew];
            end
        end

        sim.rstate.n = size(sim.rstate.r,1);

        sim.rstate.next = nans(sim.rstate.n,sim.rew.n);
        for si = 1:sim.rstate.n
            if sim.rstate.nexp(si) < sim.obj.n
                curr = sim.rstate.r(si,:);
                curr = curr(~isnan(curr));

                for ri = 1:sim.rew.n
                    nextr = sort([curr sim.rew.all(ri)],2,'descend');
                    nextstateid = find(sim.rstate.nexp == sim.rstate.nexp(si)+1 & rowis(sim.rstate.r(:,1:numel(nextr)),nextr));
                    assert(numel(nextstateid)==1,'must find unique next state for each (current state,explored reward) pair');
                    sim.rstate.next(si,ri) = nextstateid;
                end
            end
        end

        % -- states = (environment, explored reward sequence)
        sim.state = struct();
        sim.state.nexp = [];
        sim.state.nexp = repmat(sim.rstate.nexp,[sim.env.n 1]);
        sim.state.r = repmat(sim.rstate.r,[sim.env.n 1]);
        sim.state.ri = repmat((1:sim.rstate.n)',[sim.env.n 1]); % index into 'rstates' for this reward history sub-state
        sim.state.env = arrayfun(@(ei) ei*ones(sim.rstate.n,1),1:sim.env.n,'uniform',0);
        sim.state.env = vertcat(sim.state.env{:});
        sim.state.next = arrayfun(@(ei) sim.rstate.next + sim.rstate.n*(ei-1),1:sim.env.n,'uniform',0);
        sim.state.next = vertcat(sim.state.next{:});

        sim.state.n = size(sim.state.r,1);


        % p(next revealed reward size | state, choice to explore)
        sim.state.pnextrew_if_exp = nans(sim.state.n,sim.rew.n);
        for ei = 1:sim.env.n
            for oi = 1:sim.obj.n
                % states in environment ei where are about to explore object oi
                okstates = sim.state.env == ei & sim.state.nexp == (oi-1);
                pnextrew = sim.env.prew_for_obj{ei}(:,oi)';
                sim.state.pnextrew_if_exp(okstates,:) = repmat(pnextrew,[sum(okstates) 1]);
            end
        end

        assert(all(all(isnan(sim.state.next) | sim.state.next > repmat((1:sim.state.n)',[1 sim.rew.n]))),'states must only transition to other states with a higher state index (i.e. MDP must not be recursive) for our current calculations to work');

        % p(reach this state during a trial | always explore)
        sim.state.preach_state_if_always_explore = zeros(sim.state.n,1);
        % initialize with probability of starting in each state
        for ei = 1:sim.env.n
            % find start state in this environment
            si = find(sim.state.env == ei & sim.state.nexp == 0);
            assert(numel(si)==1,['could not find unique starting state for environment ' num2str(ei)]);
            sim.state.preach_state_if_always_explore(si) = sim.env.penter(ei);
        end
        % for each state, add its outgoing probability of transitioning to each other state
        for si = 1:sim.state.n
            snext = sim.state.next(si,:);
            num_ok_snext = sum(~isnan(snext));
            if num_ok_snext > 0
                assert(num_ok_snext == numel(snext),'this state has at least one valid successor but oddly also has some invalid successors');
                % calc overall probability of being in this state,
                % and then reaching each possible next state (if always explore)
                pnext = sim.state.preach_state_if_always_explore(si) .* (sim.state.pnextrew_if_exp(si,:))';
                % add it to those states' probabilities of being reached
                sim.state.preach_state_if_always_explore(snext) = sim.state.preach_state_if_always_explore(snext) + pnext;
            end
        end


        % Bayesian beliefs about current environment if its identity is not
        % known (computed under assumption that you have always explored 
        % thus far this trial). Used for POMDP implementation of the
        % 'r and nexp' model.

        % p(environment = ei | reward history = ri)
        for ri = 1:sim.rstate.n
            oks = sim.state.ri == ri;
            ps = sim.state.preach_state_if_always_explore(oks);
            es = sim.state.env(oks);

            pe = arrayfun(@(ei) sum(ps .* (es==ei)),1:sim.env.n);
            pe = pe ./ sum(pe);

            sim.rstate.penv(ri,:) = pe;
            sim.state.penv_given_r(oks,:) = repmat(pe,[sum(oks) 1]);
        end

        % p(next revealed reward size | state with reward history = ri, environment = unknown)
        sim.state.pnextrew_if_exp_and_unknown_env = nans(sim.state.n,sim.rew.n);
        for ri = 1:sim.rstate.n
            si = find(sim.state.ri == ri);
            assert(isequal(sim.state.env(si),(1:sim.env.n)'),'each reward history must have exactly one corresponding state in each environment');

            % p(env) {1 x Nenv matrix} x p(rew | env) {Nenv x Nrew matrix}
            % = p(rew) {1 x Nrew matrix}
            prew = sim.rstate.penv(ri,:) * sim.state.pnextrew_if_exp(si,:);

            sim.state.pnextrew_if_exp_and_unknown_env(si,:) = repmat(prew,[sim.env.n 1]);
        end
        
        % save resulting rstate and state settings for the current
        % environment parameter set
        sim.env.param.result{epsi} = struct('rstate',sim.rstate,'state',sim.state);


        % find each model's optimal policy to handle this environment
        for mi = 1:sim.model.n

            sim.calc = struct();
            
            sim.calc.maxiter = 40;
            sim.calc.stoppingrule.min_delta_rewrate = 10^-7; %10^-8;
            sim.calc.niter = nan;

            sim.calc.initial_rewrate = 0;
            sim.calc.calcdur = nans(1,sim.calc.maxiter);
            
            sim.calc.assumed_rewrate = nans(1,sim.calc.maxiter); % average reward rate assumed in calculations
            sim.calc.rewrate = nans(1,sim.calc.maxiter); % actual average reward rate if follow its policy

            sim.calc.vcon = nans(sim.state.n,sim.calc.maxiter); % value if consume
            sim.calc.ncon = nans(sim.state.n,sim.calc.maxiter); % number of rewards consumed if consume
            sim.calc.rcon = nans(sim.state.n,sim.calc.maxiter); % amount of reward consumed if consume
            sim.calc.rnotcon = nans(sim.state.n,sim.calc.maxiter); % amount of reward found but NOT consumed if consume

            sim.calc.vexp = nans(sim.state.n,sim.calc.maxiter); % value if explore

            sim.calc.v = nans(sim.state.n,sim.calc.maxiter); % state value
            sim.calc.besta = nans(sim.state.n,sim.calc.maxiter); %optimal axn (1=consume, 2=explore)

            sim.calc.preach_state = nans(sim.state.n,sim.calc.maxiter); % p(reach state during trial)
            sim.calc.pstop_state = nans(sim.state.n,sim.calc.maxiter); % p(stop trial in this state), i.e. p(choose to consume)

            sim.calc.fieldnames_to_save = {'assumed_rewrate','rewrate','vcon','ncon','rcon','rnotcon','vexp','v','besta','preach_state','pstop_state'};

            % order of visiting states when calculating Vexp using dynamic
            % programming. This is different for different models to help
            % them apply their respective constraints.
            state_visit_order = sim.state.n:-1:1;
            state_visits_to_apply_constraints = false(sim.state.n,1);
            switch sim.model.name{mi}
                case 'nexp'
                    state_visit_order = [];
                    state_visits_to_apply_constraints = false(sim.state.n,1);
                    for nexp = sim.obj.n:-1:0
                        si = find(sim.state.nexp == nexp);
                        state_visit_order = [state_visit_order si(end:-1:1)'];
                        state_visits_to_apply_constraints(si(1)) = true;
                    end
                case 'e and nexp'
                    higher_nexp = [false ; (sim.state.nexp(2:end) > sim.state.nexp(1:(end-1)))];
                    same_env = [false ; (sim.state.env(2:end) == sim.state.env(1:(end-1)))];
                    state_visits_to_apply_constraints = higher_nexp & same_env;
            end
            assert(isequal(1:sim.state.n,sort(state_visit_order)),'in DP calculation, all states must be visited exactly once!');

            for ci = 1:sim.calc.maxiter
                curtimer = tic;

                % set 'average reward rate' parameter
                if ci == 1
                    % for first iter, set to initial setting
                    sim.calc.assumed_rewrate(ci) = sim.calc.initial_rewrate;
                else
                    % otherwise, set to result of last iteration
                    sim.calc.assumed_rewrate(ci) = sim.calc.rewrate(ci-1);
                end

                % calculate Vconsume (value of consumption)
                for oi = 0:sim.obj.n
                    okstate = sim.state.nexp == oi;
                    nokstate = sum(okstate);

                    % for each number of found rewards we might potentially consume...
                    % (including 0)

                    % time spent by the time finished consuming them
                    time_spent = (oi.*sim.t.explore + (0:oi).*sim.t.consume);

                    % decay factor for each potentially consumed reward
                    rdecay = sim.t.reward_size_decay_factor.^time_spent;
                    % sum of (partially decayed) reward
                    rew_from_trial = cumsum([zeros(nokstate,1) sim.state.r(okstate,1:oi)] .* repmat(rdecay,[nokstate 1]),2);

                    % time remaining from the fixed-duration time segment
                    % that we're considering, after travel time is taken into account
                    time_remaining_after_travel = sim.t.max_trial_duration - time_spent - sim.t.travel;

                    % expected reward gained during remaining time
                    rew_from_time_remaining = time_remaining_after_travel .* sim.calc.assumed_rewrate(ci);

                    % expected reward rate during the fixed-duration time segment
                    rewrate = (rew_from_trial + repmat(rew_from_time_remaining,[nokstate 1])) ./ sim.t.max_trial_duration;

                    % find optimal reward rate, and the number of rewards we should 
                    % choose to consume to achieve it
                    % NOTE: it is "optimal_nconsumed_plus_one" because if it is 1, 
                    %  (i.e. indexing the first column of rewrate), that correponds to
                    %  the case where we consume 0 of the rewards
                    [optimal_rewrate,optimal_nconsumed_plus_one] = max(rewrate,[],2);

                    sim.calc.vcon(okstate,ci) = optimal_rewrate;
                    sim.calc.ncon(okstate,ci) = optimal_nconsumed_plus_one - 1;
                    sim.calc.rcon(okstate,ci) = rew_from_trial(sub2ind(size(rew_from_trial),(1:nokstate)',optimal_nconsumed_plus_one));
                    sim.calc.rnotcon(okstate,ci) = max(rew_from_trial,[],2) - sim.calc.rcon(okstate,ci);
                end

                % calculate Vexplore (value of exploration)
                % via dynamic programming
                % and then find each state's optimal action 

                % end states: already explored all objects, so cannot explore
                % (enforced by setting vexp = -inf). Optimal action is to consume. 
                okstate = sim.state.nexp == sim.obj.n;
                sim.calc.vexp(okstate,ci) = -inf;
                sim.calc.v(okstate,ci) = sim.calc.vcon(okstate,ci);
                sim.calc.besta(okstate,ci) = 1;

                % all other states... 
                for si = state_visit_order
                    % if it is an end state, we already calculated its value of
                    % exploration and optimal policy (because we always
                    % consume)
                    if sim.state.nexp(si) == sim.obj.n
                        continue;
                    end


                    if strcmp(sim.model.name{mi},'r and nexp')
                        % POMDP ALGORITHM for 'r and nexp' model
                        % which does not know the current environment
                        % and hence has to predict p(next reward) by
                        % marginalizing over its current belief about the
                        % environment's identity
                        pnextrew = sim.state.pnextrew_if_exp_and_unknown_env(si,:);
                    else
                        % STANDARD ALGORITHM
                        pnextrew = sim.state.pnextrew_if_exp(si,:);
                    end

                    % calc expected value of exploration, considering all possible
                    % reward sizes that could be revealed
                    vnextrew = sim.calc.v(sim.state.next(si,:),ci);
                    sim.calc.vexp(si,ci) = pnextrew * vnextrew;

                    if sim.calc.vcon(si,ci) >= sim.calc.vexp(si,ci)
                        % consumption is better (or equal), so set it as optimal policy
                        sim.calc.besta(si,ci) = 1;
                        sim.calc.v(si,ci) = sim.calc.vcon(si,ci);
                    else
                        % exploration is better, so set it as optimal policy
                        sim.calc.besta(si,ci) = 2;
                        sim.calc.v(si,ci) = sim.calc.vexp(si,ci);
                    end

                    % in certain models, we constrain certain sets of states 
                    % to always use the same actions
                    if state_visits_to_apply_constraints(si)
                        switch sim.model.name{mi}
                            case 'e and nexp'
                                % enforce model constraint to only consider
                                %  environment and # of explored fractals, 
                                %  ignoring rewards seen so far
                                % by forcing all states in this environment with 
                                %  the same # explored actions to take the same
                                %  action
                                oks = sim.state.env == sim.state.env(si) & sim.state.nexp == sim.state.nexp(si);
                            case 'nexp'
                                % enforce model constraint to only consider
                                %  # of explored fractals, ignoring environment
                                %  identity and rewards seen so far
                                oks = sim.state.nexp == sim.state.nexp(si);
                            otherwise 
                                error(['do not know how to apply constraints to model "' sim.model.name{mi} '"']);
                        end

                        cur_vcon = sim.calc.vcon(oks,ci);
                        cur_vexp = sim.calc.vexp(oks,ci);
                        cur_ps = sim.state.preach_state_if_always_explore(oks);
                        cur_ps = cur_ps ./ sum(cur_ps);

                        % E[V(action)] = sum_s p(state = s | environment,#explored)*V(action | state = s)
                        cur_evcon = sum(cur_ps .* cur_vcon);
                        cur_evexp = sum(cur_ps .* cur_vexp);

                        if cur_evcon >= cur_evexp
                            % consumption is better (or equal), so consume
                            sim.calc.besta(oks,ci) = 1;
                            sim.calc.v(oks,ci) = cur_vcon;
                        else
                            % exploration is better, so explore
                            sim.calc.besta(oks,ci) = 2;
                            sim.calc.v(oks,ci) = cur_vexp;
                        end

                        % save these calculations
                        if ~isfield(sim.calc,'mod')
                            sim.calc.mod = cell(1,sim.model.n);
                        end
                        if ~isfield(sim.calc.mod{mi},'evcon')
                            sim.calc.mod{mi}.evcon = nans(sim.state.n,sim.calc.maxiter);
                            sim.calc.mod{mi}.evexp = nans(sim.state.n,sim.calc.maxiter);
                        end
                        sim.calc.mod{mi}.evcon(oks,ci) = cur_evcon;
                        sim.calc.mod{mi}.evexp(oks,ci) = cur_evexp;
                    end
                end

                % calculate true average reward rate if follow this policy
                sim.calc.preach_state(:,ci) = 0;
                sim.calc.pstop_state(:,ci) = 0;

                % initialize with probability of starting in each state
                for ei = 1:sim.env.n
                    % find start state in this environment
                    si = find(sim.state.env == ei & sim.state.nexp == 0);
                    assert(numel(si)==1,['could not find unique starting state for environment ' num2str(ei)]);
                    sim.calc.preach_state(si,ci) = sim.env.penter(ei);
                end

                for si = 1:sim.state.n
                    % stop if choose to consume rewards found so far
                    % otherwise, continue exploring
                    pstop = (sim.calc.besta(si,ci)==1);

                    % calc overall probability of being in this state, 
                    % and stopping here
                    sim.calc.pstop_state(si,ci) = sim.calc.preach_state(si,ci).*pstop;

                    if pstop < 1
                        % calc overall probability of being in this state,
                        % choosing to explore,
                        % and reaching each possible next state
                        snext = sim.state.next(si,:);
                        pnext = sim.calc.preach_state(si,ci) .* (1-pstop) .* (sim.state.pnextrew_if_exp(si,:))';
                        % add it to those states' probabilities of being reached
                        sim.calc.preach_state(snext,ci) = sim.calc.preach_state(snext,ci) + pnext;
                    end

                end

                assert(abs(sum(sim.calc.pstop_state(:,ci)) - 1) < 10^-9,'probability distribution of end states does not quite sum to 1');

                % true reward rate = expected reward rate based on probability of
                % reaching & then stopping in each state, and value gained from doing
                % so
                sim.calc.rewrate(ci) = sum(sim.calc.pstop_state(:,ci) .* sim.calc.vcon(:,ci));

                sim.calc.calcdur(ci) = toc(curtimer);
                disp([' envparam ' num2str(epsi) '/' num2str(sim.env.param.n) ' model ' num2str(mi) '/' num2str(sim.model.n) ' iter ' num2str(ci) '/' num2str(sim.calc.maxiter) ', rewrate = ' roundstr(.00001,sim.calc.rewrate(ci)) ', elapsed = ' roundstr(.1,sum(sim.calc.calcdur(1:ci))) ' s']);

                % test stopping rule

                sim.niter = ci;

                delta_rewrate = abs(sim.calc.rewrate(ci) - sim.calc.assumed_rewrate(ci));
                if delta_rewrate < sim.calc.stoppingrule.min_delta_rewrate
                    disp([' stopping, since this iteration rewrate changed by less than ' num2str(sim.calc.stoppingrule.min_delta_rewrate)]);
                    break;
                end
            end

            % save results of 'optimal' policy reached in final iteration
            for fni = 1:numel(sim.calc.fieldnames_to_save)
                fn = sim.calc.fieldnames_to_save{fni};
                fsize = size(sim.calc.(fn));
                assert(fsize(end) == sim.calc.maxiter,'expected data field to have one entry for each iteration');
                if mi == 1 && epsi == 1
                    assert(~isfield(sim.model,fn),['calculation data field "' fn '" is already in the "model" structure!']);
                    sim.model.(fn) = nans([fsize(1:(end-1)) sim.model.n sim.env.param.n]);
                end
                switch numel(fsize) % there must be a more elegant way to do this...
                    case 2
                        sim.model.(fn)(:,mi,epsi) = sim.calc.(fn)(:,sim.niter);
                    case 3
                        sim.model.(fn)(:,:,mi,epsi) = sim.calc.(fn)(:,:,sim.niter);
                    case 4
                        sim.model.(fn)(:,:,:,mi,epsi) = sim.calc.(fn)(:,:,:,sim.niter);
                    otherwise
                        error('only have code to select & copy data fields with 2 <= # dimensions <= 4');
                end
            end

            % estimate results for specific reward sequences in the task
            if mi == 1 && epsi == 1
                assert(~isfield(sim.model,'rewseq'),'field "rewseq" is already in the "model" structure');
                sim.model.rewseq.nexp = nans(sim.env.n,sim.task.rewseq.n,sim.model.n,sim.env.param.n);
            end
            for rsi = 1:sim.task.rewseq.n
                for ei = 1:sim.env.n
                    for oi = 0:sim.obj.n
                        si = find(sim.state.env == ei & sim.state.nexp == oi & rowis(sim.state.r(:,1:oi),sim.task.rewseq.r{rsi}(1:oi)));
                        assert(numel(si)==1,['could not find unique state in env ' num2str(ei) ' for reward sequence ' num2str(rsi)]);

                        if sim.model.besta(si,mi,epsi) == 1
                            sim.model.rewseq.nexp(ei,rsi,mi,epsi) = oi;
                            break;
                        end
                    end
                end
            end
        end
    end
end

% if requested, save the results to a file
if do_save
    disp(['saving ' savefilename]);
    save(savefilename,'sim');
end

% if requested, plot the results
if do_plot
    figuren;
    
    plot_model_names = {'Policy 1','Policy 2','Policy 3','Policy 4 (optimal)'};
    for epsi = 1:sim.env.param.n
        epset = sim.env.param.set{epsi};
        

        ncol_param = 2;
        
        nrow = sim.env.param.n;
        ncol = ncol_param + 1 + sim.model.n;

        rowstart = (epsi-1);
        
        color_rich = [245 127 32]./255;
        color_poor = [155 48 146]./255;
        
        nsubplot(nrow,ncol,rowstart+1,1);
        envcolor = arrayfun(@(ei) interpcolor(color_poor,color_rich,(ei-0.5)./sim.env.n),1:sim.env.n,'uniform',0);
        envtext = arrayfun(@(ei) [texcolor(envcolor{ei}) 'Env ' num2str(ei)],1:sim.env.n,'uniform',0);
        envtext = envtext(end:-1:1);
        
        curall = vertcat(epset.pbait_for_obj{:});
        if all(rowis(curall,curall(1,:)))
            plot(1:sim.obj.n,epset.pbait_for_obj{1},'.-','color','k','linewidth',2,'markersize',14);
        else
            for ei = 1:sim.env.n
                plot(1:sim.obj.n,epset.pbait_for_obj{ei},'.-','color',envcolor{ei},'linewidth',2,'markersize',14);
            end
        end
        set(gca,'xtick',1:sim.obj.n,'xlim',[1 sim.obj.n] + [-1 1]*.5);
        ylim([0 1]);
        if epsi == 1
            etextn('rt',envtext);
        end
        xlabel('# objects inspected');
        ylabel('p(reward)');
        if epsi == 1
            title({'Reward','probabilities'});
        end

        nsubplot(nrow,ncol,rowstart+1,2);
        curall = vertcat(epset.prew_given_bait{:});
        okx = sim.rew.all ~= 0;
        if all(rowis(curall,curall(1,:)))
            plot(sim.rew.all(okx),epset.prew_given_bait{1}(okx),'.-','color','k','linewidth',2,'markersize',14);
        else
            for ei = 1:sim.env.n
                plot(sim.rew.all(okx),epset.prew_given_bait{ei}(okx),'.-','color',envcolor{ei},'linewidth',2,'markersize',14);
            end
        end
        set(gca,'xtick',sim.rew.all(okx),'xlim',minmax(sim.rew.all(okx))' + [-1 1]*.5,'ylim',[0 1]);
        xlabel('Reward size');
        ylabel('p(reward size | reward)');
        if epsi == 1
            title({'Reward size','distributions'});
        end

        nsubplot(nrow,ncol,rowstart+1,ncol_param+1);
        plot(1:sim.model.n,sim.model.rewrate(1,:,epsi),'.-','color',interpcolor('k','w'),'linewidth',3,'markersize',20);
        setlim(gca,'ylim','tight',.1);
        
        yticks = minmax(sim.model.rewrate(1,:,epsi))';
        set(gca,'ytick',round(yticks,3),'yticklabel',arrayfun(@(z) sprintf('%.3f',z),round(yticks,3),'uniform',0));
        set(gca,'xtick',1:sim.model.n,'xticklabel',1:sim.model.n,'xlim',[1 sim.model.n] + [-1 1]*.4);
        xlabel('Policy');
        if epsi == 1
            title({'Reward rates','of policies'});
        end


        for mi = 1:sim.model.n
            hmod(2,mi)=nsubplot(nrow,ncol,rowstart+1,ncol_param+1+mi);
            for rsi = 1:sim.task.rewseq.n
                plot(sim.model.rewseq.nexp(:,rsi,mi,epsi),'o-','color',sim.task.rewseq.color{rsi},'linestyle',sim.task.rewseq.linestyle{rsi},'linewidth',2);
            end
            set(gca,'xtick',1:sim.env.n,'xlim',[1 sim.env.n] + [-1 1]*.5,'ylim',[0 sim.obj.n]);
            if epsi == 1
                etextn('lt',arrayfun(@(rsi) [texcolor(sim.task.rewseq.color{rsi}) sim.task.rewseq.name{rsi}],1:sim.task.rewseq.n,'uniform',0));
            end
            xlabel('Environment');
            if mi == 1
                ylabel('# objects explored');
            end
            set(gca,'ytick',1:2:9,'ylim',[1 9],'yticklabel',{'1','','','','9'});
            if epsi == 1
                title(plot_model_names{mi});
            end
        end
    end
    drawnow;
    scale(gcf,'hvscale',[1.75 3]);
    pause(0.05);
    scale('eachax',gcf,'hscale',1.2);
    
    
    %% plot decision boundary for specified simulation

    env_colors = arrayfun(@(ei) interpcolor(interpcolor('m','k',.5),interpcolor([1 .5 0],'k',.15),(ei-0.5)./sim.env.n),1:sim.env.n,'uniform',0);
    
    % Note: change this setting if you want to plot the optimal decision
    % boundaries from other environment parameter sets
    environment_param_sets_to_plot = 1;

    figuren;
    mi = sim.model.n;
    for epstpi = 1:numel(environment_param_sets_to_plot)
        epsi = environment_param_sets_to_plot(epstpi);

        rstate = sim.env.param.result{epsi}.rstate;
        state = sim.env.param.result{epsi}.state;

        bestcon = sim.model.besta(:,mi,epsi)==1;
        bestexp = sim.model.besta(:,mi,epsi)==2;

        vcon = sim.model.vcon(:,mi,epsi);

        policy_bound = nans(sim.env.n,sim.obj.n);
        policy_bound_accuracy = nans(sim.env.n,sim.obj.n);
        policy_bound_error = nans(sim.env.n,sim.obj.n);

        for ei = 1:sim.env.n
            for nexp = 1:sim.obj.n
                oks = state.nexp == nexp & state.env == ei;
                okscon = oks & bestcon;
                oksexp = oks & bestexp;
                if any(okscon) && any(oksexp)
                    vcon_min_con = min(vcon(okscon));
                    vcon_max_exp = max(vcon(oksexp));
                    vcon_bound = vcon_max_exp;
                    vcon_linestyle = '-';
                    if vcon_max_exp >= vcon_min_con
                        fprintf('warning, imperfect bound! model %d epsi %d ei %d nexp %d\n',mi,epsi,ei,nexp);
                        vcon_linestyle = ':';
                    end
                    policy_bound(ei,nexp) = vcon_bound;

                    predicted_a = nans(state.n,1);
                    optimal_a = sim.model.besta(:,mi,epsi);
                    predicted_a(okscon) = 1 + (vcon(okscon) <= policy_bound(ei,nexp));
                    predicted_a(oksexp) = 1 + (vcon(oksexp) <= policy_bound(ei,nexp));
                    pstate = sim.model.preach_state(:,mi,epsi);
                    policy_bound_accuracy(ei,nexp) = sum(pstate(oks) .* (predicted_a(oks) == optimal_a(oks))) ./ sum(pstate(oks));
                    policy_bound_error(ei,nexp) = sum(pstate(oks) .* (predicted_a(oks) ~= optimal_a(oks)));
                elseif ~any(okscon) && any(oksexp)
                    policy_bound(ei,nexp) = inf;
                    policy_bound_accuracy(ei,nexp) = 1;
                    policy_bound_error(ei,nexp) = 0;
                elseif any(okscon) && ~any(oksexp)
                    policy_bound(ei,nexp) = -inf;
                    policy_bound_accuracy(ei,nexp) = 1;
                    policy_bound_error(ei,nexp) = 0;
                end
            end
        end

        nsubplot(1,numel(environment_param_sets_to_plot),1,epstpi);
        curtext = {};
        for ei = 1:sim.env.n
            curx = 1:sim.obj.n;
            ok_to_plot = ~isinf(policy_bound(ei,:));
            plot(curx(ok_to_plot),policy_bound(ei,ok_to_plot),'.-','color',env_colors{ei},'linewidth',4,'markersize',32);
            curtext{end+1} = [texcolor(env_colors{ei}) 'Environment ' num2str(ei)];
        end
        curtext{end+1} = [];
        curtext{end+1} = [texcolor('k') 'The depicted policy'];
        curtext{end+1} = 'is accurate for';
        curtext{end+1} = [roundstr(.001,100 - 100*nansum(policy_bound_error(:))) '% of decisions'];

        etextn('rt',curtext,'fontsize',12);
        title(sprintf('Optimal policy decision boundaries\nfor environment parameter set\n"%s"',sim.env.param.set{epsi}.paramset_name));
        liney(sim.model.rewrate(1,mi,epsi),'k:');
        set(gca,'xtick',1:sim.obj.n);
        xlabel('# objects inspected');
        if epsi == 1
            ylabel({'Value of consuming revealed rewards','(R_{consume})','(reward/s))'});
        end

        set(gca,'fontsize',14);
        setlim(gca,'xlim','tight',.1);
        etext('rb',sim.obj.n,sim.model.rewrate(1,mi,epsi),{'avg rew rate'});
        axis square;
    end
    drawnow;
    setlim('all','ylim','tight',.1);
end
