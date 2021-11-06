clear all; %close all; clc;

addpath('functions'); %keep hepers functions here
addpath('MatlabToolbox-2.8');

do_load = 1;
do_calc = 1;
do_save = 1; savefilename = 'tranersearch_paper_script_v13.mat';
do_plot = 1;
NumberofPerm=10000;
load('Upload_storage_General.mat');
if do_calc
    % dense eye search parameters
    despar = struct();
    
    despar.bin.n_rc = [size(storage.denseEyeSearch{1},1) size(storage.denseEyeSearch{1},2)];
    despar.bin.center_rc = (despar.bin.n_rc./2) + 1;
    
    [despar.bin.col,despar.bin.row] = meshgrid(1:despar.bin.n_rc(2),1:despar.bin.n_rc(1));
    
    despar.bin.per_deg = (100/61); % based on paper Methods saying that a 61x61 deg space was divided into a 100x100 grid.
    
    despar.t.sample_dur_s = 1./1000; % one sample = 1 millisecond
    despar.t.gaze_dur_in_win_to_count_as_viewed_s = 50./1000; % must have gaze in window for at least 50 ms to count as viewing an object
    
    fractal_names = {'center','targ1','targ2','targ3','targ4','inn1','inn2','inn3','inn4','mid1','mid2','mid3','mid4','out1','out2','out3','out4', ...
        'everything_else'};
    
    % id, is center, is target, is environmental, environmental ring, angle (deg), eccentricity
    fractal_info_varnames = {'id','is_center','is_target','is_envir','is_ctrl_envir','is_everything_else','ring','angle_deg','ecc_deg'};
    fractal_info = [ ...
        1, 1, 0, 0, 0, 0, nan,     -0,  0; ...
        2, 0, 1, 0, 0, 0, nan,    -40, 12; ...
        3, 0, 1, 0, 0, 0, nan,   -130, 12; ...
        4, 0, 1, 0, 0, 0, nan,   -220, 12; ...
        5, 0, 1, 0, 0, 0, nan,   -310, 12; ...
        6, 0, 0, 1, 0, 0,   1,  -22.5,  6; ...
        7, 0, 0, 1, 0, 0,   1, -112.5,  6; ...
        8, 0, 0, 1, 0, 0,   1, -202.5,  6; ...
        9, 0, 0, 1, 0, 0,   1, -292.5,  6; ...
        10, 0, 0, 1, 0, 0,   2,  -67.5, 15; ...
        11, 0, 0, 1, 0, 0,   2, -157.5, 15; ...
        12, 0, 0, 1, 0, 0,   2, -247.5, 15; ...
        13, 0, 0, 1, 0, 0,   2, -337.5, 15; ...
        14, 0, 0, 1, 0, 0,   3,  -22.5, 17; ...
        15, 0, 0, 1, 0, 0,   3, -112.5, 17; ...
        16, 0, 0, 1, 0, 0,   3, -202.5, 17; ...
        17, 0, 0, 1, 0, 0,   3, -292.5, 17; ...
        18, 0, 0, 0, 0, 1, nan, nan, nan; ...
        ];
    
    fractal_info_varnames{end+1} = 'loc_x_deg';
    fractal_info_varnames{end+1} = 'loc_y_deg';
    fractal_info(:,end+1) = fractal_info(:,9).*cos(fractal_info(:,8).*(pi/180));
    fractal_info(:,end+1) = fractal_info(:,9).*sin(fractal_info(:,8).*(pi/180));
    
    fractal_info_varnames{end+1} = 'loc_r';
    fractal_info_varnames{end+1} = 'loc_c';
    fractal_info(:,end+1) = fractal_info(:,11).*despar.bin.per_deg + despar.bin.center_rc(1);
    fractal_info(:,end+1) = fractal_info(:,10).*despar.bin.per_deg + despar.bin.center_rc(2);
    
    tmp = fractal_info(:,end);
    fractal_info(:,end) = fractal_info(:,end-1);
    fractal_info(:,end-1) = tmp;
    
    despar.frac.info = array2table(fractal_info,'VariableNames',fractal_info_varnames);
    despar.frac.name = fractal_names;
    despar.frac.n = size(despar.frac.info,1);
    assert(despar.frac.n == numel(despar.frac.name));
    
    despar.frac.win.ok_rc = arrayfun(@(z) false(despar.bin.n_rc),(1:despar.frac.n)','uniform',0);
    despar.frac.win.nok_rc = nans(despar.frac.n,1);
    despar.frac.win.all.ok_rc = false(despar.bin.n_rc);
    despar.frac.win.all.nok_rc = nans(despar.frac.n,1);
    
    despar.frac.win.style = 'euclidean_match';
    
    % Manhattan distance (square(ish) windows)
    despar.frac.win.manhattan.max_dist_deg = 3;
    
    % Euclidean distance (circular(ish) windows)
    despar.frac.win.euclidean.max_dist_deg = 3;
    
    % Euclidean distance, matched bin
    % treat gaze as being on a fractal, if gaze is in a bin
    % that is one of THIS MANY closest bins to the center of
    % the fractal
    despar.frac.win.euclidean_match.n_bins = 77; % ~ 3 deg

    despar.frac.win.max_bin_dist_from_center_deg = nans(despar.frac.n,1);
    despar.frac.win.all_bin_dist_from_center_deg = cell(despar.frac.n,1);
    
    
    for fi = 1:despar.frac.n
        if strcmp(despar.frac.name{fi},'everything_else')
            assert(fi == despar.frac.n,'the "everything_else" object must be the last object');
            ok = ~despar.frac.win.all.ok_rc;
            despar.frac.win.ok_rc{fi} = ok;
            despar.frac.win.nok_rc(fi) = sum(ok(:));
            continue;
        end
        
        dist_from_center_bin_manhattan = max( ...
            abs(despar.bin.col - despar.frac.info{fi,'loc_c'}), ...
            abs(despar.bin.row - despar.frac.info{fi,'loc_r'}) );
        dist_from_center_bin_euclidean = sqrt( ...
            (despar.bin.col - despar.frac.info{fi,'loc_c'}).^2 + ...
            (despar.bin.row - despar.frac.info{fi,'loc_r'}).^2 );
        
        switch despar.frac.win.style
            case 'manhattan'
                dist_from_center_bin = dist_from_center_bin_manhattan;
            case {'euclidean','euclidean_match'}
                dist_from_center_bin = dist_from_center_bin_euclidean;
            otherwise
                error('unknown window style');
        end
        dist_from_center_deg = dist_from_center_bin ./ despar.bin.per_deg;
        
        ok = false(despar.bin.n_rc);
        switch despar.frac.win.style
            case {'manhattan','euclidean'}
                ok = dist_from_center_deg <= despar.frac.win.(despar.frac.win.style).max_dist_deg;
            case {'euclidean_match'}
                % sort all bins in the denseEyeSearch matrix based on their distance
                % from the center location of the fractal
                cur_bin_info = [despar.bin.row(:) despar.bin.col(:) dist_from_center_deg(:)];
                [sorted_bin_info,sid] = sortrows(cur_bin_info,+3);
                
                % pick the first n of those pixels and treat them as the window for
                % gazing at the fractal
                sorted_ok = sorted_bin_info(1:despar.frac.win.euclidean_match.n_bins,:);
                for okpixeli = 1:size(sorted_ok,1)
                    ok(sorted_ok(okpixeli,1),sorted_ok(okpixeli,2)) = true;
                end
            otherwise
                error('unknown window style');
        end
        
        okdist = dist_from_center_deg(ok);
        despar.frac.win.max_bin_dist_from_center_deg(fi) = max(okdist(:));
        despar.frac.win.all_bin_dist_from_center_deg{fi} = sort(okdist(:));

        despar.frac.win.ok_rc{fi} = ok;
        despar.frac.win.nok_rc(fi) = sum(ok(:));
        
        despar.frac.win.all.ok_rc = despar.frac.win.all.ok_rc | ok;
        despar.frac.win.all.nok_rc = sum(despar.frac.win.all.ok_rc(:));
    end
    
    
    
    des = struct();
    des.ntrials = numel(storage.trialnumber);
    des.search_dur_s = storage.searchDuration(:,1);
    
    des.rewid = nans(des.ntrials,1);
    des.envid = nans(des.ntrials,1);
    des.prew = nans(des.ntrials,1);
    
    des.rewid(storage.rewDur == 0) = 1;
    des.rewid(storage.rewDur > 0) = 2;
    
    des.envid(ismember(storage.envID,[1 7])) = 1;
    des.envid(ismember(storage.envID,[2 8])) = 2;
    des.envid(ismember(storage.envID,[3 9])) = 3;
    des.envid(ismember(storage.envID,[4 10])) = 4;
    des.envid(ismember(storage.envID,[5 11])) = 5;
    
    des.prew = (des.envid-1)*0.25;
    
    des.rew.name = {'NoR','Rew'};
    des.rew.delivered = [0 1];
    des.rew.n = numel(des.rew.name);
    
    des.env.name = {'0%','25%','50%','75%','100%'};
    des.env.prew = [0 .25 .50 .75 1.00];
    des.env.n = numel(des.env.name);
    
    des.monk.name = {'Animal M','Animal B'};
    des.monk.monkID = [1 2];
    des.monk.n = numel(des.monk.name);
    
    des.task.name = {'search','single fractal'};
    des.task.oktr = { ...
        ~isnan(des.envid)    & (storage.finType'==0.2 | storage.finType' >= 2.0), ...
        storage.envID' == -1 & (storage.finType'==0.3 | storage.finType' >= 2.0)};
    des.task.n = numel(des.task.name);
    
    des.correct.name = {'search fix1','search','single fractal'};
    des.correct.finType = {0.1,0.2,0.3};
    des.correct.n = numel(des.correct.name);
    
    des.error.name = {'fix1 timeout','fix1 break','fix2 timeout','fix2 break','fix2 other','search timeout'};
    des.error.finType = {1.0,1.1,2.0,2.1,2.2,3};
    des.error.n = numel(des.error.name);
    
    des.monkid = nans(des.ntrials,1);
    des.taskid = nans(des.ntrials,1);
    des.correctid = nans(des.ntrials,1);
    des.errorid = nans(des.ntrials,1);
    
    for mi = 1:des.monk.n
        des.monkid(storage.monkID == des.monk.monkID(mi)) = mi;
    end
    for ti = 1:des.task.n
        des.taskid(des.task.oktr{ti}) = ti;
    end
    for cori = 1:des.correct.n
        des.correctid(ismember(storage.finType(:),des.correct.finType{cori})) = cori;
    end
    for erri = 1:des.error.n
        des.errorid(ismember(storage.finType(:),des.error.finType{erri})) = erri;
    end
    
    % correct and error trials
    des.cor = ~isnan(des.correctid);
    des.err = ~isnan(des.errorid);
    assert(~any(des.cor(des.err)) && ~any(des.err(des.cor)));
    
    % which fractal ID is the target on each trial
    despar.frac.all_targ_fractal_ids = find(despar.frac.info{:,'is_target'});
    des.targ_fractal_id = nans(numel(storage.targDir),1);
    for fid = 1:numel(despar.frac.all_targ_fractal_ids)
        fi = despar.frac.all_targ_fractal_ids(fid);
        curangle = despar.frac.info{fi,'angle_deg'};
        while abs(curangle) > 360
            curangle = curangle - sign(curangle)*360;
        end
        curok = abs(storage.targDir - curangle) < 0.0001;
        curok = curok | abs(storage.targDir - (360 + curangle)) < 0.0001;
        des.targ_fractal_id(curok) = fi;
    end
    
    % OK trials for analysis
    % must have a valid outcome (reward or no reward)
    des.oktr = ~isnan(des.rewid);
    % if is search task, must have an environment ID
    des.oktr = des.oktr & (des.taskid~=1 | (~isnan(des.envid) & storage.envFracAnalGood));
    % must not be a control trial with novel fractal
    des.oktr = des.oktr & storage.novelTrial(:)~=1;
    
    % for each trial, calculate:
    % - total number of samples anywhere in denseEyeSearch matrix
    % - total number of samples in ANY window
    % - number of data samples in EACH fractal's window
    % - total search duration (in sec)
    des.n_samples = nans(des.ntrials,1);
    des.n_in_any_win = nans(des.ntrials,1);
    des.n_in_win = nans(des.ntrials,despar.frac.n);
    
    curtimer = tic;
    last_timer_printout = toc(curtimer);
    fprintf('starting per-trial quantification of gaze in windows\n');
    for tr = 1:des.ntrials
        if toc(curtimer) > last_timer_printout + 30
            last_timer_printout = toc(curtimer);
            fprintf(' trial %d/%d, %d sec\n',tr,des.ntrials,round(toc(curtimer)));
        end
        
        cur = storage.denseEyeSearch{tr};
        des.n_samples(tr) = sum(cur,'all');
        
        tmp = cur .* despar.frac.win.all.ok_rc;
        des.n_in_any_win(tr) = sum(tmp,'all');
        
        for fi = 1:despar.frac.n
            tmp = cur .* despar.frac.win.ok_rc{fi};
            des.n_in_win(tr,fi) = sum(tmp,'all');
        end
    end
    fprintf(' ...done!\n');
    
    % viewing of center, target, and
    % other possible-but-not-shown-this-trial target locations
    des.n_in_envir_wins = sum(des.n_in_win(:,despar.frac.info{:,'is_envir'}==1),2);
    des.n_in_ctrle_wins = sum(des.n_in_win(:,despar.frac.info{:,'is_ctrl_envir'}==1),2);
    des.n_in_all_targ_wins = sum(des.n_in_win(:,despar.frac.info{:,'is_target'}==1),2);
    des.n_in_center_win = sum(des.n_in_win(:,despar.frac.info{:,'is_center'}==1),2);
    des.n_in_targ_win = nans(size(des.n_in_center_win));
    des.n_in_nontarg_wins = nans(size(des.n_in_center_win));
    des.n_in_else_win = sum(des.n_in_win(:,despar.frac.info{:,'is_everything_else'}==1),2);
    
    % number of objects of each type that were viewed this trial
    min_samples_to_count_as_viewed = round(despar.t.gaze_dur_in_win_to_count_as_viewed_s .* 1000);
    des.n_envir_frac_viewed = sum(des.n_in_win(:,despar.frac.info{:,'is_envir'}==1)>=min_samples_to_count_as_viewed,2);
    des.n_ctrle_frac_viewed = sum(des.n_in_win(:,despar.frac.info{:,'is_ctrl_envir'}==1)>=min_samples_to_count_as_viewed,2);
    des.was_center_viewed = sum(des.n_in_win(:,despar.frac.info{:,'is_center'}==1)>=min_samples_to_count_as_viewed,2);
    des.was_targ_viewed = false(size(des.n_in_center_win));
    des.n_nontarg_viewed = false(size(des.n_in_center_win));
    des.was_else_viewed = sum(des.n_in_win(:,despar.frac.info{:,'is_everything_else'}==1)>=min_samples_to_count_as_viewed,2);
    
    for fid = 1:numel(despar.frac.all_targ_fractal_ids)
        fi = despar.frac.all_targ_fractal_ids(fid);
        ok = des.targ_fractal_id == fi;
        
        otherfi = setdiff(despar.frac.all_targ_fractal_ids,fi);
        
        des.n_in_targ_win(ok) = sum(des.n_in_win(ok,fi),2);
        des.n_in_nontarg_wins(ok) = sum(des.n_in_win(ok,otherfi),2);
        
        des.was_targ_viewed(ok) = sum(des.n_in_win(ok,fi)>=min_samples_to_count_as_viewed,2);
        des.n_nontarg_viewed(ok) = sum(des.n_in_win(ok,otherfi)>=min_samples_to_count_as_viewed,2);
    end
    
    % BASIC STATISTICS FOR EACH MONK
    monkstats = {};
    monkstats{end+1} = struct('name','all','mi',3,'oktr',true(size(des.monkid)));
    monkstats{end+1} = struct('name',des.monk.name{2},'mi',2,'oktr',des.monkid==2);
    monkstats{end+1} = struct('name',des.monk.name{1},'mi',1,'oktr',des.monkid==1);
    
    % REWARD RATE CALCULATIONS
    % which error types will we analyze?
    cor_okerrids = find(ismember(des.error.name,{}));
    cor_okcorids = find(ismember(des.correct.name,{'search','single fractal'}));
    err_okerrids = find(ismember(des.error.name,{'search timeout'}));
    err_okcorids = find(ismember(des.correct.name,{}));
    
    oktr_cor = (ismember(des.errorid,cor_okerrids) | ismember(des.correctid,cor_okcorids)) & ~isnan(des.search_dur_s);
    oktr_err = ismember(des.errorid,err_okerrids) | ismember(des.correctid,err_okcorids);
    
    monkstats_task = struct();
    monkstats_task.name = {'1fract','search','searchC','searchU'};
    monkstats_task.oktr = {des.taskid==2,des.taskid==1,des.taskid==1 & ismember(des.prew,[0 1]),des.taskid==1 & ismember(des.prew,[0.25 0.50 0.75])};
    monkstats_task.n = numel(monkstats_task.name);
    
    for pmi = 1:numel(monkstats)
        mi = monkstats{pmi}.mi;
        oktr_monk = monkstats{pmi}.oktr;
        monkname = monkstats{pmi}.name;
        
        taskstat = struct();
        taskstat.monk.monkid = mi;
        taskstat.monk.oktr = oktr_monk;
        taskstat.monk.name = monkname;
        
        % col 1: stats for norew trials
        % col 2: stats for rew trials
        % col 3: stats for all trials
        % col 4: stats for what all trials WOULD be if p(error) and
        %        E[search dur | correct] on norew trials were the same
        %        as rew trials
        taskstat.col.name = {'norew','rew','all','all_theory'};
        taskstat.col.n = numel(taskstat.col.name);
        
        % comparison of actual reward rate vs. theoretical reward rate
        % that would be achieved if always behaved 'as if' the target
        % was a rewarded target
        taskstat.comp.name = {'ratio','diff'};
        taskstat.comp.n = numel(taskstat.comp.name);
        
        taskstat.boot.n = 2000;
        
        taskstat.n= nans(des.task.n,taskstat.col.n);
        
        taskstat.perr = nans(des.task.n,taskstat.col.n);
        taskstat.erew_cor = nans(des.task.n,taskstat.col.n);
        taskstat.cor_dur_s = nans(des.task.n,taskstat.col.n);
        taskstat.edur_s = nans(des.task.n,taskstat.col.n);
        taskstat.rew_per_s = nans(des.task.n,taskstat.col.n);
        taskstat.rew_vs_theory = nans(des.task.n,taskstat.comp.n);
        
        taskstat.boot.perr = nans(des.task.n,taskstat.col.n,taskstat.boot.n);
        taskstat.boot.erew_cor = nans(des.task.n,taskstat.col.n,taskstat.boot.n);
        taskstat.boot.cor_dur_s = nans(des.task.n,taskstat.col.n,taskstat.boot.n);
        taskstat.boot.edur_s = nans(des.task.n,taskstat.col.n,taskstat.boot.n);
        taskstat.boot.rew_per_s = nans(des.task.n,taskstat.col.n,taskstat.boot.n);
        taskstat.boot.rew_vs_theory = nans(des.task.n,taskstat.comp.n,taskstat.boot.n);
        
        taskstat.boot.se.perr = nans(des.task.n,taskstat.col.n);
        taskstat.boot.se.erew_cor = nans(des.task.n,taskstat.col.n);
        taskstat.boot.se.cor_dur_s = nans(des.task.n,taskstat.col.n);
        taskstat.boot.se.edur_s = nans(des.task.n,taskstat.col.n);
        taskstat.boot.se.rew_per_s = nans(des.task.n,taskstat.col.n);
        taskstat.boot.se.rew_vs_theory = nans(des.task.n,taskstat.comp.n);
        
        for pti = 1:monkstats_task.n
            for booti = 0:taskstat.boot.n
                oktr = oktr_monk & (oktr_cor | oktr_err) & monkstats_task.oktr{pti} & ~isnan(des.rewid) & (~isnan(des.search_dur_s) | des.err);
                
                % which OK trials are from norew vs rew trials
                assert(strcmp(des.rew.name{1},'NoR') && strcmp(des.rew.name{2},'Rew'));
                oktr_nor = oktr & des.rewid == 1;
                oktr_rew = oktr & des.rewid == 2;
                assert(all(oktr == (oktr_nor | oktr_rew)));
                
                % do in terms of trial IDs instead of a logical vector,
                % so can resample with replacement
                oktrid_nor = find(oktr_nor);
                oktrid_rew = find(oktr_rew);
                oktrid = [oktrid_nor ; oktrid_rew];
                
                % if is a bootstrap, resample trials with replacement
                % separately for nonrewarded trials and rewarded trials
                if booti > 0
                    oktrid_nor = oktrid_nor(randsample(numel(oktrid_nor),numel(oktrid_nor),true));
                    oktrid_rew = oktrid_rew(randsample(numel(oktrid_rew),numel(oktrid_rew),true));
                    oktrid = [oktrid_nor ; oktrid_rew];
                end
                
                % which of them are correct vs. error trials
                oktrid_cor = oktrid(des.cor(oktrid));
                oktrid_err = oktrid(des.err(oktrid));
                oktrid_nor_cor = oktrid_nor(des.cor(oktrid_nor));
                oktrid_nor_err = oktrid_nor(des.err(oktrid_nor));
                oktrid_rew_cor = oktrid_rew(des.cor(oktrid_rew));
                oktrid_rew_err = oktrid_rew(des.err(oktrid_rew));
                
                if booti == 0
                    fprintf('pmi %d, pti %d\n %5d %5d %5d\n %5d %5d %5d\n %5d %5d\n %5d %5d\n %5d %5d\n\n', ...
                        pmi, pti, ...
                        sum(oktr), sum(oktr_nor), sum(oktr_rew),  ...
                        numel(oktrid),  numel(oktrid_nor),  numel(oktrid_rew),  ...
                        numel(oktrid_cor),  numel(oktrid_err),  ...
                        numel(oktrid_nor_cor),  numel(oktrid_nor_err),  ...
                        numel(oktrid_rew_cor),  numel(oktrid_rew_err));
                end
                
                assert( ...
                    all(des.errorid([oktrid_nor_err ; oktrid_rew_err]) == find(strcmp(des.error.name,'search timeout'))) ...
                    && all(isnan(des.errorid([oktrid_nor_cor ; oktrid_rew_cor]))), ...
                    'expected all errors to be timeouts, and all correct trials to be successfully completed');
                
                ntrials = [numel(oktrid_nor) numel(oktrid_rew) numel(oktrid) numel(oktrid)];
                
                prew_if_correct = numel(oktrid_rew_cor)/(numel(oktrid_rew_cor) + numel(oktrid_nor_cor));
                erew_cor = [0 1 prew_if_correct prew_if_correct];
                
                perr = [mean(des.err(oktrid_nor)) mean(des.err(oktrid_rew)) mean(des.err(oktrid)) mean(des.err(oktrid_rew))];
                pcor = 1 - perr;
                
                cor_dur_s = [mean(des.search_dur_s(oktrid_nor_cor)) mean(des.search_dur_s(oktrid_rew_cor)) mean(des.search_dur_s(oktrid_cor)) mean(des.search_dur_s(oktrid_rew_cor))];
                
                % add in other time periods after search starts
                % and before the trial ends:
                % - required time to hold fixation on the target
                % - required time to wait after fixating the target before the outcome
                % - inter-trial interval
                search_fixation_requirement_s = 0.75;
                pre_outcome_period_s = 0.5;
                mean_iti_duration_s = 1.5; %1-2 sec, per Methods
                
                cor_dur_s = cor_dur_s + search_fixation_requirement_s + pre_outcome_period_s + mean_iti_duration_s;
                err_dur_s = 5 + mean_iti_duration_s;
                
                % expected total search duration to complete a
                % single trial, including time taken during error
                % trials before the final completion of the trial
                % (num errors) = pcor*0 + perr*pcor*1 + perr^2*pcor*2 + ... + perr^N*pcor*N
                nerr = (0:100)';
                enum_errors = sum((nerr*pcor).*(perr.^nerr),1);
                edur_s = cor_dur_s + enum_errors .* err_dur_s;
                
                rew_per_s = erew_cor ./ edur_s;
                
                % calculate directly based on separate reward rates and
                % durations for reward trials vs. noreward trials
                prew = prew_if_correct;
                
                % alternately, assume scheduled prew is truly 50/50. In 
                % actuality, seems slightly below 50/50 because animal may 
                % end the session more often by refusing to complete a 
                % no-reward trial
                %prew = 0.5;
                
                rew_per_s(3) = ((1-prew)*erew_cor(1) + (prew)*erew_cor(2)) ./ ((1-prew)*edur_s(1) + (prew)*edur_s(2));
                rew_per_s(4) = ((1-prew)*erew_cor(1) + (prew)*erew_cor(2)) ./ ((1-prew)*edur_s(2) + (prew)*edur_s(2)); % effect of error rates and durations being 'as if' the trial was always rewarded
                
                rew_diff_vs_theory = rew_per_s(3) - rew_per_s(4);
                rew_ratio_vs_theory = rew_per_s(3) ./ rew_per_s(4);
                rew_vs_theory = [rew_ratio_vs_theory rew_diff_vs_theory];
                
                if booti == 0
                    taskstat.n(pti,:) = ntrials;
                    taskstat.perr(pti,:) = perr;
                    taskstat.erew_cor(pti,:) = erew_cor;
                    taskstat.cor_dur_s(pti,:) = cor_dur_s;
                    taskstat.edur_s(pti,:) = edur_s;
                    taskstat.rew_per_s(pti,:) = rew_per_s;
                    taskstat.rew_vs_theory(pti,:) = rew_vs_theory;
                else
                    taskstat.boot.perr(pti,:,booti) = perr;
                    taskstat.boot.erew_cor(pti,:,booti) = erew_cor;
                    taskstat.boot.cor_dur_s(pti,:,booti) = cor_dur_s;
                    taskstat.boot.edur_s(pti,:,booti) = edur_s;
                    taskstat.boot.rew_per_s(pti,:,booti) = rew_per_s;
                    taskstat.boot.rew_vs_theory(pti,:,booti) = rew_vs_theory;
                end
            end
            % calc bootstrap summary stats
            taskstat.boot.se.perr(pti,:) = std(taskstat.boot.perr(pti,:,:),[],3);
            taskstat.boot.se.erew_cor(pti,:) = std(taskstat.boot.erew_cor(pti,:,:),[],3);
            taskstat.boot.se.cor_dur_s(pti,:) = std(taskstat.boot.cor_dur_s(pti,:,:),[],3);
            taskstat.boot.se.edur_s(pti,:) = std(taskstat.boot.edur_s(pti,:,:),[],3);
            taskstat.boot.se.rew_per_s(pti,:) = std(taskstat.boot.rew_per_s(pti,:,:),[],3);
            taskstat.boot.se.rew_vs_theory(pti,:) = std(taskstat.boot.rew_vs_theory(pti,:,:),[],3);
        end
        
        % comparison of "diff between actual vs. theoretical rewrate"
        % between different task conditions
        % row 1: search vs 1fract
        % row 2: searchU vs searchC
        taskstat.boot.rew_vs_theory_diff_ci95 = cell(2,taskstat.comp.n);
        taskstat.boot.rew_vs_theory_diff_pseudo_pvalue = cell(2,taskstat.comp.n);
        cur = taskstat.boot.rew_vs_theory;
        for ci = 1:taskstat.comp.n
            curdiff = cur(2,ci,:) - cur(1,ci,:);
            taskstat.boot.rew_vs_theory_diff_ci95{1,ci} = quantile(curdiff,[.025 .975]);
            taskstat.boot.rew_vs_theory_diff_pseudo_pvalue{1,ci} = 2*min([1-mean(curdiff < 0) 1-mean(curdiff > 0)]);
            
            curdiff = cur(4,ci,:) - cur(3,ci,:);
            taskstat.boot.rew_vs_theory_diff_ci95{2,ci} = quantile(curdiff,[.025 .975]);
            taskstat.boot.rew_vs_theory_diff_pseudo_pvalue{2,ci} = 2*min([1-mean(curdiff < 0) 1-mean(curdiff > 0)]);
        end
        
        monkstats{pmi}.rewrate = taskstat;
    end
end

if do_plot
    
    rewcolor = {[157 155 204]./255,[243 118 132]./255};
    
    if 1
        % row 1, col 1: log(gaze), 50% rew, R targ
        % row 1, col 2: log(gaze), 50% rew, N targ
        % row 1, col 3-5: # fractals explored
        % row 2, cols 1-5: dur total, env, targ, start pos, all else
        % row 3, cols 1-3: % RR var expl, beta r, beta p(r)
        
        plot_taskid = 1;
        plot_monkid = 3;
        plot_oktr_all = des.oktr & des.cor & des.taskid == plot_taskid;
        if ismember(plot_monkid,1:des.monk.n)
            plot_oktr_all = plot_oktr_all & des.monkid == plot_monkid;
            plot_monkname = des.monk.name{plot_monkid};
        else
            plot_monkname = 'all';
        end
        
        
        
        %% HEATMAP PLOTS
        allhmp = {};
        figuren;
        allhmp{end+1} = struct('rcargs',{{4,2,1,1}},'monkid',2,'taskid',2,'rewid',2,'targ_fractal_id',5);
        allhmp{end+1} = struct('rcargs',{{4,2,1,2}},'monkid',2,'taskid',2,'rewid',1,'targ_fractal_id',5);
        
        for hmpi = 1:numel(allhmp)
            hmp = allhmp{hmpi};
            
            oktypes = {'monkid','taskid','rewid','targ_fractal_id'};
            oktr_all = true(numel(des.n_samples),numel(oktypes));
            okname = '';
            for okti = 1:numel(oktypes)
                on = oktypes{okti};
                on_short = on(1:(end-2));
                oktr_all(:,okti) = oktr_all(:,okti) & ~isnan(des.(on));
                if ~isnan(hmp.(on))
                    oktr_all(:,okti) = oktr_all(:,okti) & des.(on) == hmp.(on);
                    if isfield(des,on_short)
                        okname = [okname ' ' des.(on_short).name{hmp.(on)}];
                    else
                        okname = [okname ' ' on ' = ' num2str(hmp.(on))];
                    end
                else
                    okname = [okname ' ' on ' = all'];
                end
            end
            oktr = all(oktr_all,2);
            oktr = oktr & des.cor;
            
            oktrid = find(oktr);

            gazemap = zeros(size(storage.denseEyeSearch{1}));
            for tri = 1:numel(oktrid)
                tr = oktrid(tri);
                gazemap = gazemap + storage.denseEyeSearch{tr};
            end
            
            % convert from ms to ms/trial
            gazemap = gazemap ./ sum(oktr);

            % convert to log2 units
            gazemap = log2(gazemap);
            gazemap(isinf(gazemap)) = nan;
            
            % swap rows and cols
            gazemap = gazemap';
            
            allhmp{hmpi}.oktr = oktr;
            allhmp{hmpi}.okname = okname;
            allhmp{hmpi}.gazemap = gazemap;
        end
        
        % get common color limits for all plots
        allgazemap = cellfun(@(z) z.gazemap,allhmp,'uniform',0);
        allgazemap = horzcat(allgazemap{:});
        allgazemap = allgazemap(:);
        allgazemap = allgazemap(~isnan(allgazemap) & ~isinf(allgazemap));
        cl = [min(allgazemap(:)) max(allgazemap(:))];
        
        
        for hmpi = 1:numel(allhmp)
            hmp = allhmp{hmpi};
            
            arg = hmp.rcargs;
            allhmp{hmpi}.h = nsubplot(arg{1},arg{2},arg{3},arg{4});
            
            gazemap = allhmp{hmpi}.gazemap;
            ntrials = sum(allhmp{hmpi}.oktr);
            
            x = ((1:despar.bin.n_rc(2))-despar.bin.center_rc(2));
            raw_bin_xl = x([1 end]) + [-1 1]*0.5;
            
            y = ((1:despar.bin.n_rc(1))-despar.bin.center_rc(1));
            raw_bin_yl = y([1 end]) + [-1 1]*0.5;
            
            xl_deg = [-20 20] + [-1 1]*0.5;
            xl_bin = xl_deg * despar.bin.per_deg;

            yl_deg = [-20 20] + [-1 1]*0.5;
            yl_bin = xl_deg * despar.bin.per_deg;
            
            xtick_deg = -20:20:20;
            ytick_deg = -20:20:20;
            xtick_bin = xtick_deg .* despar.bin.per_deg;
            ytick_bin = ytick_deg .* despar.bin.per_deg;
            set(gca,'xtick',xtick_bin,'xticklabel',xtick_deg,'ytick',ytick_bin,'yticklabel',ytick_deg);

            xlabel('Horizontal position (deg)');
            xlim(xl_bin);
            ylabel('Vertical position (deg)');
            ylim(yl_bin);
            
            image(x,y,colormapify(gazemap,cl,'k','r','y','k'));
            axis square;
            
            curtext = {};
            
            frac_colors = {};
            fracloc_id = hmp.targ_fractal_id;
            if isnan(fracloc_id)
                fracloc_id = find(despar.frac.info{:,'is_target'});
            end
            if isnan(hmp.rewid)
                frac_colors{end+1} = struct('name','Target','id',fracloc_id,'color',interpcolor('m','k',.25),'style','-');
            elseif des.rew.delivered(hmp.rewid) == 0
                frac_colors{end+1} = struct('name','No reward target','id',fracloc_id,'color',rewcolor{1},'style','-');
            elseif des.rew.delivered(hmp.rewid) == 1
                frac_colors{end+1} = struct('name','Reward target','id',fracloc_id,'color',rewcolor{2},'style','-');
            else
                error('unknown target');
            end
            frac_colors{end+1} = struct('name','Starting pos','id',find(despar.frac.info{:,'is_center'}),'color',[1 1 1]*.5,'style','-');
            
            for fci = 1:numel(frac_colors)
                curname = frac_colors{fci}.name;
                curcolor = frac_colors{fci}.color;
                curstyle = frac_colors{fci}.style;
                curtext{end+1} = [texcolor(curcolor) curname];
                
                for curfi = 1:numel(frac_colors{fci}.id)
                    fi = frac_colors{fci}.id(curfi);
                    [currow,curcol] = find(despar.frac.win.ok_rc{fi});

                    currow = currow - despar.bin.center_rc(2);
                    curcol = curcol - despar.bin.center_rc(1);
                    k = convhull(curcol,currow);
                    if contains(curname,'target')
                        % swap rows and cols, since data is in (x,y) format
                        plot(currow(k),curcol(k),'linestyle',curstyle,'color','w','linewidth',4);
                    end
                    % swap rows and cols, since data is in (x,y) format
                    plot(currow(k),curcol(k),'linestyle',curstyle,'color',curcolor,'linewidth',2);
                end
            end
            etextn('lt',curtext);
        end
        
        pos = get(gca,'pos');
        hscalebar = axes('position',[(pos(1)+pos(3)+.005) pos(2) pos(3)*.05 pos(4)]);
        sb_z = linspace(cl(1),cl(2),101)';
        sb_y = 1:numel(sb_z);
        sb_yl = [min(sb_y) max(sb_y)] + [-1 +1]*.5;
        
        
        sb_ztick = -10:2:10;
        sb_zticklab = 2.^sb_ztick;
        sb_ytick = min(sb_y) + (max(sb_y) - min(sb_y))*((sb_ztick - cl(1)) ./ (cl(2) - cl(1)));
        sb_yticklab = cell(size(sb_zticklab));
        for yi = 1:numel(sb_yticklab)
            cur = sb_zticklab(yi);
            assert(cur > 0);
            if cur >= 0.99
                sb_yticklab{yi} = sprintf('%d',round(cur));
            else
                sb_yticklab{yi} = sprintf('1/%d',round(1./cur));
            end
        end
        
        oky = inbounds(sb_ytick,[min(sb_y) max(sb_y)]);
        sb_ytick = sb_ytick(oky);
        sb_yticklab = sb_yticklab(oky);
        
        [sb_ytick,sid] = sort(sb_ytick);
        sb_yticklab = sb_yticklab(sid);
        
        imagesc(colormapify(sb_z,cl,'k','r','y','k'));
        ylim(sb_yl);
        xlim([.5 1.5]);
        set(hscalebar,'box','off','tickdir','out','ydir','normal','xtick',[],'ylim',sb_yl,'ytick',sb_ytick,'yticklabel',sb_yticklab,'yaxislocation','right');
        
        axis_scale_factor = [1.2 1.2];
        for hmpi = 1:numel(allhmp)
            scale(allhmp{hmpi}.h,'cm','hvscale',axis_scale_factor);
        end
        scale(hscalebar,'cm','hvscale',axis_scale_factor);
        
        
        
    end
end