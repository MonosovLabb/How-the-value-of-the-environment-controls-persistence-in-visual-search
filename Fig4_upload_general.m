clear all;
rng(1);
addpath('functions'); %keep helper functions here
addpath('MatlabToolbox-2.8');

do_load = 1;
do_calc = 1;
do_save = 1; savefilename = 'tranersearch_paper_script_v13.mat';
do_plot = 1;
NumberofPerm=10000;

tic
load('Upload_storage_General.mat');
toc
if do_calc
    % dense eye search parameters
    despar = struct();
    
    despar.bin.n_rc = [size(storage.denseEyeSearch{1},1) size(storage.denseEyeSearch{1},2)];
    despar.bin.center_rc = (despar.bin.n_rc./2) + 1;
    
    [despar.bin.col,despar.bin.row] = meshgrid(1:despar.bin.n_rc(2),1:despar.bin.n_rc(1));
    
    despar.bin.per_deg = (100/61); 

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
end

if do_plot
    
    rewcolor = {[157 155 204]./255,[243 118 132]./255};

    if 1
        % row 1, col 1: log(gaze), 50% rew, R targ
        % row 1, col 2: log(gaze), 50% rew, N targ
        % row 1, col 3-5: # fractals explored
        % row 2, cols 1-5: dur total, env, targ, start pos, all else
        % row 3, cols 1-3: % RR var expl, beta r, beta p(r)
        figuren;
        scale(gcf,'hvscale',[1 2.2]);
        
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
        allhmp{end+1} = struct('rcargs',{{4,2,1,1}},'monkid',2,'taskid',1,'envid',3,'rewid',2,'targ_fractal_id',5);
        allhmp{end+1} = struct('rcargs',{{4,2,1,2}},'monkid',2,'taskid',1,'envid',3,'rewid',1,'targ_fractal_id',5);
        
        for hmpi = 1:numel(allhmp)
            hmp = allhmp{hmpi};
            
            oktypes = {'monkid','taskid','rewid','envid','targ_fractal_id'};
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

            title(hmp.okname,'interpreter','none');
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
            frac_colors{end+1} = struct('name','Environment','id',find(despar.frac.info{:,'is_envir'}),'color',[1 1 1]*.9,'style','-');
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
        

        
        %% NUM ENVIR EXPLORED PLOT
        nsubplot(4,1,2,1);
        ylabel({'# Environmental fractals','inspected'});
        xlabel('Environment reward probability');
        set(gca,'xtick',1:des.env.n,'xticklabel',des.env.name);
        curtext = {};
        for ri = 1:des.rew.n
            barwidth = 0.4;
            xoff = (ri-1.5)*barwidth;
            xall = [];
            yall = [];
            xymean = [];
            hline = plot(1,1,'color',rewcolor{ri},'linewidth',2);
            for ei = 1:des.env.n
                oktr = plot_oktr_all & des.rewid==ri & des.envid==ei;
                x = ei + xoff;
                y = mean(des.n_envir_frac_viewed(oktr));
                yse = sem(des.n_envir_frac_viewed(oktr));
                plot_errorbar(x,y,yse,{'bar',rewcolor{ri},.4,'edgecolor','k'},{'color','k','linewidth',2});
                
                xymean = [xymean ; [x y]];
                xall = [xall ; des.prew(oktr)];
                yall = [yall ; des.n_envir_frac_viewed(oktr)];
            end
            set(hline,'xdata',xymean(:,1),'ydata',xymean(:,2));
            [p, rho] = permutation_pair_test_fast(xall,yall,NumberofPerm,'rankcorr');
            curtext{end+1} = [texcolor(rewcolor{ri}) sprintf('r = %.3f\np = %.3f',rho,p)];
        end
        etextn('lt',curtext);
        setlim(gca,'ylim','tight');
        yl = ylim(gca);
        yl(1) = 0;
        yl(2) = ceil(yl(2)*2)/2;
        ylim(yl);
        
        %% SEARCH DUR PORTIONS PLOT
        
        plotvar = {};
        plotvar{end+1} = struct('name','search_dur_s','y',des.search_dur_s);
        plotvar{end+1} = struct('name','env_gaze_dur_s','y',des.n_in_envir_wins ./ 1000);
        plotvar{end+1} = struct('name','tar_gaze_dur_s','y',des.n_in_targ_win ./ 1000);
        plotvar{end+1} = struct('name','fix_gaze_dur_s','y',des.n_in_center_win ./ 1000);
        plotvar{end+1} = struct('name','all_other_time_dur_s','y',clamp(des.search_dur_s - ((des.n_in_envir_wins + des.n_in_targ_win + des.n_in_center_win) ./ 1000),0,des.search_dur_s));

        plotmonk = {};
        switch plot_monkid
            case {1,2}
                plotmonk{end+1} = struct('name',des.monk.name{plot_monkid},'oktr',des.monkid==plot_monkid);
            case 3
                plotmonk{end+1} = struct('name','all','oktr',des.monkid==1 | des.monkid==2);
        end
        
        plotstats = cell(numel(plotmonk),numel(plotvar));
        
        assert(numel(plotmonk)==1);
        
        h = [];
        for mi = 1:numel(plotmonk)
            for pvi = 1:numel(plotvar)
                pv = plotvar{pvi};

                y_env_mean = nans(des.env.n,des.rew.n);
                y_env_se = nans(des.env.n,des.rew.n);

                raw_xy = cell(1,des.rew.n);
                raw_search_dur_s = cell(1,des.rew.n);

                for ei = 1:des.env.n
                    for ri = 1:des.rew.n
                        oktr = des.oktr & des.cor & des.taskid == plot_taskid & plotmonk{mi}.oktr & des.envid == ei & des.rewid == ri & ~isnan(pv.y);
                        if ~any(oktr)
                            continue;
                        end

                        y = pv.y(oktr);
                        y_env_mean(ei,ri) = mean(y);
                        y_env_se(ei,ri) = sem(y);

                        raw_xy{ri} = [raw_xy{ri} ; [des.env.prew(ei)*ones(numel(y),1) y]];
                        raw_search_dur_s{ri} = [raw_search_dur_s{ri} ; des.search_dur_s(oktr)];
                    end
                end
                

                h(mi,pvi)=nsubplot(4,numel(plotvar),3,pvi);

                curtext = {};
                for ri = 1:des.rew.n
                    [p, rho] = permutation_pair_test_fast(raw_xy{ri}(:,1),raw_xy{ri}(:,2),NumberofPerm,'rankcorr');

                    curcolor = rewcolor{ri};
                    plot_errorbar(des.env.prew,y_env_mean(:,ri),y_env_se(:,ri),'color',curcolor,'linewidth',3);
                end
                
                glm_x = [ ...
                    des.rew.delivered(1)*ones(size(raw_xy{1},1),1) raw_xy{1}(:,1) ; ...
                    des.rew.delivered(2)*ones(size(raw_xy{2},1),1) raw_xy{2}(:,1) ; ...
                    ];

                glm_y = [ ...
                    raw_xy{1}(:,2) ; ...
                    raw_xy{2}(:,2) ; ...
                    ];
                search_dur_y = [ ...
                    raw_search_dur_s{1} ; ...
                    raw_search_dur_s{2} ; ...
                    ];
                
                [b,dev,stats] = glmfit(glm_x,glm_y,'normal','link','identity');
                
                bootstor{mi, pvi} = bootstrp(2000, @(x) {glmfit(glm_x,glm_y,'normal','link','identity')}, [glm_x,glm_y]);
                
                glm_y_pred = glmval(b,glm_x,'identity');

                [urows,uia,uic] = unique(glm_x,'rows');
                full_y_pred = nans(size(glm_y_pred));
                full_sdur_pred = nans(size(glm_y_pred));
                for ui = 1:size(urows,1)
                    okcur = uic == ui;
                    full_y_pred(okcur) = mean(glm_y(okcur));
                    full_sdur_pred(okcur) = mean(search_dur_y(okcur));
                end
                
                nboot = 2000;
                var_explained_fun = @(x,y) 1 - (var(x - y) / var(x));
                
                % what fraction of variance in search duration ACROSS
                % CONDITIONS can be explained by subtracting the effect of
                % this variable?
                umean_sdur = nans(size(urows,1),1);
                umean_resid = nans(size(urows,1),1);
                
                boot_umean_sdur = nans(size(urows,1),nboot);
                boot_umean_resid = nans(size(urows,1),nboot);
                for ui = 1:size(urows,1)
                    okcur = uic == ui;

                    cur_sdur = search_dur_y(okcur);
                    cur_glmy = glm_y(okcur);

                    umean_sdur(ui) = mean(cur_sdur);
                    umean_resid(ui) = mean(cur_sdur - cur_glmy);
                    
                    allbootid = randsample(numel(cur_sdur),numel(cur_sdur)*nboot,true);
                    allbootid = reshape(allbootid,numel(cur_sdur),nboot);
                    for booti = 1:nboot
                        curbootid = allbootid(:,booti);
                        boot_umean_sdur(ui,booti) = mean(cur_sdur(curbootid));
                        boot_umean_resid(ui,booti) = mean(cur_sdur(curbootid) - cur_glmy(curbootid));
                    end
                end
                fract_var_across_conds_explained = 1 - (var(umean_resid) / var(umean_sdur));
                boot_fract_var_across_conds_explained = 1 - (var(boot_umean_resid,[],1) ./ var(boot_umean_sdur,[],1));
                bootse_fract_var_across_conds_explained = std(boot_fract_var_across_conds_explained);
                bootci95_fract_var_across_conds_explained = quantile(boot_fract_var_across_conds_explained,[.025 .975]);

                bootstor2{mi, pvi} = boot_fract_var_across_conds_explained;
                
                etextn('lt',curtext);
                
                curstats = struct();
                
                curstats.dur_s.mean = mean(glm_y);
                curstats.dur_s.se = sem(glm_y);
                curstats.dur_s.sd = std(glm_y);
                curstats.dur_s.med = median(glm_y);
                
                curstats.dur_sd_s.mean = std(glm_y);
                curboot = bootstrp(2000,@std,glm_y);
                curstats.dur_sd_s.se = std(curboot);
                
                curstats.beta_r.mean = stats.beta(2);
                curstats.beta_r.se = stats.se(2);
                curstats.beta_r.p = stats.p(2);

                curstats.beta_p.mean = stats.beta(3);
                curstats.beta_p.se = stats.se(3);
                curstats.beta_p.p = stats.p(3);
                
                curstats.var_exp_across.mean = fract_var_across_conds_explained;
                curstats.var_exp_across.se = bootse_fract_var_across_conds_explained;
                curstats.var_exp_across.ci95 = bootci95_fract_var_across_conds_explained;
                
                plotstats{mi,pvi} = curstats;

                xlim([0 1] + [-1 +1]*0.025);
                if pvi == round(mean(1:numel(plotvar)))
                    etextn('ct',sprintf('monk %s',plotmonk{mi}.name));
                end
                title(pv.name,'interpreter','none');
                if pvi == 1
                    xlabel('Environment reward probability');
                    ylabel('Duration (s)');
                    set(gca,'xtick',0:.25:1,'xticklabel',{'0%','','','','100%'});
                else
                    set(gca,'xtick',0:.25:1,'xticklabel',[]);
                end
                if pvi > 1
                    set(gca,'ycolor','none');
                end
            end
        end
        
        scale('eachax',h(:),'cm','hvscale',[1.1 1]);
        
        setlim(h(:),'ylim','tight');
        yl = ylim(h(1));
        yl(1) = 0;
        yl(2) = ceil(yl(2)*5)/5;
        setlim(h(:),'ylim',yl);
        
        %% Evaluate significance of bootstat comparisons
        for mi = 1:numel(plotmonk)
            tempVar = [];
            for pvi = 3:numel(plotvar)
                tempVar(pvi, 1:2) = (quantile(bootstor2{mi, 2} - bootstor2{mi, pvi}, [0.001 .999]));
            end
            test = sign(tempVar(:,1))==sign(tempVar(:,2));
            if all(test(3:5))
                disp(['For animal "' plotmonk{mi}.name '": bootstrap of the difference between all other factors and environment is significant for "Contribution to reward-related bias", indicated by the 99.9% CI of the difference excluding zero.']); 
            end
        end
        
        % Prep data for E
        betaRBootstor = {};
        betaPBootstor = {};
        for outerIter = 1:5
            for iter = 1:2000
                betaRBootstor{outerIter}(iter) = bootstor{outerIter}{iter}(2);
                betaPBootstor{outerIter}(iter) = bootstor{outerIter}{iter}(3);
            end
        end
        
        for mi = 1:numel(plotmonk)
            tempVar = [];
            for pvi = 3:numel(plotvar)
                tempVar(pvi, 1:2) = (quantile(betaRBootstor{mi, 2} - betaRBootstor{mi, pvi}, [0.001 .999]));
            end
            test = sign(tempVar(:,1))==sign(tempVar(:,2));
            if all(test(3:5))
                disp(['For animal "' plotmonk{mi}.name '": bootstrap of the difference between all other factors and environment is significant for "Reward effect on duration", indicated by the 99.9% CI of the difference excluding zero.']); 
            end
        end
        
        for mi = 1:numel(plotmonk)
            tempVar = [];
            for pvi = 3:numel(plotvar)
                tempVar(pvi, 1:2) = (quantile(betaPBootstor{mi, 2} - betaPBootstor{mi, pvi}, [0.001 .999]));
            end
            test = sign(tempVar(:,1))==sign(tempVar(:,2));
            if all(test(3:5))
                disp(['For animal "' plotmonk{mi}.name '": bootstrap of the difference between all other factors and environment is significant for "p(reward) effect on duration", indicated by the 99.9% CI of the difference excluding zero.']); 
            end
        end
        
        %% Plotting
        h = [];
        
        snames = {'var_exp_across','beta_r','beta_p'};
        sfullnames = {{'% reward-related search duration','variance explained'},{'Reward effect','on duration (s)'},{'p(reward) effect','on duration (s)'}};
        
        pvnames = {'env_gaze_dur_s','tar_gaze_dur_s','fix_gaze_dur_s','all_other_time_dur_s'};
        pvfullnames = {'Envir.','Target','Start pos.','Other'};
        
        for mi = 1:numel(plotmonk)
            for si = 1:numel(snames)
                sn = snames{si};
                
                h(mi,si)=nsubplot(4,numel(snames),4,si);
                if si == 1
                    ylabel(['monk ' plotmonk{mi}.name]);
                end
                if mi == 1
                    title(sfullnames{si},'interpreter','none');
                end
                x = [];
                xlab = {};
                for pvi = 1:numel(plotvar)
                    pv = plotvar{pvi};
                    if strcmp(pv.name,'search_dur_s')
                        continue;
                    end
                    
                    % convert fraction of variance explained to percent of
                    % variance explained
                    if ismember(sn,{'var_exp','var_exp_across'})
                        curmult = 100;
                    else
                        curmult = 1;
                    end
                    
                    if strcmp(pv.name,'env_gaze_dur_s')
                        curcolor = 'k';
                    else
                        curcolor = interpcolor('k','w',.5);
                    end
                    
                    cur = plotstats{mi,pvi}.(sn);
                    plot_errorbar(pvi,curmult.*cur.mean,curmult.*cur.se,{'bar',interpcolor(curcolor,'w',.8),'edgecolor',curcolor,'linewidth',2},{'color',curcolor,'linewidth',2});
                    
                    
                    curpvi = find(strcmp(pv.name,pvnames));
                    if numel(curpvi)==1
                        curpvname = pvfullnames{curpvi};
                    else
                        curpvname = pv.name;
                    end
                    
                    x(end+1) = pvi;
                    xlab{end+1} = curpvname;
                end
                setlim(gca,'ylim','tight',.1);
                if ismember(sn,{'var_exp','var_exp_across'})
                    ylim([0 100]);
                elseif ismember(sn,{'dur_s','dur_sd_s'})
                    setlim(gca,'ymin',0);
                end
                set(gca,'xtick',x,'xticklabel',xlab);
                xtickangle(gca,45);
            end
            
            cur_ids = find(ismember(snames,{'beta_r','beta_p'}));
            if ~isempty(cur_ids)
                setlim(h(mi,cur_ids),'ylim','tight',0,0,'square');
                yl = ylim(h(mi,cur_ids(1)));
                yl(2) = ceil(yl(2)*10)/10;
                yl(1) = -yl(2);
                set(h(mi,cur_ids),'ytick',[yl(1) 0 yl(2)],'yticklabel',[yl(1) 0 yl(2)]);
                setlim(h(mi,cur_ids),'ylim',yl + [-1 +1]*0.01);
            end
        end
    end
end

