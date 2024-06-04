function sx_sing = shot_stack(sufile,check,t0_shift,n_stack_theo,dS_theo,xseis_theo,sx_theo)

%%% S. Pasquet - V22.05.04
% Stack shot files
% sx_sing = shot_stack(sufile,check,t0_shift,n_stack_theo,dS_theo,xseis_theo,sx_theo)

wsl = ispc_wsl;

close all; drawnow;
stack_norm = 0;
j = 0;
start_i = 1;


if check == 1
    check_pos = 1;
    check_t0 = 0;
elseif check == 2
    check_pos = 0;
    check_t0 = 1;
elseif check == 3
    check_pos = 1;
    check_t0 = 1;
else
    check_pos = 0;
    check_t0 = 0;
end

% Get acquisition settings
[~,hdrs] = unix_cmd(['sugethw < ',sufile,' key=fldr,delrt,dt,ns,sx,gx,tracf output=geom'],wsl);
hdrs = str2num(hdrs);
fldrall = hdrs(:,1);
fldr_raw = unique(fldrall);

%%
ind_sx = NaN*ones(size(fldr_raw));
seismomat = cell(size(fldr_raw));
sx = ind_sx;
n_stack_raw = 0;
ax1 = []; ax2 = []; ax3 = []; fig1=[]; fig2=[];
flag = 0;

fprintf('\n   --------------------\n');
%%
for i = start_i:length(fldr_raw)

    % Read file
    [seismomat{i},xseis,~,tseis] = shotselec(sufile,fldr_raw(i),1,1,1,0);
    
    if exist('xseis_theo') && ~isempty(xseis_theo)
        xseis = xseis_theo;
    end
    
    tseis = tseis*1000;
    xseis = round(xseis*1000)/1000;
    dt = unique(round(diff(tseis)*1000)/1000);
    dx = unique(round(diff(xseis)*1000)/1000);
    
    ind_sx_std = (std(seismomat{i},[],2) == max(std(seismomat{i},[],2)));
    ind_sx(i) = find(ind_sx_std == max(ind_sx_std));
    
    std_ratio = std(seismomat{i},[],2)/max(std(seismomat{i},[],2));
    ind_std_ratio_lim = find(std_ratio > 0.9);
    
    % Get source position
    if exist('sx_theo') && ~isempty(sx_theo)
        sx(i) = sx_theo(i);
    else
        if length(ind_std_ratio_lim) <= 2
            sx(i) = round(xseis(ind_sx(i))*100)/100;
        else
            n_stack_raw = n_stack_raw + 1;
            if n_stack_raw == 1
                sx(i) = sx(find(~isnan(sx),1,'last')) + dS_theo;
            elseif n_stack_raw < n_stack_theo
                sx(i) = sx(find(~isnan(sx),1,'last'));
            else
                sx(i) = sx(find(~isnan(sx),1,'last'));
                n_stack_raw = 0;
            end
        end
    end
    sx_old(i) = sx(i);
    
    if isnan(sx(i))
        keyboard
    end

    file_dat = strcat(num2str(fldr_raw(i)),'.dat');
    file_su = strcat(num2str(fldr_raw(i)),'.su');
    
    seismomat_tmp = seismomat{i};
    seismomat_norm = bsxfun(@rdivide,seismomat{i},max(abs(seismomat{i}),[],2));
    trace_zero = seismomat_norm(abs(xseis-sx(i))<=1e-6,:)';
    
    if check > 0
        if ~isempty(trace_zero)
            [fig1,ax1] = plot_wiggle(1,-[0*trace_zero trace_zero 0*trace_zero],[sx(i)-dx sx(i) sx(i)+dx],tseis,1,1,99,16,'X (m)','Time (ms)',[sx(i)-5*dx/8 sx(i)+5*dx/8],[-50*dt +dt*250],[],[],[44 0 10 27],[]);
            set(1,'name',file_dat,'numbertitle','off');
        else
            if ishghandle(fig1)
                close(fig1);
            end
        end
        delete(ax2);
        [fig2,ax2] = plot_wiggle(2,-seismomat_norm',xseis,tseis,1,1,97,16,'X (m)','Time (ms)',[],[max([min(tseis),-25]) max(tseis)/3],[],[],[0 0 42 27],[]);
        hold on; plot(sx(i),0,'r.','markersize',15);
        set(2,'name',file_dat,'numbertitle','off');
        drawnow; hold off;
    end

    % Check if shot position is OK
    if check_pos == 1
        answer1 = inputdlg({'Shot position (m)'},'Cancel to discard shot',[1 35],{num2str(sx(i))});
        if ~isempty(answer1)
            if str2double(answer1(1)) ~= sx(i)
                n_stack_raw = n_stack_raw - 1;
                sx(i) = str2double(answer1(1)); % Shot position
                seismomat_tmp = seismomat{i};
                seismomat_norm = bsxfun(@rdivide,seismomat{i},max(abs(seismomat{i}),[],2));
                trace_zero = seismomat_norm(abs(xseis-sx(i))<=1e-6,:)';
                if ishghandle(fig1)
                    close(fig1);
                end
                if isempty(trace_zero)
                    trace_zero = seismomat_norm(abs(xseis-sx(i)) == min(abs(xseis-sx(i))),:)';
                    fprintf('\n   No trace at shot position - Display closest trace');
                    flag = 1;
                end
                [fig1,ax1] = plot_wiggle(1,-[0*trace_zero trace_zero 0*trace_zero],[sx(i)-dx sx(i) sx(i)+dx],tseis,1,1,99,16,'X (m)','Time (ms)',[sx(i)-5*dx/8 sx(i)+5*dx/8],[-50*dt +dt*250],[],[],[44 0 10 27],[]);
                set(1,'name',file_dat,'numbertitle','off');
                
                delete(ax2);
                [fig2,ax2] = plot_wiggle(2,-seismomat_norm',xseis,tseis,1,1,97,16,'X (m)','Time (ms)',[],[max([min(tseis),-25]) max(tseis)/3],[],[],[0 0 42 27],[]);
                hold on; plot(sx(i),0,'r.','markersize',15);
                set(2,'name',file_dat,'numbertitle','off');
                drawnow; hold off;
            end
            sx(i) = str2double(answer1(1)); % Shot position
            fprintf('\n   Keep shot %s at %1.2f m\n',file_dat,sx(i));
            sx_old(i) = sx(i);
        else
            sx_old(i) = sx(i);
            sx(i) = NaN;
            fprintf('\n   Discard shot %s at %1.2f m\n',file_dat,sx(i));
            delete(file_su);
        end
    else
        fprintf('\n   Keep shot %s at %1.2f m\n',file_dat,sx(i));
    end
    
    if i>start_i
        sx_prev = sx(1:i-1);
        sx_last = sx_prev(find(~isnan(sx_prev),1,'last'));
    else
        sx_prev = [];
        sx_last = [];
    end
    
    stack_cond = ((i>start_i && sx(i) ~= sx_last && ~isnan(sx(i)) && sx_old(i)~=sx_old(i-1)) || i == length(fldr_raw));
    
    if ~stack_cond
        if t0_shift == 1 && ~isnan(sx(i))
            if isempty(trace_zero)
                fprintf('\n   No trace at shot position - Keep original t0');
                flag = 1;
            end
            seismomat_shift{i} = zeros(size(seismomat{i}));
            if flag == 0
                t_zero = autopick_trace(trace_zero,tseis);
                if check_t0 == 1
                    fprintf('   Pick t0 for %s\n',file_dat);
                    
                    figure(fig1);
                    hold on; plot(sx(i),t_zero,'b+','markersize',8); drawnow;
                    a2 = magnify3(fig1,10,0.5,0.1); hold on; pause(1); % Magnifier
                    [~,tpick] = matpick_trace(sx(i),t_zero);
                    hold off;
                    
                    if isnan(tpick)
                        fprintf('   Discard shot %s at %1.2f m\n',file_dat,sx(i));
                        sx(i) = NaN;
                        delete(file_su);
                    elseif ~isempty(tpick)
                        t_zero = tpick;
                    end
                end
            else
                t_zero = 0;
                flag = 0;
            end
            
            ind_shift = round(t_zero/dt);
            
            if ind_shift > 0
                seismomat_ok = seismomat_tmp(:,ind_shift:end);
                ind_end = size(seismomat_ok,2);
                seismomat_shift{i}(:,1:ind_end) = seismomat_ok;
            else
                seismomat_ok = seismomat_tmp(:,1:end+ind_shift);
                ind_end = size(seismomat_ok,2);
                seismomat_shift{i}(:,-ind_shift+1:end) = seismomat_ok;
            end
            
        end

    else
        j = j + 1;
        filestack = strcat('shot',num2str(j,'%03d'),'_stack.su');
        fprintf('   Stack previous shots in %s\n',filestack);
        
        sx_sing(j) = sx_prev(find(~isnan(sx_prev),1,'last'));
        ind_seismomat = find(sx == sx_sing(j));
        n_stack(j) = length(ind_seismomat);
        seismomat_stack{j} = zeros(size(seismomat{i}));
        
        if i == length(fldr_raw)
            seismomat_shift{i} = zeros(size(seismomat{i}));
        end
                
        fldr_ok = [];
        for k = 1:n_stack(j)
            if t0_shift == 1
                seismomat_tmp2 = seismomat_shift{ind_seismomat(k)};
                seismomat_norm = bsxfun(@rdivide,seismomat_tmp2,max(abs(seismomat_tmp2),[],2));
            else
                seismomat_tmp2 = seismomat{ind_seismomat(k)};
                seismomat_norm = bsxfun(@rdivide,seismomat_tmp2,max(abs(seismomat_tmp2),[],2));
            end
            
            if stack_norm == 1
                seismomat_tmp2 = seismomat_norm;
            end
            
            seismomat_stack{j} = seismomat_stack{j} + seismomat_tmp2;
            seismomat_stack_norm = bsxfun(@rdivide,seismomat_stack{j},max(abs(seismomat_stack{j}),[],2));
            fldr_ok = [fldr_ok fldr_raw(ind_seismomat(k))];
            
            trace_zero_stack = seismomat_stack_norm(abs(xseis-sx_sing(j))<=1e-6,:)';
            if isempty(trace_zero_stack)
               t_zero_final(j) = 0;
            else
                t_zero_final(j) = autopick_trace(trace_zero_stack,tseis);
            end
        end
        
        if check > 0
            delete(ax3);
            [fig3, ax3] = plot_wiggle(3,-seismomat_stack_norm',xseis,tseis,1,1,99,16,'X (m)','Time (ms)',[],[max([min(tseis),-25]) max(tseis)/3],[],[],[0 0 42 27],[]);
            hold on; plot(sx_sing(j),t_zero_final(j),'r.','markersize',20);
            set(3,'name',sprintf('Stacked shot at %1.2f m',sx_sing(j)),'numbertitle','off');
            hold off; drawnow;
        end

        file_raw = strcat(num2str(fldr_ok(1)),'.su');
        filestack_tmp = strcat('shot',num2str(j,'%03d'),'_stack_tmp.su');
        file1_tmp_txt = strcat('shot',num2str(j,'%03d'),'_stack_tmp.txt');
        file1_tmp_bin = strcat('shot',num2str(j,'%03d'),'_stack_tmp.bin');
        file_head = strcat(num2str(ind_seismomat(1)),'_head.su');
        dlmwrite(file1_tmp_txt,seismomat_stack{j},'delimiter',' ');
        [~,~] = unix_cmd(sprintf('a2b < %s > %s n1=%d',file1_tmp_txt,file1_tmp_bin,length(xseis)),wsl);
        [~,~] = unix_cmd(sprintf('sustrip < %s head=%s > strip.su',file_raw,file_head),wsl);
        [~,~] = unix_cmd(sprintf('supaste < %s ns=%i head=%s > %s',file1_tmp_bin,length(tseis),file_head,filestack),wsl);
        delete('strip.su',file1_tmp_txt,file1_tmp_bin,file_head);

        for k = 1:length(fldr_ok)
            file_raw = strcat(num2str(fldr_ok(k)),'.su');
            delete(file_raw);
        end
        
        if stack_norm == 1
            unix_cmd(sprintf('sushw < %s key=sx,fldr a=%i,%i | suop op=norm > %s',filestack,sx_sing(j),j,filestack_tmp),wsl);
        else
            unix_cmd(sprintf('sushw < %s key=sx,fldr a=%i,%i > %s',filestack,sx_sing(j),j,filestack_tmp),wsl);
        end
        movefile(filestack_tmp,filestack);

        if t0_shift == 1 && ~isnan(sx(i)) && i ~= length(fldr_raw)
            if isempty(trace_zero)
                fprintf('\n   No trace at shot position - Keep original t0');
                flag = 1;
            end
            seismomat_shift{i} = zeros(size(seismomat{i}));
            if flag == 0
                t_zero = autopick_trace(trace_zero,tseis);
                
                if check_t0 == 1
                    fprintf('   Pick t0 for %s\n',file_dat);
                    
                    figure(fig1)
                    hold on; plot(sx(i),t_zero,'b+','markersize',8); drawnow;
                    a2 = magnify3(fig1,10,0.5,0.1); hold on; pause(1); % Magnifier
                    [~,tpick] = matpick_trace(sx(i),t_zero);
                    hold off;
                    
                    if isnan(tpick)
                        fprintf('   Discard shot %s at %1.2f m\n',file_dat,sx(i));
                        sx(i) = NaN;
                        delete(file_su);
                    elseif ~isempty(tpick)
                        t_zero = tpick;
                    end
                end
            else
                t_zero = 0;
                flag = 0;
            end
 
            ind_shift = round(t_zero/dt);
            if ind_shift > 0
                seismomat_ok = seismomat_tmp(:,ind_shift:end);
                ind_end = size(seismomat_ok,2);
                seismomat_shift{i}(:,1:ind_end) = seismomat_ok;
            else
                seismomat_ok = seismomat_tmp(:,1:end+ind_shift);
                ind_end = size(seismomat_ok,2);
                seismomat_shift{i}(:,-ind_shift+1:end) = seismomat_ok;
            end
        end
        
        %%
        if i==length(fldr_raw) && length(fldr_raw)-sum(n_stack)==1

            sx_prev = sx(1:i);
            sx_last = sx_prev(find(~isnan(sx_prev),1,'last'));
            
            j = j + 1;
            filestack = strcat('shot',num2str(j,'%03d'),'_stack.su');
            fprintf('   Stack previous shots in %s\n',filestack);
            
            sx_sing(j) = sx_prev(find(~isnan(sx_prev),1,'last'));
            ind_seismomat = find(sx == sx_sing(j));
            n_stack(j) = length(ind_seismomat);
            seismomat_stack{j} = zeros(size(seismomat{i}));
            
            if i == length(fldr_raw)
                seismomat_shift{i} = zeros(size(seismomat{i}));
            end
            
            fldr_ok = [];
            for k = 1:n_stack(j)
                if t0_shift == 1
                    seismomat_tmp2 = seismomat_shift{ind_seismomat(k)};
                    seismomat_norm = bsxfun(@rdivide,seismomat_tmp2,max(abs(seismomat_tmp2),[],2));
                else
                    seismomat_tmp2 = seismomat{ind_seismomat(k)};
                    seismomat_norm = bsxfun(@rdivide,seismomat_tmp2,max(abs(seismomat_tmp2),[],2));
                end
                
                if stack_norm == 1
                    seismomat_tmp2 = seismomat_norm;
                end
                
                seismomat_stack{j} = seismomat_stack{j} + seismomat_tmp2;
                seismomat_stack_norm = bsxfun(@rdivide,seismomat_stack{j},max(abs(seismomat_stack{j}),[],2));
                fldr_ok = [fldr_ok fldr_raw(ind_seismomat(k))];
                
                trace_zero_stack = seismomat_stack_norm(abs(xseis-sx_sing(j))<=1e-6,:)';
                if isempty(trace_zero_stack)
                    t_zero_final(j) = 0;
                else
                    t_zero_final(j) = autopick_trace(trace_zero_stack,tseis);
                end
            end
            
            if check > 0
                delete(ax3);
                [fig3, ax3] = plot_wiggle(3,-seismomat_stack_norm',xseis,tseis,1,1,99,16,'X (m)','Time (ms)',[],[max([min(tseis),-25]) max(tseis)/3],[],[],[0 0 42 27],[]);
                hold on; plot(sx_sing(j),t_zero_final(j),'r.','markersize',20);
                set(3,'name',sprintf('Stacked shot at %1.2f m',sx_sing(j)),'numbertitle','off');
                hold off; drawnow;
            end
            
            file_raw = strcat(num2str(fldr_ok(1)),'.su');
            filestack_tmp = strcat('shot',num2str(j,'%03d'),'_stack_tmp.su');
            file1_tmp_txt = strcat('shot',num2str(j,'%03d'),'_stack_tmp.txt');
            file1_tmp_bin = strcat('shot',num2str(j,'%03d'),'_stack_tmp.bin');
            file_head = strcat(num2str(ind_seismomat(1)),'_head.su');
            dlmwrite(file1_tmp_txt,seismomat_stack{j},'delimiter',' ');
            [~,~] = unix_cmd(sprintf('a2b < %s > %s n1=%d',file1_tmp_txt,file1_tmp_bin,length(xseis)),wsl);
            [~,~] = unix_cmd(sprintf('sustrip < %s head=%s > strip.su',file_raw,file_head),wsl);
            [~,~] = unix_cmd(sprintf('supaste < %s ns=%i head=%s > %s',file1_tmp_bin,length(tseis),file_head,filestack),wsl);
            delete('strip.su',file1_tmp_txt,file1_tmp_bin,file_head);
            
            for k = 1:length(fldr_ok)
                file_raw = strcat(num2str(fldr_ok(k)),'.su');
                delete(file_raw);
            end
            
            if stack_norm == 1
                unix_cmd(sprintf('sushw < %s key=sx,fldr a=%i,%i | suop op=norm > %s',filestack,sx_sing(j),j,filestack_tmp),wsl);
            else
                unix_cmd(sprintf('sushw < %s key=sx,fldr a=%i,%i > %s',filestack,sx_sing(j),j,filestack_tmp),wsl);
            end
            movefile(filestack_tmp,filestack);
            
            if t0_shift == 1 && ~isnan(sx(i)) && i ~= length(fldr_raw)
                if isempty(trace_zero)
                    fprintf('\n   No trace at shot position - Keep original t0');
                    flag = 1;
                end
                seismomat_shift{i} = zeros(size(seismomat{i}));
                if flag == 0
                    t_zero = autopick_trace(trace_zero,tseis);
                    
                    if check_t0 == 1
                        fprintf('   Pick t0 for %s\n',file_dat);
                        
                        figure(fig1)
                        hold on; plot(sx(i),t_zero,'b+','markersize',8); drawnow;
                        a2 = magnify3(fig1,10,0.5,0.1); hold on; pause(1); % Magnifier
                        [~,tpick] = matpick_trace(sx(i),t_zero);
                        hold off;
                        
                        if isnan(tpick)
                            fprintf('   Discard shot %s at %1.2f m\n',file_dat,sx(i));
                            sx(i) = NaN;
                            delete(file_su);
                        elseif ~isempty(tpick)
                            t_zero = tpick;
                        end
                    end
                else
                    t_zero = 0;
                    flag = 0;
                end
                
                ind_shift = round(t_zero/dt);
                if ind_shift > 0
                    seismomat_ok = seismomat_tmp(:,ind_shift:end);
                    ind_end = size(seismomat_ok,2);
                    seismomat_shift{i}(:,1:ind_end) = seismomat_ok;
                else
                    seismomat_ok = seismomat_tmp(:,1:end+ind_shift);
                    ind_end = size(seismomat_ok,2);
                    seismomat_shift{i}(:,-ind_shift+1:end) = seismomat_ok;
                end
            end
        end
        
    end
    if check > 0 && ~isempty(fig1) && ishghandle(fig1)
        close(fig1);
    end
    fprintf('\n   --------------------\n');
end
unix_cmd(sprintf('cat shot*_stack.su > %s',sufile),wsl);
delete('shot*_stack.su');
close all;