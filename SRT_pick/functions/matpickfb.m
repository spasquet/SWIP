function [xpick,tpick,closefig,shotprev,percnext,tshift,clean_air,scalnext]=matpickfb(xseis,tseis,seismomat,xpick,tpick,tshift,sx)
%%% S. Pasquet - V17.10.24

assist = 0;
assist_win = 10;

hold on;
button='';
i=0; closefig=0; shotprev=0; percnext=[]; clean_air=[]; scalnext=[];
ha=findobj(gcf,'type','axes');
% h=get(ha(1),'children');
% dt = mean(diff(tseis));

if nargin<3
    xpick=[]; tpick=[];
else
    hold on;
    h0=line(xpick,tpick,'color','r','linestyle','none','marker','+','markersize',12,'linewidth',2,'parent',ha(1));
    hold on;
    h1=line(xpick,tpick,'color','r','linestyle','none','marker','+','markersize',8,'linewidth',2,'parent',ha(2));
    i=length(xpick);
    h2=line([min(xseis) max(xseis)],[tshift tshift],'color','b','linewidth',1.5,'parent',ha(1));
    h3=line([min(xseis) max(xseis)],[tshift tshift],'color','b','linewidth',1,'parent',ha(2));
    
end
drawnow

fprintf('\n  ENTER or close window : Save current picks and go to next Xmid');
fprintf('\n  BACKSPACE : Save current picks and go to previous Xmid');
fprintf('\n  W : Save current picks and stop script');
fprintf('\n  A : Discard current picks and go to next Xmid (keep previous picks)');
fprintf('\n  Z : Discard current picks and go to previous Xmid (keep previous picks)');
fprintf('\n  X : Discard current picks and stop script (keep previous picks)');
fprintf('\n  D : Delete one or several points');
fprintf('\n  R : Reset all picks');
fprintf('\n  G : Change gain value');
fprintf('\n  T : Set t0 for delayed triggering');
fprintf('\n  H : Help to determine t0 with airwave');
fprintf('\n  C : Airwave cleaning\n');


fprintf('\n  Picking mode\n');

while strcmp (button , '') == 1
    i=i+1;
    try
        [Xpick, Tpick, button] = pick_fig;
        
        if isempty(Xpick)==0
            test=abs(xseis-Xpick);
            Xpick=xseis(test==min(test));
            
            if assist == 1
                ind_x = xseis==Xpick;
                trace_amp = seismomat(ind_x,:);
                trace_grad = gradient(trace_amp);
                max_grad = max(abs(trace_grad(tseis>0)));
                trace_grad = trace_grad/max_grad;
                ind_amp = (trace_amp>0 & tseis >= Tpick-assist_win & tseis <= Tpick+assist_win);
                ind_grad = (abs(trace_grad) > 0.02 & tseis >= Tpick-assist_win & tseis <= Tpick+assist_win);
                ind_ok = find(ind_amp & ind_grad);
                keyboard
                if ~isempty(tseis(ind_ok(1)))
                    Tpick = tseis(ind_ok(1));
                end
                
                %                 figure
                %                 plot(tseis,trace_amp,'rx');
                %                 hold on;
                %                 plot(tseis(ind_ok),trace_grad(ind_ok),'bx');
                %                 plot(Tpick,0,'rx')
                %
                
            end
            %             if Tpick - tshift < 0
            %                 fprintf('\n  !!! Cannot pick negative time - First arrival set to %1.0f ms !!!\n',0);
            %                 Tpick = tshift;
            %             end
        end
    catch
        closefig=1;
        fprintf('\n  Saving picked first arrivals - Go to next shot\n');
        return
    end
    
    if strcmp (button , 'return') == 1 % Press ENTER to end picking, save pick file and go to next shot
        if isempty(tpick)==1
            fprintf('\n  No first arrival picked - Go to next shot\n');
        else
            fprintf('\n  Saving picked first arrivals - Go to next shot\n');
        end
        break
    elseif strcmp (button , 'backspace') == 1 % Press BACKSPACE to end picking, save pick file and go to previous shot
        if isempty(tpick)==1
            fprintf('\n  No first arrival picked - Go to previous shot\n');
        else
            fprintf('\n  Saving picked first arrivals - Go to previous shot\n');
        end
        shotprev=1;
        break
    elseif strcmp (button , 'w') == 1 % Press W to end picking, save pick file and stop script
        if isempty(tpick)==1
            fprintf('\n  No first arrival picked - Stopping script\n');
        else
            fprintf('\n  Saving picked first arrivals - Stopping script\n');
        end
        shotprev=-1;
        break
    elseif strcmp (button , 'a') == 1 % Press A to abort, keep existing pick file and go to next shot
        fprintf('\n  Abort and keep existing picks - Go to next shot\n');
        tpick=0;
        break
    elseif strcmp (button , 'x') == 1 % Press X to abort, keep existing pick file and stop script
        fprintf('\n  Abort and keep existing picks - Stopping script\n');
        tpick=0; shotprev=-1;
        break
    elseif strcmp (button , 'z') == 1 % Press Z to abort, keep existing pick file and go to previous shot
        fprintf('\n  Abort and keep existing picks - Go to previous shot\n');
        tpick=0; shotprev=1;
        break
    elseif strcmp (button , 'p') == 1 % Press P to change percentile
        if isempty(tpick)==1
            fprintf('\n  No first arrival picked - Changing clip percentile\n');
        else
            fprintf('\n  Saving picked first arrivals - Changing clip percentile\n');
        end
        answer1=inputdlg({'New clip percentile'},'',1);
        if isempty(answer1)==1
            percnext=[];
        else
            percnext=str2double(answer1(1));
        end
        break
    elseif strcmp (button , 'g') == 1 % Press G to change gain
        if isempty(tpick)==1
            fprintf('\n  No first arrival picked - Changing gain\n');
        else
            fprintf('\n  Saving picked first arrivals - Changing gain\n');
        end
        answer1=inputdlg({'New gain'},'',1);
        if isempty(answer1)==1
            scalnext=[];
        else
            scalnext=str2double(answer1(1));
        end
        break
    elseif strcmp (button , 'd') == 1 % Press D to delete points
        if i==1
            i=i-1; button='';
            continue
        end
        indselec=[]; button='';
        fprintf('\n  Deletion mode');
        fprintf('\n  Click on the points you want to delete, press d or ENTER to go back to picking mode\n');
        while strcmp (button , 'return') ~= 1 && strcmp (button , 'd') ~= 1
            try
                [Xdel,~,button]=pick_fig;
                if isempty(Xdel)==0
                    test=abs(xseis-Xdel);
                    Xdel=xseis(test==min(test));
                    indselec=find(xpick==Xdel);
                end
                %                 indselec=selectdata('action','delete',...
                %                     'selectionmode','closest','ignore',h);
            catch
                closefig=1;
                return
            end
            if isempty(indselec)==0
                tpick(xpick==xpick(indselec))=[];
                xpick(xpick==xpick(indselec))=[];
                set(h0,'xdata',xpick,'ydata',tpick);
                set(h1,'xdata',xpick,'ydata',tpick);
                indselec=[];
            end
        end
        button=''; fprintf('\n  Exiting deletion mode - Back to picking mode\n');
        continue
    elseif strcmp (button , 'r') == 1; % Press R to reset picks
        fprintf('\n  Reset all picks\n');
        button=''; i=1;
        xpick=[]; tpick=[]; tshift=0;
        hold on
        if exist('h0','var')==1
            delete(h0,h1); clear h0;
        end
        if exist('h2','var')==1
            delete(h2,h3); clear h2 h3;
        end
        hold on;
        continue
    elseif strcmp (button , 't') == 1; % Press T to set t0
        fprintf('\n  Set t0 for delayed triggering\n');
        [~,tshift] = pick_fig;
        if isempty(tshift)
            tshift = 0;
        end
        %         if ismember(sx,xpick)
        %             tshift = tpick(xpick==sx)-dt;
        %             fprintf('\n  tshift = %1.3f ms\n',tshift);
        %         end
        if exist('h2','var')==1
            delete(h2,h3);
        end
        button=''; i=i-1;
        continue
    elseif strcmp (button , 'h') == 1; % Press H to pick airwave
        fprintf('\n  Pick 4 airwave arrivals\n');
        
        fprintf('\n  Pick first airwave arrival');
        [xair1,tair1]=pick_fig;
        test=abs(xseis-xair1);
        xair1=xseis(test==min(test));
        
        fprintf('\n  Pick second airwave arrival');
        [xair2,tair2]=pick_fig;
        test=abs(xseis-xair2);
        xair2=xseis(test==min(test));
        
        fprintf('\n  Pick third airwave arrival');
        [xair3,tair3]=pick_fig;
        test=abs(xseis-xair3);
        xair3=xseis(test==min(test));
        
        fprintf('\n  Pick fourth airwave arrival');
        [xair4,tair4]=pick_fig;
        test=abs(xseis-xair4);
        xair4=xseis(test==min(test));
        
        [xair_all,I] = unique([xair1 xair2 xair3 xair4]);
        tair_all = [tair1 tair2 tair3 tair4];
        tair_all = tair_all(I);
        
        %         p = polyfit(xair_all,tair_all,1);
        p = mmpolyfit(xair_all,tair_all,1,'equal',[1 0 1000/343]);
        tair_interp = polyval(p,xseis);
        hold on;
        plot(xseis,tair_interp,'c-','linewidth',2,'parent',ha(1));
        ha(2); hold on;
        plot(xseis,tair_interp,'c-','linewidth',1.5,'parent',ha(2));
        
        %         if ismember(sx,xseis)
        tshift = polyval(p,sx);
        fprintf('\n  tshift = %1.3f ms\n',tshift);
        if exist('h2','var')==1
            delete(h2,h3);
        end
        %         end
        
        button=''; i=i-1;
        continue
    elseif strcmp (button , 'c') == 1; % Press C to clean airwave
        fprintf('\n  Selected airwave start\n');
        [xair,tair1]=pick_fig;
        fprintf('\n  Selected airwave end\n');
        [~,tair2]=pick_fig;
        test=abs(xseis-xair);
        xair=xseis(test==min(test));
        clean_air = [xair,tair1,tair2];
        percnext = 0;
        scalnext = 0;
        break
    elseif strcmp (button , '') ~= 1;
        button=''; i=i-1;
        continue
    end
    if ismember(Xpick,xpick)==0
        try
            xpick=[xpick,Xpick];
            tpick=[tpick,Tpick];
        catch
            xpick=[xpick;Xpick];
            tpick=[tpick;Tpick];
        end
    else
        try
            tpick(xpick==Xpick)=Tpick;
            xpick(xpick==Xpick)=Xpick;
        catch
            keyboard
        end
    end
    hold on
    if exist('h0','var')==1
        delete(h0,h1);
    end
    hold on;
    h0=line(xpick,tpick,'color','r','linestyle','none','marker','+','markersize',12,'linewidth',2,'parent',ha(1));
    hold on;
    h1=line(xpick,tpick,'color','r','linestyle','none','marker','+','markersize',8,'linewidth',2,'parent',ha(2));
    h2=line([min(xseis) max(xseis)],[tshift tshift],'color','b','linewidth',1.5,'parent',ha(1));
    h3=line([min(xseis) max(xseis)],[tshift tshift],'color','b','linewidth',1,'parent',ha(2));
end
if length(tpick)>1
    tpick = tpick-tshift;
end
if length(tpick)==1 && tpick~=0
    tpick = tpick-tshift;
end
end
