function plot_targ_multi(dir_targ,mappick,nmodemax,fMIN,fMAX,fticks,freqtitle_long,axetop,Flogscale,VphMIN,VphMAX,fs,imgform,imgres)

%%% S. Pasquet - V16.11.30
% Quick plot of multiple dispersion curves
% plot_targ_multi(dir_targ,mappick,nmodemax,fMIN,fMAX,fticks,freqtitle_long,axetop,Flogscale,VphMIN,VphMAX,fs,imgform,imgres)

if exist('dir_targ','var')==0 || isempty(dir_targ)==1
    dir_targ=uigetdir('./','Select folder containing target files');
    if dir_targ==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select a target folder');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    targstruct=dir(fullfile(dir_targ,'*.target'));
end

if exist('mappick','var')==0 || isempty(mappick)==1
    mappick=jet(length(targstruct));
end

if exist('nmodemax','var')==0 || isempty(nmodemax)==1
    nmodemax=[];
end

run('SWIP_defaultsettings')

for i=1:length(targstruct)
    targetfile=targstruct(i).name;
    % Read target file to get picked dispersion curves
    [freqresamp,vresamp,deltaresamp,modes]=targ2pvc(fullfile(dir_targ,targetfile));
    npvc=length(modes);
    if npvc>nmodemax
        npvc=nmodemax;
    end
    for ip=1:npvc
        if mod(ip,2)==1 % non even
            symb='+-';
        else
            symb='.-';
        end
        if i==1 && ip==1
            if eb==1
                [fig,han]=plot_curv([],freqresamp{modes(ip)+1},vresamp{modes(ip)+1},...
                    deltaresamp{modes(ip)+1},symb,mappick(i,:),[],axetop,axerev,...
                    0,fs,freqtitle_long,'Phase velocity (m/s)',[],[fMIN fMAX],[VphMIN VphMAX],[],...
                    fticks,Vphticks,[],[],[],[1 1 24 18],[],[],Flogscale);
                xlimits=xlim;
                tick_length=diff(xlimits)/100;
                errorbar_tick(han,tick_length,'units');
            else
                fig=plot_curv([],freqresamp{modes(ip)+1},vresamp{modes(ip)+1},[],symb,mappick(i,:),[],axetop,axerev,...
                    0,fs,freqtitle_long,'Phase velocity (m/s)',[],[fMIN fMAX],[VphMIN VphMAX],[],...
                    fticks,Vphticks,[],[],[],[1 1 24 18],[],[],Flogscale);
            end
            hold on
        else
            figure(fig);
            plot(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},symb,'Color',mappick(i,:),...
                'linewidth',1.5,'markersize',10);
            if eb==1
                if str2double(matrelease(1:4))>2014
                    han=terrorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},1,'units');
                    set(han,'LineWidth',1.5,'Color',mappick(i,:))
                else
                    han=errorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},...
                        symb,'Color',mappick(i,:),'linewidth',1.5,'markersize',10);
                    xlimits=xlim;
                    tick_length=diff(xlimits)/100;
                    errorbar_tick(han,tick_length,'units');
                end
            end
        end
    end
end
hold off

file1=['target.curv.',imgform];
save_fig(fig,file1,imgform,imgres,1);

end