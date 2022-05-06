function plot_targ(targetfile,pickcol1,fMIN,fMAX,fticks,freqtitle_long,axetop,Flogscale,VphMIN,VphMAX,fs,imgform,imgres)

%%% S. Pasquet - V16.11.30
% Quick plot of dispersion curves
% plot_targ(targetfile,pickcol1,fMIN,fMAX,fticks,freqtitle_long,axetop,Flogscale,VphMIN,VphMAX,fs,imgform,imgres)

if exist('targetfile','var')==0 || isempty(targetfile)==1
    [targetfile,targetpath]=uigetfile('*.target','Select dispersion curve or cancel');
    if targetfile==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select a target file');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
else
    targetpath=[];
end

if exist('pickcol1','var')==0 || isempty(pickcol1)==1
   pickcol1='k'; 
end

run('SWIP_defaultsettings')


% Read target file to get picked dispersion curves
[freqresamp,vresamp,deltaresamp,modes]=targ2pvc(fullfile(targetpath,targetfile));
npvc=length(modes);
for ip=1:npvc
    if ip==1
        if eb==1
            [fig,han]=plot_curv([],freqresamp{modes(ip)+1},vresamp{modes(ip)+1},...
                deltaresamp{modes(ip)+1},'.-',pickcol1,[],axetop,axerev,...
                0,fs,freqtitle_long,'Phase velocity (m/s)',[],[fMIN fMAX],[VphMIN VphMAX],[],...
                fticks,Vphticks,[],[],[],[1 1 24 18],[],[],Flogscale);
            xlimits=xlim;
            tick_length=diff(xlimits)/100;
            errorbar_tick(han,tick_length,'units');
        else
            fig=plot_curv([],freqresamp{modes(ip)+1},vresamp{modes(ip)+1},[],'.-',pickcol1,[],axetop,axerev,...
                0,fs,freqtitle_long,'Phase velocity (m/s)',[],[fMIN fMAX],[VphMIN VphMAX],[],...
                fticks,Vphticks,[],[],[],[1 1 24 18],[],[],Flogscale);
        end
        hold on
    else
        figure(fig);
        plot(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},'.-','Color',pickcol1,...
            'linewidth',1.5,'markersize',10);
        if eb==1
            if str2double(matrelease(1:4))>2014
                han=terrorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},1,'units');
                set(han,'LineWidth',1.5,'Color',pickcol1)
            else
                han=errorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},...
                    '.-','Color',pickcol1,'linewidth',1.5,'markersize',10);
                xlimits=xlim;
                tick_length=diff(xlimits)/100;
                errorbar_tick(han,tick_length,'units');
            end
        end
    end
end
hold off

file1=[targetfile(1:end-6),'curv.',imgform];
save_fig(fig,file1,imgform,imgres,1);

end