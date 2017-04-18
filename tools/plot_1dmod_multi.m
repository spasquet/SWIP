function plot_1dmod_multi(dir_mod,mappick,doi,vsMIN,vsMAX,fs,imgform,imgres)

%%% S. Pasquet - V17.03.29
% Quick plot of multiple 1D velocity models
% plot_1dmod_multi(dir_mod,mappick,doi,vsMIN,vsMAX,fs,imgform,imgres)

if exist('dir_mod','var')==0 || isempty(dir_mod)==1
    dir_mod=uigetdir('./','Select folder containing velocity files');
    if dir_mod==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select a velocity folder');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    modstruct=[dir(fullfile(dir_mod,'*.smooth')); dir(fullfile(dir_mod,'*.layered'));...
        dir(fullfile(dir_mod,'*.ridge')); dir(fullfile(dir_mod,'*.best')) ; dir(fullfile(dir_mod,'*vptomo')); dir(fullfile(dir_mod,'*.txt')); dir(fullfile(dir_mod,'*.dat'))];
end

if exist('doi','var')==0 || isempty(doi)==1
    doi=[];
end

if exist('mappick','var')==0 || isempty(mappick)==1
    mappick=jet(length(modstruct));
end

run('SWIP_defaultsettings')

for i=1:length(modstruct)
    modfile=modstruct(i).name;
    modvel=dlmread(fullfile(dir_mod,modfile),'',1,0);
    moddepth=[0;cumsum(modvel(:,1))];
    vssw=modvel(:,3);
    [~,str_tmp] = fileparts(modfile);
    if i==1
        str = {str_tmp};
        [VSplot,Zplot]=stair2plot(vssw,moddepth);
        [fig,han(i)]=plot_curv([],VSplot,Zplot,[],'-',mappick(i,:),[],1,1,...
            0,fs,'Vs (m/s)',depthtitle,[],[vsMIN vsMAX],[dpMIN dpMAX],[],...
            vsticks,dticks,[],doi,[],[0 0 24 18],[],0);
        hold on
    else
        str = [str ; {str_tmp}];
        [VSplot,Zplot]=stair2plot(vssw,moddepth);
        figure(fig);
        han(i)=plot(VSplot,Zplot,'-','Color',mappick(i,:),...
            'linewidth',2,'markersize',10);
        hold on
    end
end
hold off
legend(han,str);

file1=['mod1d_all.',imgform];
save_fig(fig,file1,imgform,imgres,1);

end
