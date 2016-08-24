function [dir_all,dir_inv_img]=dir_create(flag,nWmin,nWmax,dW,dSmin,dSmax,side)

%%% S. Pasquet - V16.6.28
% Create directories for surface-wave analysis

dir_inv_img=[];
dir_all.dir_start=pwd; % Current working directory

if flag~=1
    % Select main folder (Wmin_max.dWb.nSmin_max.side)
    dir_all.dir_main=uigetdir('./','Select main folder');
    if dir_all.dir_main==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select a SWIP folder');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    [~,dir_all.dir_main,ext]=fileparts(dir_all.dir_main);
    dir_all.dir_main=[dir_all.dir_main,ext];
else
    dir_all.dir_main=['W',num2str(nWmin),'_',num2str(nWmax),'.dW'...
        ,num2str(dW),'.dS',num2str(dSmin),'_',num2str(dSmax),'.',side];
    if exist(dir_all.dir_main,'dir')~=7
        mkdir(dir_all.dir_main);
    end
end

% Path of file.dat folder (contains processing files)
dir_all.dir_dat=fullfile(dir_all.dir_main,'/file.dat/');
% Path of file.pick folder (contains picked dispersion curves ASCII files)
dir_all.dir_pick=fullfile(dir_all.dir_main,'/file.pick/');
% Path of file.img folder to store .img files
dir_all.dir_img=fullfile(dir_all.dir_main,'/file.img/');
% Path of file.xmid folder to store .img files
dir_all.dir_img_xmid=fullfile(dir_all.dir_img,'/1D_data/');
% Path of file.targ folder to store .targ files
dir_all.dir_targ=fullfile(dir_all.dir_main,'/file.targ/');
% Path of file.param folder to store global .param files
dir_all.dir_param='file.param';
% Path of file.inv folder to store inversion results
dir_all.dir_inv=fullfile(dir_all.dir_main,'/file.inv/');
% Path of file.xzv folder to store final 2D models
dir_all.dir_xzv=fullfile(dir_all.dir_main,'/file.xzv/');

if exist(dir_all.dir_dat,'dir')~=7
    if flag==1
        mkdir(dir_all.dir_dat);
    else
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Missing file.dat folder');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        dir_all.dir_main=0;
        return
    end
end
if exist(dir_all.dir_pick,'dir')~=7
    mkdir(dir_all.dir_pick);
end
if exist(dir_all.dir_img,'dir')~=7
    mkdir(dir_all.dir_img);
end
if exist(dir_all.dir_img_xmid,'dir')~=7
    mkdir(dir_all.dir_img_xmid);
end
if exist(dir_all.dir_targ,'dir')~=7
    mkdir(dir_all.dir_targ);
end
if exist(dir_all.dir_param,'dir')~=7
    mkdir(dir_all.dir_param);
end
if exist(dir_all.dir_inv,'dir')~=7
    mkdir(dir_all.dir_inv);
end
if exist(dir_all.dir_xzv,'dir')~=7
    mkdir(dir_all.dir_xzv);
end

if flag==2
    % Select inversion folder containing report folders
    dir_inv_img.dir_rep_inv=uigetdir([dir_all.dir_main,'/file.inv/'],...
        'Select inversion folder');
    if dir_inv_img.dir_rep_inv==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select an inversion folder');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    % Get name of inversion folder
    if isunix==1
        dir_inv_img.dir_inv_name=dir_inv_img.dir_rep_inv...
            (strfind(dir_inv_img.dir_rep_inv,'/file.inv/')+10:end);
    else
        dir_inv_img.dir_inv_name=dir_inv_img.dir_rep_inv...
        (strfind(dir_inv_img.dir_rep_inv,'\file.inv\')+10:end);
    end
    % Folder to store inversion .img files
    dir_inv_img.dir_img_inv=[dir_all.dir_img,dir_inv_img.dir_inv_name];
    if exist(dir_inv_img.dir_img_inv,'dir')~=7
        mkdir(dir_inv_img.dir_img_inv);
    end
    % Folder to store inversion .xzv files
    dir_inv_img.dir_xzv_inv=[dir_all.dir_xzv,dir_inv_img.dir_inv_name];
    if exist(dir_inv_img.dir_img_inv,'dir')~=7
        mkdir(dir_inv_img.dir_img_inv);
    end
end

end