function dat2dsp(dspmat,f,v,flip,dspfile_new,dspfile_old)

%%% S. Pasquet - V22.05.04
% Convert matrix in .dsp SU file
% dat2dsp(dspmat,f,v,flip,dspfile_new,dspfile_old)

wsl = ispc_wsl;

dspfile_mat=[dspfile_new,'.dat'];
dspfile_tmp=[dspfile_new,'.tmp'];
dspfile_tmp2=[dspfile_new,'.tmp2'];
dspfile_head=[dspfile_new,'.head'];

nf = length(f);
nv = length(v);
if flip==0
    n1=nf;
    dlmwrite(dspfile_mat,dspmat','delimiter','\t','precision','%.6f');
else
    n1=nv;
    dlmwrite(dspfile_mat,dspmat,'delimiter','\t','precision','%.6f');
end
com1=sprintf('a2b < %s > %s n1=%d',dspfile_mat,dspfile_tmp,n1);
[~,~]=unix_cmd(com1,wsl);
delete(dspfile_mat);

com1=sprintf('sustrip < %s > %s head=%s',dspfile_old,dspfile_tmp2,dspfile_head);
[~,~]=unix_cmd(com1,wsl);
delete(dspfile_tmp2);

com1=sprintf('supaste < %s > %s ns=%d head=%s',dspfile_tmp,dspfile_new,n1,dspfile_head);
[~,~]=unix_cmd(com1,wsl);
delete(dspfile_tmp,dspfile_head);
