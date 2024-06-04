function matop(infile,op,flip,outfile)

%%% S. Pasquet - V22.05.04
% SUOP for matlab
% matop(infile,op,flip,outfile)

if nargin==3
    outfile=infile;
end

wsl = ispc_wsl;

infile_unix = unix_wsl_path(infile,wsl);
outfile_unix = unix_wsl_path(outfile,wsl);

dsppath=fileparts(infile);
dspfile=fullfile(dsppath,'tmp.dsp');
dspfile_unix = unix_wsl_path(dspfile,wsl);
if flip==0
    com1=sprintf('suflip < %s flip=1 > %s',infile_unix,dspfile_unix);
    unix_cmd(com1,wsl);
    com1=sprintf('suop < %s op=%s > %s',dspfile_unix,op,outfile_unix);
    unix_cmd(com1,wsl);
    com1=sprintf('suflip < %s flip=-1 > %s',outfile_unix,dspfile_unix);
    unix_cmd(com1,wsl);
else
    com1=sprintf('suop < %s op=%s > %s',outfile_unix,op,dspfile_unix);
    unix_cmd(com1,wsl);
end
movefile(fullfile(dsppath,'tmp.dsp'),outfile);
end