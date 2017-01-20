function matop(infile,op,flip,outfile)

%%% S. Pasquet - V17.01.13
% SUOP for matlab
% matop(infile,op,flip,outfile)

if nargin==3
    outfile=infile;
end
dsppath=fileparts(infile);
if flip==0
    com1=sprintf('suflip < %s flip=1 > %s',infile,fullfile(dsppath,'tmp.dsp'));
    unix(com1);
    com1=sprintf('suop < %s op=%s > %s',fullfile(dsppath,'tmp.dsp'),op,outfile);
    unix(com1);
    com1=sprintf('suflip < %s flip=-1 > %s',outfile,fullfile(dsppath,'tmp.dsp'));
    unix(com1);
else
    com1=sprintf('suop < %s op=%s > %s',outfile,op,fullfile(dsppath,'tmp.dsp'));
    unix(com1);
end
movefile(fullfile(dsppath,'tmp.dsp'),outfile);
end