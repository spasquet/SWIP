function matop(infile,op,flip,outfile)

%%% S. Pasquet - V16.4.14
% SUOP for matlab

if nargin==3
    outfile=infile;
end
if flip==0
    com1=sprintf('suflip < %s flip=1 > tmp.dsp',infile);
    unix(com1);
    com1=sprintf('suop < tmp.dsp op=%s > %s',op,outfile);
    unix(com1);
    com1=sprintf('suflip < %s flip=-1 > tmp.dsp',outfile);
    unix(com1);
else
    com1=sprintf('suop < %s op=%s > tmp.dsp',outfile,op);
    unix(com1);
end
movefile('tmp.dsp',outfile);
end