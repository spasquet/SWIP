function matshift(infile,tshift,outfile)

%%% S. Pasquet - V23.02.22
% Shift time in SU files
% matshift(infile,tshift,outfile)
% tshift in ms (ex: -100 to shift upward by 100ms)

if nargin==2
    outfile=infile;
end
supath=fileparts(infile);
com1=sprintf('sushw < %s a=%s,0 key=delrt,tstat | sustatic hdrs=1 sign=1 > %s',infile,num2str(tshift),fullfile(supath,'tmp.su'));
unix(com1);

movefile(fullfile(supath,'tmp.su'),outfile);
end