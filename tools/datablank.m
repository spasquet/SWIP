function datablank
%%% S. Pasquet - V16.2.17
% Blank X,Y,Z file and save blanked file

[filedata,pathdata]=uigetfile({'*.model;*.dat;*.xzv'},'Select XYZ data file');
[filebln,pathbln]=uigetfile({'*.bln'},'Select XYZ data file');

datafile=fullfile(pathdata,filedata);
blnfile=fullfile(pathbln,filebln);
[~,datablnfilename,extension]=fileparts(datafile);
datablnfile=fullfile(pathdata,[datablnfilename,'_bln',extension]);

[ZI,XI,YI]=bln2xyz(datafile,blnfile);

X=reshape(XI,1,size(XI,1)*size(XI,2));
Y=reshape(YI,1,size(XI,1)*size(XI,2));
Z=reshape(ZI,1,size(XI,1)*size(XI,2));

X=X(isnan(Z)==0)';
Y=Y(isnan(Z)==0)';
Z=Z(isnan(Z)==0)';

dlmwrite(datablnfile,[X,Y,Z],'delimiter',' ');

