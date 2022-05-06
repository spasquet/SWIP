function model=readmodel(file,nmaxmod,runnb)

%%% S. Pasquet - V20.03.23
% Read velocity and density models calculated by dinver inversion and
% extracted using gpdcreport
% model=readmodel(file,nmaxmod,runnb)
% Handles dinver > 2

version = geopsy_version();

% Read model
fileID=fopen(file);
formatSpec = '%f %f';
fgets(fileID); % Read line
nn=0;
model=cell(nmaxmod,2);
for n=1+(runnb-1)*nmaxmod:(runnb-1)*nmaxmod+nmaxmod
    nn=nn+1;
    if version == 2
        fgets(fileID); % Read line
        model(nn,1:2)=textscan(fileID,formatSpec,'HeaderLines',1);
    else
        model(nn,1:2)=textscan(fileID,formatSpec,'HeaderLines',1);
    end
end
fclose(fileID);