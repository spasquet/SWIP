function model=readmodel(file,nmaxmod,runnb)

% S. Pasquet - V16.1.21
% Read velocity and density models calculated by dinver inversion and
% extracted using gpdcreport

% Read model
fileID=fopen(file);
formatSpec = '%f %f';
fgets(fileID); % Read line
nn=0;
model=cell(nmaxmod,2);
for n=1+(runnb-1)*nmaxmod:(runnb-1)*nmaxmod+nmaxmod
    nn=nn+1;
    fgets(fileID); % Read line
    model(nn,1:2)=textscan(fileID,formatSpec,'HeaderLines',1);
end
fclose(fileID);