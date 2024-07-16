function c = makegrid(x,z)

%function to make a grid of x and z points or x and y whichever you prefer.
% returns a two by (x*z) matrix of x,z pairs


[mx,nx]=size(x);
[mz,nz]=size(z);

if mx < nx;
    x=x';
end

if mz < nz;
    z=z';
end

tmpx=[];
tmpz=[];

for i=1:length(z);
    tmpx=[tmpx; x];
end
    
for i=1:length(z);
    tmpz=[tmpz; z(i)*ones(length(x),1)];
end

tmpc=[tmpx tmpz];

tmpc=sortrows(tmpc,2);
tmpc=sortrows(tmpc,1);


c=tmpc;