function meanvel=meanvel_find(dspmat_weight_win,v)

meanvel=zeros(size(dspmat_weight_win,3),1);
for i=1:min(size(dspmat_weight_win,3))
    dspmat_temp=dspmat_weight_win(:,:,i);
    [~,I]=max(dspmat_temp');
    meanvel(i)=mean(v(I));
end