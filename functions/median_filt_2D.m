function [data_filt, mean_data_std] = median_filt_2D(data,m,n,test)

if nargin < 4
    test = 0;
end

[xlength, ylength] = size(data);
X = 1:xlength;
Y = 1:ylength;

if m > 1
    win_size_x = fix(m/2);
else
    win_size_x = 1;
end
if n > 1
    win_size_y = fix(n/2);
else
    win_size_y = 1;
end

data_filt = data;
data_std = data;

for i=1:xlength
    for j=1:ylength
        
        if isnan(data(i,j))
            continue
        end
        dat_select = data(X>i-win_size_x & X<i+win_size_x,Y>j-win_size_y & Y<j+win_size_y);
        data_std(i,j) = nanstd(dat_select(:));
        if test == 1
            if abs(median(dat_select(:))-data(i,j))>abs(median(dat_select(:))-std(dat_select(:)))
                data_filt(i,j) = nanmedian(dat_select(:));
            end
        else
            data_filt(i,j) = nanmedian(dat_select(:));
        end
    end
end
    mean_data_std = nanmean(data_std);