function [X,Y,Z] = xyz2plot(X,Y,Z)

if min(size(X)) == 1
    Xnew = X;
    Xnew(1:end-1) = Xnew(1:end-1) - diff(Xnew)/2;
    Xnew(end) = Xnew(end) - diff(X(end-1:end))/2;
    Xnew(end+1) = Xnew(end) + diff(X(end-1:end));
    X = Xnew;
else
    ind_x = find(size(X) == length(unique(X)));
    if ind_x == 1
        Xnew = X';
    else
        Xnew = X;
    end
    Xnew(:,1:end-1) = Xnew(:,1:end-1) - diff(Xnew,1,ind_x)/2;
    Xnew(:,end) = Xnew(:,end) - diff(X(:,end-1:end),1,ind_x)/2;
    Xnew(:,end+1) = Xnew(:,end) + diff(X(:,end-1:end),1,ind_x);
    Xnew(end+1,:) = Xnew(end,:);
    if ind_x == 1
        X = Xnew';
    else
        X = Xnew;
    end
end
if min(size(Y)) == 1
    Ynew = Y;
    Ynew(1:end-1) = Ynew(1:end-1) - diff(Ynew)/2;
    Ynew(end) = Ynew(end) - diff(Y(end-1:end))/2;
    Ynew(end+1) = Ynew(end) + diff(Y(end-1:end));
    Y = Ynew;
else
    if ind_x == 1
        Ynew = Y';
        ind_y = 2;
    else
        Ynew = Y;
        ind_y = 1;
    end
    Ynew(1:end-1,:) = Ynew(1:end-1,:) - diff(Ynew,1,ind_y)/2;
    Ynew(end,:) = Ynew(end,:) - diff(Y(end-1:end,:),1,ind_y)/2;
    Ynew(end+1,:) = Ynew(end,:) + diff(Y(end-1:end,:),1,ind_y);
    Ynew(:,end+1) = Ynew(:,end);
    if ind_x == 1
        Y = Ynew';
    else
        Y = Ynew;
    end
end
Z = [Z Z(:,end)];
Z = [Z;Z(end,:)];