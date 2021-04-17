
function [distance_val_SR,distances_val,D_sort,val_sort] = calculate_distance_from_tip_corr_ch12(tiploc,Ydata,ch1,ch2,histedge,maxdistlimit,sortoption)

BW=Ydata;
BW(~isnan(Ydata))=1;
ch1_seg=ch1.*BW;
ch2_seg=ch2.*BW;

ch1data=ch1(~isnan(ch1_seg));
ch2data=ch2(~isnan(ch2_seg));
ch12=[ch1data,ch2data];
mean_ch12=mean(ch12,1);


X=[tiploc(2,:);tiploc(1,:)]';

D_sort=[]; val_sort=[];

[r,c]= find(~isnan(Ydata));% >= histedge_Y(n) & Ydata< histedge_Y(n+1));
[indY]= find(~isnan(Ydata));% >= histedge_Y(n) & Y

if ~isequal(length(ch1data),length(ch2data),length(indY))
    length(ch1data)
    length(ch2data)
    length(indY)
    error('the size should be matched')
end


Y=[c,r]; %x, y;

D = pdist2(X,Y,'euclidean','Smallest',1); %% I is index of X

Val=Ydata(indY);
ch1val=ch1(indY);
ch2val=ch2(indY);

% remove distances are lareger that max distance limit
removeind=find(D>=maxdistlimit);
D(removeind)=[];
Val(removeind)=[];
ch1val(removeind)=[];
ch2val(removeind)=[];

% ordered by Y
D_val(:,1)=D;
D_val(:,2)=Val;
D_val(:,3)=ch1val;
D_val(:,4)=ch2val;


%order by Distance
D_val = sortrows(D_val,1) ;
D_sort=D_val(:,1);
val_sort=D_val(:,2);
ch1_sort=D_val(:,3);
ch2_sort=D_val(:,4);
%scatter (D_val_sort(:,1),D_val_sort(:,2))

clear D Val;
sortdata=[];
switch sortoption
    case 'distance'
        sortdata=D_sort;
    case 'value'
        sortdata=val_sort;
end

if isempty(sortdata)
    error('choose between distance or ratio')
end
%

for n=1:(length(histedge)-2)
    [indD]= find(sortdata >= histedge(n) &  sortdata< histedge(n+1));
    
    if isempty(indD)%size(r,1)>2
        D_hist=[];
        Val_hist=[];
        ch1_hist=[];
        ch2_hist=[];
        SR_hist=nan;
    else
        D_hist=D_sort(indD);
        Val_hist=val_sort(indD);
        ch1_hist=ch1_sort(indD);
        ch2_hist=ch2_sort(indD);
        s_ch1_hist=sqrt((sum((ch1_hist-mean(ch1data,1)).^2,1))/length(ch1_hist));
        s_ch2_hist=sqrt((sum((ch2_hist-mean(ch2data,1)).^2,1))/length(ch2_hist));
        
        Xhist=[ch1_hist, ch2_hist] ; %(part of X data)
        xc = Xhist - mean_ch12;
        c_hist = (xc' * xc) ./ size(Xhist,1);
        coeff_hist=c_hist/(s_ch1_hist*s_ch2_hist);
        SR_hist=coeff_hist(2);
        
    end
    
    
    distances_val{n,1}=D_hist;
    distances_val{n,2}=Val_hist;
    distances_val{n,3}=ch1_hist;
    distances_val{n,4}=ch2_hist;
    distance_val_SR(n)=SR_hist;
    indD=[]; coeff_hist=[];SR_hist=[];
end

n=(length(histedge)-1);

[indD]= find(sortdata>= histedge(n) & sortdata<= histedge(n+1));

if isempty(indD)
    D_hist=[];
    Val_hist=[];
    ch1_hist=[];
    ch2_hist=[];
    SR_hist=nan;
    
else
    D_hist=D_sort(indD);
    Val_hist=val_sort(indD);
    ch1_hist=ch1_sort(indD);
    ch2_hist=ch2_sort(indD);
    s_ch1_hist=sqrt((sum((ch1_hist-mean(ch1data,1)).^2,1))/length(ch1_hist));
    s_ch2_hist=sqrt((sum((ch2_hist-mean(ch2data,1)).^2,1))/length(ch2_hist));
    
    
    Xhist=[ch1_hist, ch2_hist] ;
    xc = Xhist - mean_ch12;
    c_hist = (xc' * xc) ./ size(Xhist,1);
    coeff_hist=c_hist/(s_ch1_hist*s_ch2_hist);
    SR_hist=coeff_hist(2);
    
    
end

distances_val{n,1}=D_hist;
distances_val{n,2}=Val_hist;
distances_val{n,3}=ch1_hist;
distances_val{n,4}=ch2_hist;
distance_val_SR(n)=SR_hist;
indD=[]; coeff_hist=[]; SR_hist=[];

end

