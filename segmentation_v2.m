function [tipBW,colBW,cbBW,excBW,analresult] = segmentation_v2(tiploc,SegBW,SegL,ch1,ch2,thresh_ch1,thresh_ch2,tiprangelength,excEdgesize,AU)
MVBW=zeros(size(SegL));
cbBW=zeros(size(SegL));
MVBW(SegL==3)=1; 
cbBW(SegL==1)=1;

ch1MV_seg=ch1.*MVBW;
ch2MV_seg=ch2.*MVBW;

se = strel('disk',12);
DilMVBW=imdilate(MVBW,se);

cbBW_new=cbBW.*(~DilMVBW);
MVBW_new=(cbBW-cbBW_new)+MVBW;

X=[tiploc(2,:);tiploc(1,:)]';

[r,c]= find(~isnan(MVBW));  

Y=[c,r]; %x, y;
D = pdist2(X,Y,'euclidean','Smallest',1); %% I is index of X not Y;



%% tip vs column area BW map
tipBW=zeros(size(MVBW_new));
tipBW(D<=tiprangelength+AU)=1;
tipBW=MVBW_new.*tipBW;


se = strel('disk',excEdgesize+AU);
excBW= (imdilate(tipBW,se)-tipBW).*MVBW_new;

excBW3=imdilate(MVBW_new,se).*cbBW_new;
excBW=double(imbinarize(excBW+excBW3));

cbBW_new=double(imbinarize(cbBW_new-excBW3));



colBW=zeros(size(MVBW_new));
colBW(D>tiprangelength+AU)=1;
colBW=double(imbinarize(colBW-excBW));

colBW=MVBW_new.*colBW;
cbBW=cbBW_new;

%%
tipBW(tipBW==0)=nan;% change 0 to nan
colBW(colBW==0)=nan;
cbBW(cbBW==0)=nan;

tip_area=sum(tipBW(:),'omitnan');
col_area=sum(colBW(:),'omitnan');
cb_area =sum(cbBW(:),'omitnan');

ch1_bgsb=ch1;
ch2_bgsb=ch2;
ch1_bgsb=ch1_bgsb-thresh_ch1;
ch2_bgsb=ch2_bgsb-thresh_ch2;

ch1_bgsb(ch1_bgsb<0)=0;
ch2_bgsb(ch2_bgsb<0)=0;

minmax_ch1=[min(ch1_bgsb(:)),max(ch1_bgsb(:))];
minmax_ch2=[min(ch2_bgsb(:)),max(ch2_bgsb(:))];

ch1_I = mat2gray(ch1_bgsb,minmax_ch1);
ch2_I= mat2gray(ch2_bgsb,minmax_ch2);
ch1_lowhigh=stretchlim(ch1_I,[0,0.99]); % [0.01 0.99], saturating the upper 1% and the lower 1%.
ch2_lowhigh=stretchlim(ch2_I,[0,0.99]);
ch1_norm = imadjust(ch1_I,[ch1_lowhigh(1),ch1_lowhigh(2)]);
ch2_norm = imadjust(ch2_I,[ch2_lowhigh(1),ch2_lowhigh(2)]);


ch1_norm=ch1_norm.*SegBW;
ch2_norm=ch2_norm.*SegBW;


ch1_tip_area = tipBW.*ch1_norm;
ch2_tip_area = tipBW.*ch2_norm;

ch1_col_area = colBW.*ch1_norm;
ch2_col_area = colBW.*ch2_norm;

ch1_cb_area = cbBW.*ch1_norm;
ch2_cb_area = cbBW.*ch2_norm;

%% result

ch1_tip_stat.ch1_tip_area=ch1_tip_area(~isnan(ch1_tip_area));
ch1_tip_stat.ch1_tiparea_mean=mean(ch1_tip_area(:),'omitnan');
ch1_tip_stat.ch1_tiparea_median=median(ch1_tip_area(:),'omitnan');
ch1_tip_stat.ch1_tiparea_std=std(ch1_tip_area(:),[],'omitnan');
ch1_tip_stat.ch1_tiparea_ste=ch1_tip_stat.ch1_tiparea_std/(tip_area^2);

ch2_tip_stat.ch2_tip_area=ch2_tip_area(~isnan(ch2_tip_area));
ch2_tip_stat.ch2_tiparea_mean=mean(ch2_tip_area(:),'omitnan');
ch2_tip_stat.ch2_tiparea_median=median(ch2_tip_area(:),'omitnan');
ch2_tip_stat.ch2_tiparea_std=std(ch2_tip_area(:),[],'omitnan');
ch2_tip_stat.ch2_tiparea_ste=ch2_tip_stat.ch2_tiparea_std/(tip_area^2);


ch1_col_stat.ch1_col_area=ch1_col_area(~isnan(ch1_col_area));
ch1_col_stat.ch1_colarea_mean=mean(ch1_col_area(:),'omitnan');
ch1_col_stat.ch1_colarea_median=median(ch1_col_area(:),'omitnan');
ch1_col_stat.ch1_colarea_std=std(ch1_col_area(:),[],'omitnan');
ch1_col_stat.ch1_colarea_ste=ch1_col_stat.ch1_colarea_std/(col_area^2);

ch2_col_stat.ch2_col_area=ch2_col_area(~isnan(ch2_col_area));
ch2_col_stat.ch2_colarea_mean=mean(ch2_col_area(:),'omitnan');
ch2_col_stat.ch2_colarea_median=median(ch2_col_area(:),'omitnan');
ch2_col_stat.ch2_colarea_std=std(ch2_col_area(:),[],'omitnan');
ch2_col_stat.ch2_colarea_ste=ch2_col_stat.ch2_colarea_std/(col_area^2);


ch1_cb_stat.ch1_cb_area=ch1_cb_area(~isnan(ch1_cb_area));
ch1_cb_stat.ch1_cbarea_mean=mean(ch1_cb_area(:),'omitnan');
ch1_cb_stat.ch1_cbarea_median=median(ch1_cb_area(:),'omitnan');
ch1_cb_stat.ch1_cbarea_std=std(ch1_cb_area(:),[],'omitnan');
ch1_cb_stat.ch1_cbarea_ste=ch1_cb_stat.ch1_cbarea_std/(cb_area^2);

ch2_cb_stat.ch2_cb_area=ch2_cb_area(~isnan(ch2_cb_area));
ch2_cb_stat.ch2_cbarea_mean=mean(ch2_cb_area(:),'omitnan');
ch2_cb_stat.ch2_cbarea_median=median(ch2_cb_area(:),'omitnan');
ch2_cb_stat.ch2_cbarea_std=std(ch2_cb_area(:),[],'omitnan');
ch2_cb_stat.ch2_cbarea_ste=ch2_cb_stat.ch2_cbarea_std/(cb_area^2);

analresult.ch1_tip=ch1_tip_stat;
analresult.ch2_tip=ch2_tip_stat;
analresult.ch1_col=ch1_col_stat;
analresult.ch2_col=ch2_col_stat;
analresult.ch1_cb=ch1_cb_stat;
analresult.ch2_cb=ch2_cb_stat;
end

