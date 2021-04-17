clear all; close all;

%-----set parameters -----%
ch1name='CD3-CF633';
ch2name='CD45-AF488';
removeisland_off=1;
Remove_Nucleus_option=0; % use if nuclues background signal is too high
pixelsize=50;
Exfactor=4; % gel-expansion factor
cellbodywidthpix=24;
closedistlimit=50;
MVlengthlimit_min_nm=50;
MVlengthlimit_max_nm=1500;
Rlim=500;
tiprangelengthnm =100; % nm scale
excEdgesizenm=75;
AUwl=650; % 650 nm emission peak for the AF633  % 600nm for emission peack for AF568

%% -----open a file -----%
[FileName, PathName] = uigetfile('*.mat','select a z-stack testfile',...
    'MultiSelect', 'on');
cd (PathName)

I=importdata(FileName);
nframe=size(I,1);
ExampleDisplayFrame=floor(nframe/2);

%% run analysis frame by frame
cnt=0;
for fr=1:5
    ch1=I{fr,1};
    ch2=I{fr,2};
    
    %----- get segmented area -----%
    level_ch1 = multithresh(ch1,2);
    multiThreshOut_ch1=level_ch1(1);
    BW1= imbinarize(ch1,multiThreshOut_ch1);
    level_ch2 = multithresh(ch2,2);
    multiThreshOut_ch2=level_ch2(1);
    BW2= imbinarize(ch2,multiThreshOut_ch2);
    BW3=BW1+BW2; BW3(BW3>0)=1;
    
    [SegBW,InBG,InEdgeBW,OutBG] = segmentation(imadd(ch1,ch2),BW3);
    level=graythresh(uint8(ch1));
    Thresh_ch1=level*255;
    BW1= imbinarize(ch1,Thresh_ch1);
    
    level=graythresh(uint8(ch2));
    Thresh_ch2=level*255;
    BW2= imbinarize(ch2,Thresh_ch2);
    BW3=BW1+BW2; BW3(BW3>0)=1;
    
    
    if Remove_Nucleus_option ==1
        [SegBW,InBG,InEdgeBW] = remove_Nucleus(OutBG,SegBW);
    end
    
    InEdgeBW_frame(:,:,fr)= InEdgeBW;
    
    [SegBW] = segmentation(imadd(ch1,ch2),BW3);
    SegBW=~InBG.*SegBW;
    
    ch1=ch1.*BW3;
    ch2=ch2.*BW3;
    ch1_seg=ch1.* SegBW;
    ch2_seg=ch2.* SegBW;
    
    se=strel('disk',cellbodywidthpix);
    BW=imdilate(InBG,se);
    BW=BW-InBG;
    CB=SegBW.*BW;
    MV=SegBW.*(~CB);
    SegL=InBG*2+CB*1+MV*3; 
    OutBG=SegL==0;

    %----- calculate RIS (normalized) -----%
    ch1_seg(SegBW==0)=nan;
    ch2_seg(SegBW==0)=nan;
    meanInt_ch1=mean(ch1_seg(:),'omitnan'); % mean intensity of ch1
    meanInt_ch2=mean(ch2_seg(:),'omitnan'); % mean intensity of ch2
    ch1_seg_norm=ch1_seg.*(meanInt_ch2/meanInt_ch1);
    RIS_seg_norm=(ch2_seg-ch1_seg_norm)./(ch1_seg_norm+ch2_seg);
    RIS_seg_norm_frame(:,:,fr)=RIS_seg_norm;
    
    %----- get the Microvilli-tip locations (tiploc) -----%
    
    if removeisland_off==1
        SegL_new=SegL;
    else
        [SegL_new] = remove_islands(RIS_seg_norm,SegL);
    end
    
    MVlengthlimit_min = (MVlengthlimit_min_nm*Exfactor/pixelsize)+cellbodywidthpix; % mimimun-length limit in pixel
    MVlengthlimit_max = (MVlengthlimit_max_nm*Exfactor/pixelsize)+cellbodywidthpix; % maximum-length limit in pixel
    [tiploc]=get_tiplocation(InEdgeBW,SegL_new,MVlengthlimit_min,MVlengthlimit_max,closedistlimit);
    
    
    %% Segmentation analysis
    if  ~isempty(tiploc)
        cnt=cnt+1;
        
        switch ch1name
            case 'L-sel-AF568'
                AUwl=600;
        end
        AU=ceil(0.61*AUwl/1.2/50/1.7);
        tiprangelength=tiprangelengthnm*Exfactor/pixelsize; % pixel (12*12.5)% 150nm
        excEdgesizepix=round(excEdgesizenm*Exfactor/pixelsize);
        
        %-------- segmentation and calculate Intensity for each segmented areas
        
        % SegBW(SegBW==0)=nan;
        [tipBW,colBW,cbBW,excBW,analresult] = segmentation_v2(tiploc,SegBW,SegL,ch1,ch2,Thresh_ch1,Thresh_ch2,tiprangelength,excEdgesizepix,AU);
        
        %-------- calculate SR
        
        ch1_tip=ch1(tipBW==1);
        ch2_tip=ch2(tipBW==1);
        [SR_tip]=calculate_SR(ch1,ch2,SegBW,ch1_tip,ch2_tip);
        
        
        ch1_col=ch1(colBW==1);
        ch2_col=ch2(colBW==1);
        [SR_col]=calculate_SR(ch1,ch2,SegBW,ch1_col,ch2_col);
        
        
        ch1_cb=ch1(cbBW==1);
        ch2_cb=ch2(cbBW==1);
        [SR_cb]=calculate_SR(ch1,ch2,SegBW,ch1_cb,ch2_cb);
        
        ch1_entire=ch1(SegBW==1);
        ch2_entire=ch2(SegBW==1);
        [SR_entire]=calculate_SR(ch1,ch2,SegBW,ch1_entire,ch2_entire);
        
        SR_seg.tip(cnt)=SR_tip; 
        SR_seg.col(cnt)=SR_col;
        SR_seg.cb(cnt)=SR_cb;
        SR_seg.entire(cnt)=SR_entire;
        
        %-------- calculate RIS for each segment
        colRIS_norm.mean(cnt)=mean(RIS_seg_norm(colBW==1));
        colRIS_norm.median(cnt)=median((RIS_seg_norm(colBW==1)));
        colRIS_norm.std(cnt)=std(RIS_seg_norm(colBW==1),[],1);
        colRIS_norm.perc25(cnt)=prctile((RIS_seg_norm(colBW==1)),25);
        colRIS_norm.perc75(cnt)=prctile((RIS_seg_norm(colBW==1)),75);
        
        tipRIS_norm.mean(cnt)=mean(RIS_seg_norm(tipBW==1));
        tipRIS_norm.median(cnt)=median((RIS_seg_norm(tipBW==1)));
        tipRIS_norm.std(cnt)=std(RIS_seg_norm(tipBW==1),[],1);
        tipRIS_norm.perc25(cnt)=prctile((RIS_seg_norm(tipBW==1)),25);
        tipRIS_norm.perc75(cnt)=prctile((RIS_seg_norm(tipBW==1)),75);
        
        cbRIS_norm.mean(cnt)=mean(RIS_seg_norm(cbBW==1));
        cbRIS_norm.median(cnt)=median((RIS_seg_norm(cbBW==1)));
        cbRIS_norm.std(cnt)=std(RIS_seg_norm(cbBW==1),[],1);
        cbRIS_norm.perc25(cnt)=prctile((RIS_seg_norm(cbBW==1)),25);
        cbRIS_norm.perc75(cnt)=prctile((RIS_seg_norm(cbBW==1)),75);
        
        entireRIS_norm.mean(cnt)=mean(RIS_seg_norm(SegBW==1));
        entireRIS_norm.median(cnt)=median((RIS_seg_norm(SegBW==1)));
        entireRIS_norm.std(cnt)=std(RIS_seg_norm(SegBW==1),[],1);
        entireRIS_norm.perc25(cnt)=prctile((RIS_seg_norm(SegBW==1)),25);
        entireRIS_norm.perc75(cnt)=prctile((RIS_seg_norm(SegBW==1)),75);
        
        %--------
        Segmentation.colBW_All{cnt}=colBW;
        Segmentation.tipBW_All{cnt}=tipBW;
        Segmentation.cbBW_All{cnt}=cbBW;
        
        Segmentation.ch1_tip_area{cnt}=analresult.ch1_tip.ch1_tip_area;
        Segmentation.ch2_tip_area{cnt}=analresult.ch2_tip.ch2_tip_area;
        Segmentation.ch1_col_area{cnt}=analresult.ch1_col.ch1_col_area;
        Segmentation.ch2_col_area{cnt}=analresult.ch2_col.ch2_col_area;
        Segmentation.ch1_cb_area{cnt}=analresult.ch1_cb.ch1_cb_area;
        Segmentation.ch2_cb_area{cnt}=analresult.ch2_cb.ch2_cb_area;
        
        Segmentation.statresult.ch1_tip_mean_frame(cnt)=analresult.ch1_tip.ch1_tiparea_mean;
        Segmentation.statresult.ch1_tip_median_frame(cnt)=analresult.ch1_tip.ch1_tiparea_median;
        Segmentation.statresult.ch2_tip_mean_frame(cnt)=analresult.ch2_tip.ch2_tiparea_mean;
        Segmentation.statresult.ch2_tip_median_frame(cnt)=analresult.ch2_tip.ch2_tiparea_median;
        
        Segmentation.statresult.ch1_col_mean_frame(cnt)=analresult.ch1_col.ch1_colarea_mean;
        Segmentation.statresult.ch1_col_median_frame(cnt)=analresult.ch1_col.ch1_colarea_median;
        Segmentation.statresult.ch2_col_mean_frame(cnt)=analresult.ch2_col.ch2_colarea_mean;
        Segmentation.statresult.ch2_col_median_frame(cnt)=analresult.ch2_col.ch2_colarea_median;
        
        Segmentation.statresult.ch1_cb_mean_frame(cnt)=analresult.ch1_cb.ch1_cbarea_mean;
        Segmentation.statresult.ch1_cb_median_frame(cnt)=analresult.ch1_cb.ch1_cbarea_median;
        Segmentation.statresult.ch2_cb_mean_frame(cnt)=analresult.ch2_cb.ch2_cbarea_mean;
        Segmentation.statresult.ch2_cb_median_frame(cnt)=analresult.ch2_cb.ch2_cbarea_median;
        
        NewSegL=zeros(size(SegL));
        NewSegL(cbBW==1)= 1;
        NewSegL(colBW==1)=2;
        NewSegL(tipBW==1)=3;
        NewSegL(InEdgeBW==1)=4;
        NewSegL(excBW==1)=5;
        NewSegL_all(:,:,cnt)=NewSegL;
    end
    
    %-------- calculate 'Intensity', 'RIS', 'SR'' from tips of MVs
    if  ~isempty(tiploc)
        max_radius=Rlim*Exfactor/pixelsize;
        radiushistedge=0:1:max_radius; 
        [SR_distance_tip,distance_RIS_tip,D_tip_sort,RIS_tip_sort] = calculate_distance_from_tip_corr_ch12(tiploc,RIS_seg_norm,ch1,ch2,radiushistedge,max_radius,'distance');
        
        distance_tip_all(fr,1:length(radiushistedge)-1)=distance_RIS_tip(:,1);
        RIS_tip_all(fr,1:length(radiushistedge)-1)=distance_RIS_tip(:,2);
        ch1_tip_all(fr,1:length(radiushistedge)-1)=distance_RIS_tip(:,3);
        ch2_tip_all(fr,1:length(radiushistedge)-1)=distance_RIS_tip(:,4);
        SR_distance_tip_all(fr,:)=SR_distance_tip;
    end
    
    %---- subtract background intensity for each channel (for intensity plot)
    for n=1:size(distance_RIS_tip,1)
        ch1_subbac=distance_RIS_tip{n,3}-Thresh_ch1;
        ch2_subbac=distance_RIS_tip{n,4}-Thresh_ch2;
        ch1_subbac(ch1_subbac<0)=0;
        ch2_subbac(ch2_subbac<0)=0;
        ch1_tip_all_subbac{fr,n}=ch1_subbac;
        ch2_tip_all_subbac{fr,n}=ch2_subbac;
        
    end
    %-----  draw figures for each frame
    figure (1)
    set(gcf,'Position',get(0,'Screensize'));
    
    minmax_ch1=[min(ch1(:)),max(ch1(:))];
    minmax_ch2=[min(ch2(:)),max(ch2(:))];
    ch1_I = mat2gray(ch1,minmax_ch1);
    ch2_I = mat2gray(ch2,minmax_ch2);
    lowhigh_ch1=stretchlim(ch1_I,[0.01,0.995]);
    lowhigh_ch2=stretchlim(ch2_I,[0.01,0.995]);
    
    Thtable{1}(1)=round(minmax_ch1(1)*lowhigh_ch1(1));
    Thtable{1}(2)=round(minmax_ch1(2)*lowhigh_ch1(2));
    Thtable{2}(1)=round(minmax_ch2(1)*lowhigh_ch2(1));
    Thtable{2}(2)=round(minmax_ch2(2)*lowhigh_ch2(2));
    
    chcolormap{1}= scaledcolormap((Thtable{1}(2)-Thtable{1}(1)),[1,0,1],0,1); %red
    chcolormap{2}= scaledcolormap((Thtable{2}(2)-Thtable{2}(1)),[0,1,0],0,1); %green
    rgbImage_ch1 = ind2rgb(ch1, chcolormap{1});
    rgbImage_ch2 = ind2rgb(ch2, chcolormap{2});
    rgbImage_mer = imadd(rgbImage_ch1,rgbImage_ch2) ;
    
    %----- < subplot 1: ch1 image > -----%
    subplot(2,3,1)
    imagesc(rgbImage_ch1); axis image
    title(ch1name)
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'box','on');
    set(gca,'YColor','k');
    set(gca,'XColor','k');
    
    %----- < subplot 2: ch2 image > -----%
    subplot(2,3,2)
    
    imagesc(rgbImage_ch2); axis image
    title(ch2name)
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'box','on');
    set(gca,'YColor','k');
    set(gca,'XColor','k');
    
    %----- < subplot 3: Merged image of ch1 and ch2> -----%
    subplot(2,3,3);
    
    imagesc(rgbImage_mer);
    axis image
    title('Merge');
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'box','on');
    set(gca,'YColor','k');
    set(gca,'XColor','k');
    
    
    %----- < subplot 4: segmented map> -----%
    subplot(2,3,5);
    mycolormap = [0.2422 0.1504 0.6603; 0.1540 0.5902 0.9218; 0.5 0.8 0.4; 0.9769 0.9839 0.0805;1,1,1;0.5 0.5 0.5];
    
    colormap(mycolormap)
    imagesc(NewSegL);
    axis image tight
    hold on
    caxis([0,5]);
    
    x0 = get(gca,'xlim') ;
    y0 = get(gca,'ylim') ;
    % draw dummy data to show legend
    plot(0,0,'Color',[1,0,0]);%,'MarkerFaceColor',[0.95,0.95,0.95])
    hold on
    plot(0,0,'Color',[0.9769    0.9839    0.0805]); % MV-tip-Area
    hold on
    plot(0,0,'Color',[0.5000    0.8000    0.4000]); % MV-col-Area
    hold on
    plot(0,0,'Color',[0.1540    0.5902    0.9218]); %CB
    hold on
    plot(0,0,'Color',[0.5    0.5    0.5]); % excluded area
    hold on
    plot(0,0,'Color',[1,1,1]); % Inbackground Edge-line
    axis([x0 y0])
    axis off %hide axis
    freezeColors
    [lgd, hobj] = legend('MV-tips(*)','MV-tip-Area','MV-col-Area','CB','excluded','InBG Edge-line');
    lgd.Position=[0.6108    0.20    0.12    0.19];
    lgd.Box='off';
    lgd.FontSize= 11;
    set(hobj,'linewidth',3);
    hold on
    if ~isempty(tiploc)
        scatter(tiploc(2,:),tiploc(1,:),'*','r')
    end
    
    title('segmented map')
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'box','on'); 
    set(gca,'YColor','k');
    set(gca,'XColor','k');
    axis image tight
    hold off
    
    %----- < subplot 5: RIS image> -----%
    subplot(2,3,4);
    imagesc(RIS_seg_norm);
    axis image tight
    flipedColormap  =flipud(colormap(parula(25)));
    colormap(flipedColormap);
    title('RIS');
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'box','on');
    set(gca,'YColor','k');
    set(gca,'XColor','k');
    colorbar
    c = colorbar;
    c.XTick=-1:0.25:1;
    caxis([-1,1]);
    freezeColors
    hold off
    add_scalebar(2000,[1,1,1],3,pixelsize,Exfactor,4);% draw scale bar
    % get frame for movie write
    if nframe>1
        F = getframe(gcf);
        Anal_image(:,:,:,fr)=F.cdata;
    end
end

%----- write movie file
if nframe>1
    
    mov = immovie(Anal_image);
    savefilename=['frames_' FileName(1:end-4)];
    v = VideoWriter([savefilename,'_RIS'],'Motion JPEG AVI'); 
    v.FrameRate=2;
    v.Quality=100;
    open(v)
    writeVideo(v,mov)
    close(v)
    close all
    clear  mov
end
%% ---combine data from all frames

statresult.ch1_tip_all_hist=hist_combine(ch1_tip_all_subbac','distance',pixelsize,Exfactor);
statresult.ch1_tip_all=ch1_tip_all_subbac;
[comb_distance_tip_ch1]=combine_multicells(statresult.ch1_tip_all);
statresult.comb_distance_tip_ch1=hist_combine(comb_distance_tip_ch1,'distance',pixelsize,Exfactor);
statresult.ch2_tip_all_hist=hist_combine(ch2_tip_all_subbac','distance',pixelsize,Exfactor);
statresult.ch2_tip_all=ch2_tip_all_subbac;
[comb_distance_tip_ch2]=combine_multicells(statresult.ch2_tip_all);
statresult.comb_distance_tip_ch2=hist_combine(comb_distance_tip_ch2,'distance',pixelsize,Exfactor);
statresult.distance_tip_hist_in_frame=hist_combine(distance_tip_all','distance',pixelsize,Exfactor);
statresult.RIS_from_tip_hist_in_frame=hist_combine(RIS_tip_all','distance',pixelsize,Exfactor);
statresult.distance_tip_all=distance_tip_all;
statresult.RIS_tip_all=RIS_tip_all;

[comb_RIS_tip]=combine_multicells(statresult.RIS_tip_all);
[comb_distance_tip]=combine_multicells(statresult.distance_tip_all);

statresult.comb_RIS_from_MVtip=hist_combine(comb_RIS_tip,'distance',pixelsize,Exfactor);
statresult.comb_distance_from_MVtip=hist_combine(comb_distance_tip,'distance',pixelsize,Exfactor);
statresult.comb_SR_from_MVtip = stat_combine_mat(SR_distance_tip_all,1,0,1);

Segmentation.tipRIS_norm=tipRIS_norm;
Segmentation.colRIS_norm=colRIS_norm;
Segmentation.cbRIS_norm=cbRIS_norm;
Segmentation.entireRIS_norm=entireRIS_norm;

[Segmentation.statresult.ch1_tiparea] = stat_combine_cell(Segmentation.ch1_tip_area,1);
[Segmentation.statresult.ch1_colarea] = stat_combine_cell(Segmentation.ch1_col_area,1);
[Segmentation.statresult.ch1_cbarea] = stat_combine_cell(Segmentation.ch1_cb_area,1);

[Segmentation.statresult.ch2_tiparea] = stat_combine_cell(Segmentation.ch2_tip_area,1);
[Segmentation.statresult.ch2_colarea] = stat_combine_cell(Segmentation.ch2_col_area,1);
[Segmentation.statresult.ch2_cbarea] = stat_combine_cell(Segmentation.ch2_cb_area,1);

Segmentation.SR_seg=SR_seg;
Segmentation.SR_seg.stat_tip=stat_combine_mat(Segmentation.SR_seg.tip,2,0,1);
Segmentation.SR_seg.stat_col=stat_combine_mat(Segmentation.SR_seg.col,2,0,1);
Segmentation.SR_seg.stat_cb=stat_combine_mat(Segmentation.SR_seg.cb,2,0,1);
Segmentation.SR_seg.stat_entire=stat_combine_mat(Segmentation.SR_seg.entire,2,0,1);

ch1_tip_mean_frame=mean(Segmentation.statresult.ch1_tip_median_frame);
ch2_tip_mean_frame=mean(Segmentation.statresult.ch2_tip_median_frame);
ch1_col_mean_frame=mean(Segmentation.statresult.ch1_col_median_frame);
ch2_col_mean_frame=mean(Segmentation.statresult.ch2_col_median_frame);
ch1_cb_mean_frame=mean(Segmentation.statresult.ch1_cb_median_frame);
ch2_cb_mean_frame=mean(Segmentation.statresult.ch2_cb_median_frame);

ch1_tip_std_frame=std(Segmentation.statresult.ch1_tip_median_frame,[],2);
ch2_tip_std_frame=std(Segmentation.statresult.ch2_tip_median_frame,[],2);
ch1_col_std_frame=std(Segmentation.statresult.ch1_col_median_frame,[],2);
ch2_col_std_frame=std(Segmentation.statresult.ch2_col_median_frame,[],2);
ch1_cb_std_frame=std(Segmentation.statresult.ch1_cb_median_frame,[],2);
ch2_cb_std_frame=std(Segmentation.statresult.ch2_cb_median_frame,[],2);

RIS_norm_tip_mean_frame=mean(Segmentation.tipRIS_norm.mean);
RIS_norm_col_mean_frame=mean(Segmentation.colRIS_norm.mean);
RIS_norm_cb_mean_frame=mean(Segmentation.cbRIS_norm.mean);
RIS_norm_entire_mean_frame=mean(Segmentation.entireRIS_norm.mean);

RIS_norm_tip_std_frame=std(Segmentation.tipRIS_norm.mean,[],2);
RIS_norm_col_std_frame=std(Segmentation.colRIS_norm.mean,[],2);
RIS_norm_cb_std_frame=std(Segmentation.cbRIS_norm.mean,[],2);
RIS_norm_entire_std_frame=std(Segmentation.entireRIS_norm.mean,[],2);

RIS_norm_tip_median_frame=mean(Segmentation.tipRIS_norm.median);
RIS_norm_col_median_frame=mean(Segmentation.colRIS_norm.median);
RIS_norm_cb_median_frame=mean(Segmentation.cbRIS_norm.median);
RIS_norm_entire_median_frame=mean(Segmentation.entireRIS_norm.median);

RIS_norm_tip_median_std_frame=std(Segmentation.tipRIS_norm.median,[],2);
RIS_norm_col_median_std_frame=std(Segmentation.colRIS_norm.median,[],2);
RIS_norm_cb_median_std_frame=std(Segmentation.cbRIS_norm.median,[],2);
RIS_norm_entire_median_std_frame=std(Segmentation.entireRIS_norm.median,[],2);

%% figure for intensity,
figure(1)
set(gcf, 'Position',  [301.8000  50.6000  979.2000  500.8000])

% -------- ch1 intensity
subplot(1,4,1)

medianbar_data = [Segmentation.statresult.ch1_tiparea.median, Segmentation.statresult.ch1_colarea.median, Segmentation.statresult.ch1_cbarea.median];
medianbar_err= [Segmentation.statresult.ch1_tiparea.ste,Segmentation.statresult.ch1_colarea.ste,Segmentation.statresult.ch1_cbarea.ste];

c = categorical({'MV-tip area','MV-column area','Cellbody'});
b= bar(c,medianbar_data);
b.FaceColor = 'flat';
b.CData(1,:) = mycolormap(2,:);
b.CData(2,:) = mycolormap(3,:);
b.CData(3,:) = mycolormap(4,:);
title([ch1name] ,'FontSize',14, 'Interpreter','none')
ylabel('Intensity (median)')
hold on
er = errorbar(c,medianbar_data,medianbar_err);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
set(gca,'FontSize',14);

% -------- ch2 intensity
subplot(1,4,2)

medianbar_data2 = [Segmentation.statresult.ch2_tiparea.median, Segmentation.statresult.ch2_colarea.median, Segmentation.statresult.ch2_cbarea.median];
medianbar_err2 = [Segmentation.statresult.ch2_tiparea.ste,Segmentation.statresult.ch2_colarea.ste,Segmentation.statresult.ch2_cbarea.ste];

c = categorical({'MV-tip area','MV-column area','Cellbody'});
b= bar(c,medianbar_data2);
b.FaceColor = 'flat';

b.CData(1,:) = mycolormap(2,:);
b.CData(2,:) = mycolormap(3,:);
b.CData(3,:) = mycolormap(4,:);

title([ch2name] ,'FontSize',14, 'Interpreter','none')
ylabel('Intensity (median)')
hold on

er = errorbar(c,medianbar_data2,medianbar_err2);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
set(gca,'FontSize',14);

% -------- RIS
subplot(1,4,3)
bardata_RIS_norm_median=[RIS_norm_entire_median_frame,RIS_norm_tip_median_frame,RIS_norm_col_median_frame,RIS_norm_cb_median_frame];
barerr_RIS_norm_median=[RIS_norm_entire_median_std_frame,RIS_norm_tip_median_std_frame,RIS_norm_col_median_std_frame,RIS_norm_cb_median_std_frame];
c = categorical({'Cell','MV-tip area','MV-column area','Cellbody'});
b=bar(c,bardata_RIS_norm_median);
b.FaceColor = 'flat';
b.CData(1,:) = mycolormap(1,:);
b.CData(2,:) = mycolormap(2,:);
b.CData(3,:) = mycolormap(3,:);
b.CData(4,:) = mycolormap(4,:);
title('RIS' ,'FontSize',14, 'Interpreter','none')
ylabel('RIS','FontSize',14)
hold on
er = errorbar(c,bardata_RIS_norm_median,barerr_RIS_norm_median);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
set(gca,'FontSize',14);
%box on

% -------- SR

subplot(1,4,4)
c = categorical({'Cell','MV-tip area','MV-column area','Cellbody'});
SRbar= [Segmentation.SR_seg.stat_entire.mean, Segmentation.SR_seg.stat_tip.mean, Segmentation.SR_seg.stat_col.mean,Segmentation.SR_seg.stat_cb.mean];
SRbar_err= [Segmentation.SR_seg.stat_entire.ste, Segmentation.SR_seg.stat_tip.ste, Segmentation.SR_seg.stat_col.ste,Segmentation.SR_seg.stat_cb.ste];
b= bar(c,SRbar);
b.FaceColor = 'flat';

b.CData(1,:) = mycolormap(1,:);
b.CData(2,:) = mycolormap(2,:);
b.CData(3,:) = mycolormap(3,:);
b.CData(4,:) = mycolormap(4,:);
title(['SR'''] ,'FontSize',14, 'Interpreter','none')
ylabel('SR''')
hold on
hold on
er = errorbar(c,SRbar,SRbar_err);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
set(gca,'FontSize',14);
savefig(gcf,'result_plot1')

%% Summary Figure for Intensity, RIS and SR'
figure(2)
set(gcf, 'Position',  [301.8000  50.6000  979.2000  500.8000])

%----  <subplot1: ch1 intensity plot  > ----%
subplot(1,4,1);
histedgelengthnm=radiushistedge*pixelsize/Exfactor;
p1=errorbar(histedgelengthnm(1:end-1),statresult.comb_distance_tip_ch1.median,statresult.comb_distance_tip_ch1.ste,'Marker','o');
xlim([0,histedgelengthnm(end)])
xlabel('distance from tips (nm)')
ylabel('Intensity')
title([ch1name,'(median)'])

%----  <subplot1: ch2 intensity plot > ----%
subplot(1,4,2);
p1=errorbar(histedgelengthnm(1:end-1),statresult.comb_distance_tip_ch2.median,statresult.comb_distance_tip_ch2.ste,'Marker','o');
xlim([0,histedgelengthnm(end)])
xlabel('distance from tips (nm)')
ylabel(['Intensity'])
title([ch2name,'(median)'])

%----  <subplot1: RIS plot > ----%
subplot(1,4,3);
p1=errorbar(histedgelengthnm(1:end-1),statresult.comb_RIS_from_MVtip.median,statresult.comb_RIS_from_MVtip.ste,'Marker','o');
hline = refline([0 0]);
hline.Color=[0.5,0.5,0.5];
ylim([-1,1])
xlim([0,histedgelengthnm(end)])
xlabel('distance from tips (nm)')
ylabel('RIS')
title('RIS (median)')
hold off

%----  <subplot4: SR' plot > ----%
subplot(1,4,4);
p1=errorbar(histedgelengthnm(1:end-1),statresult.comb_SR_from_MVtip.mean,statresult.comb_SR_from_MVtip.ste,'Marker','o');
hline = refline([0 0]);
hline.Color=[0.5,0.5,0.5];
xlabel('distance from tips (nm)')
ylabel('SR''')
title('SR''')
savefig(gcf,'result_plot2')

