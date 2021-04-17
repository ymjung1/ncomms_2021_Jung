%%% for full ROI image anlysis : ROI setting: 'X (pix)':'1,427','Y (pix)': '1,684' (defalut)
%%% for single cell: choose specific x, y range; input 'X (pix)':'130,210','Y (pix)': '180,260' 
%%% choose ch0, choose ch1;
clc; clear all;
close all

%% channel and threshold parameter setting
ch=2; % ch0=0; ch1=1, both=2
SetThreshold=1; % yes=1, no=0
displayplot=0;% yes=1, no=0 data regression plot for thresholding steps
%% figure parameter setting
drawFig=1; % yes=1, no=0
writepixsizeimage=0;
C_range=[-400 400 ]; % color axis range for z-caxis; []==> [min max]; %C_range=[100,200];;
dotsize=1;% marker size on image
%% statistic analysis
do_stat=1; % yes=1, no=0
%% merge both channel
merge_ch=1;% merge two channel images
displaychange=0; % show before and after the shift correction
%% shift correction parametter setting
do_driftcorrect=1;% yes=1, no=0
manual_driftcorrect=0;
shiftRef='STORMPixImage';%'LocImage'; %'STORMPixImage';
shiftfindoption=2; % 1=weightedcentroid 2=2D gaussian fit
ws=11;% for xcoor-window
filtertype='GaussFilter';% filter type ['MeanFilter', 'MedianFilter','LogFilter','GaussFilter']
showpeakchoice=0; % display on the peak of the 2D cross correlation image
ccresultsave=1; % display xcoor2 result yes=1, no=1
ApplyfilterCC=0; % cmean filter
applyfilteronimage=0; %% apply blur filter sigma default 1;
kernelsizeCC=7; % X-corr filter kernelsize
filtersigma=2; % X-coor filter sigmasize
pixelsize =117; % ONI 100x objective, fixed value
scalebarsize=5000; % nm scale
% annotation =1; % scalebar annotation display, yes =1, no=0;]
%% Rendering Image parameter setting %%
Supres=1; % yes=1, no=0 ; rendering image in nm scale
RendSize=10; % if Supres==1, input desired pixel-size (nm)
RendGaussfilterSigma=1; %% Render Gaussian filtler
drawSTORMpiximg=1;

%%  Segmentation and 3D-Pair-correlation analysis parameter setting 
MVlengthlimit_min_nm=50;
MVlengthlimit_max_nm=1500;
cellbodywidthpix=24;
closedistlimit=50;
resolvenm=10;
tipws=9; 

%% import csv file

[file_ch0, pathname_ch0] = uigetfile( ...
    {'*.csv','csv files (*.csv)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Pick a csvfile for ch0', ...
    'MultiSelect', 'off');

motherpath = [pathname_ch0 '\'];
cd (motherpath)

[file_ch1, pathname_ch1] = uigetfile( ...
    {'*.csv','csv files (*.csv)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Pick a file for ch1', ...
    'MultiSelect', 'off');

motherpath = [pathname_ch1 '\'];
cd (motherpath)


%% theshold filter setting

if SetThreshold==1
    prompt = {'Frame','X range','Y range','Z (nm)','X precision (nm)'...
        ,'Y precision (nm)','Photons','Background','PSF Sigma X (pix)','PSF Sigma Y (pix)','Simga ratio (mean std)'};
    nametitle = (['Set Threshold-parameters']);
    dlg_title = nametitle;
    num_lines=1;

    defaultans = {'34, 100000','130,210','180,260','-400, 400','0,14'...
        ,'0, 14','0, 20000','0, 100','0 ,2.5','0, 2.5', '1.05, 0.15, 3'}; 
    ThTable=inputdlg(prompt,dlg_title,num_lines,defaultans);
    if ~isempty (ThTable{11})
        ratio=str2num(ThTable{11});
        ThTable{11}=[num2str(ratio(1)-ratio(2)*ratio(3)),',',num2str(ratio(1)+ratio(2)*ratio(3))];
    end
else
    ThTable=cell(1,11);

end

%% get filename and ONI_file_ID

fileID_ch0 = [pathname_ch0 file_ch0];
fileID_ch1 = [pathname_ch1 file_ch1];
disp(file_ch0)
disp(file_ch1)

%%  single CSV file and seperate data for each channel

data_ch0=[];
data_ch1=[];
data_th_ch0=[];
data_th_ch1=[];
data_th_ch1_shift=[];
[data_ch0, ~]=readcsv_ch_multifile(fileID_ch0,0);
[~, data_ch1]=readcsv_ch_multifile(fileID_ch1,1);
if ~isempty (data_ch0)
    data_ch0(:,17)=data_ch0(:,15)./data_ch0(:,16);
end
if ~isempty (data_ch1)
    data_ch1(:,17)=data_ch1(:,15)./data_ch1(:,16);
end

PixSizeimage_ch0=[];
PixSizeimage_ch1=[];
Stat_th_ch0=[];
Stat_th_ch1=[];
Stat_ch0=[];
Stat_ch1=[];

Scell={char(file_ch0),char(file_ch1)};
common_filename = Scell{1}(1:find(any(diff(char(Scell(:))),1),1,'first')-1);
if isempty(common_filename)
    error ('check the files')
end

commonname=[common_filename];
newfoldername=commonname;

if exist (newfoldername,'dir')==0
    mkdir(motherpath,newfoldername)
    cd(motherpath);
    cd (newfoldername);
    filepath=[pwd '\'];
else
    cd(motherpath);
    cd(newfoldername)
end

%% ROI subfolder

if isempty(ThTable{2,1})
    range_x=[1,427];
else
    xrange=ThTable{2,1};
    range_x=str2num(xrange);
end

if isempty(ThTable{3,1})
    range_y=[1,684];
else
    yrange=ThTable{3,1};
    range_y=str2num(yrange);
end


ROIfoldername=['ROI_',num2str(range_x(1)),'_',num2str(range_x(2)),'_',num2str(range_y(1)),'_',num2str(range_y(2)),'_',common_filename];

if exist (ROIfoldername,'dir')==0
    mkdir(ROIfoldername);
    cd(motherpath);
    cd (newfoldername);
    cd(ROIfoldername);
    filepath=[pwd '\'];
else
    cd(motherpath);
    cd(newfoldername);
    cd(ROIfoldername);
end

%% Data_threshold and STORM image

T0=[];
% file header
head={'Channel','Frame','X_nm','Y_nm','Z_nm',...
    'X_precision_nm','Y_precision_nm','X_pix','Y_pix','Z_pix',...
    'X_precision_pix','Y_precision_pix','Photons','Background',...
    'PSF_Sigma_X_pix','PSF_Sigma_Y_pix','sigmaRatio_sXoversY'};

switch SetThreshold
    case 0
        if ch==0 || ch==2 && ~isempty(data_ch0)
            savefilename=['ch0all_',file_ch0(1:end-4)];
            save(savefilename,'data_ch0')
            if drawFig==1
                savefilename=[savefilename];
                [PixSizeimage_ch0]=ONI_Img_pix(data_ch0,ThTable,dotsize,savefilename,0,C_range,writepixsizeimage);
                save(['STORMpiximg_',savefilename],'PixSizeimage_ch0')
            end
        end
        
        if ch==1 || ch==2 && ~isempty(data_ch1)
            savefilename=['ch1all_',file_ch1(1:end-4)];
            save(savefilename,'data_ch1')
            if drawFig==1
                savefilename=[savefilename];
                [PixSizeimage_ch1]=ONI_Img_pix(data_ch1,ThTable,dotsize,savefilename,1,C_range,writepixsizeimage);
            end
        end
        
    case 1
        if ch==0 || ch==2 &&~isempty(data_ch0)
            savefilename=['ch0Th_',file_ch0(1:end-4)];
            
            [data_th_ch0]=Threshold_ONI(data_ch0,savefilename,ThTable,prompt,1);
            
            save(savefilename,'data_th_ch0')
            if ~isempty(data_th_ch0) && drawSTORMpiximg==1
                T0 = array2table(data_th_ch0,'VariableNames',head);
                writetable(T0,[savefilename,'.csv']);
                if drawFig==1
                    savefilename=[savefilename];
                    [PixSizeimage_ch0]=ONI_Img_pix(data_th_ch0,ThTable,dotsize,savefilename,0,C_range,writepixsizeimage);               
                end
            end
            close all
            
        end
        
        close all
        if ch==1 || ch==2 &&~isempty(data_ch1)
            savefilename=['ch1Th_',file_ch1(1:end-4)];
            [data_th_ch1]=Threshold_ONI(data_ch1,savefilename,ThTable,prompt,1); 
            if ~isempty(data_th_ch1)
                save(savefilename,'data_th_ch1')
                T1 = array2table(data_th_ch1,'VariableNames',head);
                writetable(T1,[savefilename,'.csv']);
                
                if drawFig==1
                    savefilename=[savefilename];
                    [PixSizeimage_ch1]=ONI_Img_pix(data_th_ch1,ThTable,dotsize,savefilename,1,C_range,writepixsizeimage);
                end
            end
            close all
        end
end


%% channel drift correctin (based on pixelized image; ch0 is fixed; ch1 is corrected)
if do_driftcorrect==1
    if ch==2 && merge_ch==1 && ~isempty(PixSizeimage_ch0) && ~isempty(PixSizeimage_ch1) && manual_driftcorrect==0
        switch shiftRef
            case 'STORMPixImage'
                [xshift,yshift,xshift_nm,yshift_nm]=shiftcorrect_2DnormXcorr(PixSizeimage_ch0, PixSizeimage_ch1,applyfilteronimage,...
                    shiftfindoption,ws,pixelsize,ApplyfilterCC,filtertype,kernelsizeCC,filtersigma,showpeakchoice,ccresultsave);
            case 'LocImage'%
                %Load LocImage
                
                TiffFile=[file(1:ind(2)-1),'_locImage.tiff'];
                Tiffid = [pathname  TiffFile];
                LocImage=imread(Tiffid);
                
                if ch==2 % for both channel
                    LocImage_ch0=double(LocImage(:,:,2));
                    LocImage_ch1=double(LocImage(:,:,1));
                    imshowpair(LocImage_ch0,LocImage_ch1);
                elseif ch==0 % for first channel
                    LocImage_ch0=double(LocImage(:,:,2));
                elseif ch==1 % for second channel
                    LocImage_ch1=double(LocImage(:,:,1));
                end
                
                [xshift,yshift,xshift_nm,yshift_nm]=shiftcorrect_2DnormXcorr(LocImage_ch0, LocImage_ch1,applyfilteronimage,...
                    shiftfindoption,ws,pixelsize,ApplyfilterCC,filtertype,kernelsizeCC,filtersigma,showpeakchoice,ccresultsave);
            
        end
        
    elseif ch==2 && merge_ch==1 && ~isempty(PixSizeimage_ch0) && ~isempty(PixSizeimage_ch1) && manual_driftcorrect==1
        manualdriftinput=inputdlg({'driftx(nm)','drifty(nm)'});
        xshift=str2num(manualdriftinput{1})/pixelsize;
        yshift=str2num(manualdriftinput{2})/pixelsize;
    end
    
    data_ch1_shift=data_ch1;
    data_ch1_shift(:,8)=data_ch1_shift(:,8)-xshift;
    data_ch1_shift(:,9)=data_ch1_shift(:,9)-yshift;

    switch SetThreshold
        case 0
            savefilename=[commonname];
            ONI_storm(data_ch0,data_ch1,xshift,yshift,ThTable,dotsize,savefilename);
            savefilename=['ch1_shift_' file_ch1(1:end-4)];
            save(savefilename,'ch1_shift')
            
        case 1
            savefilename=[commonname];
            [data_th_ch1_shift]=Threshold_ONI(data_ch1_shift,savefilename,ThTable,prompt,0); 
            savefilename=['ch1Th_shift_' file_ch1(1:end-4)];
            save(savefilename,'data_th_ch1_shift')
            
            if ~isempty(data_th_ch1_shift)
                T0=[];
                T0 = array2table(data_th_ch1_shift,'VariableNames',head);
                writetable(T0,[savefilename,'.csv']);
                savefilename=[commonname];
                ONI_storm(data_th_ch0,data_th_ch1,xshift,yshift,ThTable,dotsize,savefilename);
            end
    end
    
end
close all
%% Merge the 2 channel STORM images

if merge_ch==1 && ch==2 && SetThreshold==0 && ~isempty(data_ch0) && ~isempty(data_ch1) && do_driftcorrect==0
    savefilename=['STORM_2D_Col_Mer_',commonname];
    stormimage2DColor_merge (data_ch0(:,8),data_ch0(:,9),data_ch1(:,8),data_ch1(:,9),range_x,range_y,dotsize,savefilename)
    
elseif merge_ch==1 && ch==2 && SetThreshold==0 && ~isempty(data_ch0) && ~isempty(data_ch1) && do_driftcorrect==1
    savefilename=['STORM_2D_Col_Mer_shift_',commonname];
    stormimage2DColor_merge(data_ch0(:,8),data_ch0(:,9),data_ch1_shift(:,8),data_ch1_shift(:,9),range_x,range_y,dotsize,savefilename)
    
elseif merge_ch==1 && ch==2 && SetThreshold==1 && ~isempty(data_th_ch0) && ~isempty(data_th_ch1_shift) && do_driftcorrect==1
    savefilename=['STORM_2D_Col_Mer_th_shift',commonname];
    stormimage2DColor_merge(data_th_ch0(:,8),data_th_ch0(:,9),data_th_ch1_shift(:,8),data_th_ch1_shift(:,9),range_x,range_y,dotsize,savefilename)
end


%% rendering images
if Supres==1 && SetThreshold==1
    if isempty(ThTable{2})
        range_x=[1,427];
        range_y=[1,684];
    else
        range_x=str2num(ThTable{2});
        range_y=str2num(ThTable{3});
    end
    if ch==0 || ch==2 && ~isempty(data_th_ch0)
        X=(data_th_ch0(:,8));Y=(data_th_ch0(:,9));
        
        savefilename=[file_ch0(1:end-4),'_',num2str(range_x(1)),'_',num2str(range_x(2)),'_',num2str(range_y(1)),'_',num2str(range_y(2))];
        [Ch0nm]=renderimage_nm(Y,X,RendSize,range_y,range_x,pixelsize,1,savefilename,RendGaussfilterSigma);
    end
    clear X Y
    if ch==1 || ch==2 && ~isempty(data_th_ch1) && do_driftcorrect==0
        savefilename=[file_ch1(1:end-4),'_',num2str(range_x(1)),'_',num2str(range_y(2)),'_',num2str(range_y(1)),'_',num2str(range_y(2))];
        if merge_ch==0
            X=(data_th_ch1(:,8));Y=(data_th_ch1(:,9));
            
            [Ch1nm]=renderimage_nm(Y,X,RendSize,range_y,range_x,pixelsize,0,savefilename,RendGaussfilterSigma);
            
        elseif merge_ch==1
            X=(data_th_ch1(:,8));Y=(data_th_ch1(:,9));
            Z=(data_th_ch1(:,5));
            [Ch1nm]=renderimage_nm(Y,X,RendSize,range_y,range_x,pixelsize,0,savefilename,RendGaussfilterSigma);
          
        end
        
    elseif ch==1 || ch==2 && ~isempty(data_th_ch1) && do_driftcorrect==1
        savefilename=[file_ch1(1:end-4),'_',num2str(range_x(1)),'_',num2str(range_y(2)),'_',num2str(range_y(1)),'_',num2str(range_y(2))];
        if merge_ch==0
            X=(data_th_ch1_shift(:,8));Y=(data_th_ch1_shift(:,9));
            
            [Ch1nm]=renderimage_nm(Y,X,RendSize,range_y,range_x,pixelsize,0,savefilename,RendGaussfilterSigma);
        
        elseif merge_ch==1
            X=(data_th_ch1_shift(:,8));Y=(data_th_ch1_shift(:,9));
            Z=(data_th_ch1_shift(:,5));
            [Ch1nm]=renderimage_nm(Y,X,RendSize,range_y,range_x,pixelsize,0,savefilename,RendGaussfilterSigma);
         
        end
    end
    clear X Y
    

end
close all
%% statistic analyisis
if do_stat==1
    if ch==0 || ch==2
        
        if SetThreshold==0 && ~isempty(data_ch0)
            [Stat_ch0]=statistic_fig_2D_3D(data_ch0,0.1,['stat_',file_ch0(1:end-4)]);
        elseif SetThreshold==1 && ~isempty(data_th_ch0)
            [Stat_th_ch0]=statistic_fig_2D_3D(data_th_ch0,0.1,['stat_th_',file_ch0(1:end-4)]);
        end
    end
    if ch==1 || ch==2
        if SetThreshold==0 && ~isempty(data_ch1) &&  do_driftcorrect==0
            [Stat_ch1]=statistic_fig_2D_3D(data_ch1,0.1,['stat_',file_ch1(1:end-4)]);
        elseif SetThreshold==0 && ~isempty(data_ch1_shift) &&  do_driftcorrect==1
            [Stat_ch1]=statistic_fig_2D_3D(data_ch1_shift,0.1,['stat_',file_ch1(1:end-4)]);
        elseif SetThreshold==1 && ~isempty(data_th_ch1) && do_driftcorrect==0
            [Stat_th_ch1]=statistic_fig_2D_3D(data_th_ch1,0.1,['stat_th_',file_ch1(1:end-4)]);
        elseif SetThreshold==1 && ~isempty(data_th_ch1_shift) && do_driftcorrect==1
            [Stat_th_ch1]=statistic_fig_2D_3D(data_th_ch1_shift,0.1,['stat_th_',file_ch1(1:end-4)]);
        end
    end
end

%pause
close all

%% Tip position

    ch1=Ch0nm;
    ch2=Ch1nm;
    BW1=ch1>0;
    BW2=ch2>0; 
    BW3=BW1+BW2; BW3(BW3>0)=1;

   [~,InBG,InEdgeBW] = segmentation(imadd(ch1,ch2),BW3);
    level=graythresh(uint8(ch1));
    Thresh_ch1=level*255;
    BW1= imbinarize(ch1,Thresh_ch1);
    
    level=graythresh(uint8(ch2));
    Thresh_ch2=level*255;
    BW2= imbinarize(ch2,Thresh_ch2);
    BW3=BW1+BW2; BW3(BW3>0)=1;

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
    SegL=InBG*2+CB*1+MV*3; % assign unique value for CB(1), MV(3) OutBg (0) and InBG (2)
    OutBG=SegL==0;

    [SegL_new] = remove_islands(SegBW,SegL);

    MVlengthlimit_min = (MVlengthlimit_min_nm/resolvenm)+cellbodywidthpix;
    MVlengthlimit_max = (MVlengthlimit_max_nm/resolvenm)+cellbodywidthpix;
    [tiploc]=get_tiplocation(InEdgeBW,SegL_new,MVlengthlimit_min,MVlengthlimit_max,closedistlimit);
  
    figure
    axis image
    Xrange=range_x;
    xm=Xrange(1)-.5; xp=Xrange(2)+.5;
    Yrange=range_y;
    ym=Yrange(1)-.5; yp=Yrange(2)+.5;

    X=data_th_ch0(:,8);
    Y=data_th_ch0(:,9);
    
    X2=data_th_ch1_shift(:,8);
    Y2=data_th_ch1_shift(:,9);
    
    scatter(X2,Y2,1,'g','filled'); axis image
    hold on
    scatter(X,Y,1,'m','filled'); axis image
    hold on 
    
    tipx=tiploc(2,:);
    tipy=tiploc(1,:);

    tipX_ori=range_x(1)+floor(tipx*resolvenm/117)+0.5;
    tipY_ori=range_y(2)-floor(tipy*resolvenm/117)+0.5;
    tiploc_xyori(1,:)=tipX_ori;
    tiploc_xyori(2,:)=tipY_ori;


    
    scatter(tipX_ori,tipY_ori,30,'y+'); axis image
    
    xlim([xm xp]);
    ylim([ym yp]);
    set(gca,'box','on')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'XLabel',[]);
    set(gca,'YLabel',[]);
    set(gca,'XTickLabel',[]);
    set(gca,'XTickLabel',[]);
    set(gca,'title',[]);
    set(gca,'Color','k')
    hold off
        
%% pair correlation analysis

    data1.x=data_th_ch0(:,8);
    data1.y=data_th_ch0(:,9);
    data1.z=data_th_ch0(:,10);
    data2.x=data_th_ch1_shift(:,8);
    data2.y=data_th_ch1_shift(:,9);
    data2.z=data_th_ch1_shift(:,10);
    savefilename='paircorr-3D';
    [paircorr_3D, set_paircorr_3D]=tipbytipanalysis_3D_Paircorr(savefilename,tiploc_xyori,data1,data2,tipws,resolvenm,pixelsize,1);

%%


function[data_ch0, data_ch1]=readcsv_ch_multifile(filename,ch)
head=[];
data_ch0=[];
data_ch1=[];

[csvdata]=csvread(filename,1,0);
numMol_ch1= nnz(csvdata(:,1));
numMol_ch0=size(csvdata,1)-numMol_ch1;

if ch== 0 || ch== 2
    if numMol_ch0==0
        
        data_ch0=[];
    else
        [ch0row]=find(~(csvdata(:,1)));
        data_ch0=csvdata(ch0row,:);
        
        if nnz(data_ch0(:,1))>0
            error ('Check the Channel 0 values')
        end
    end
end

if ch==1 || ch==2
    if numMol_ch1==0
        data_ch1=[];
    else
        [ch1row]=find(csvdata(:,1));
        data_ch1=csvdata(ch1row,:);
        if nnz(~data_ch1(:,1))>0
            error ('Check the Channel 1 values')
        end
    end
end
end

function[data_ch0, data_ch1]=readcsv_ch(filename,ch)
head=[];
data_ch0=[];
data_ch1=[];

[csvdata]=csvread(filename,1);
numMol_ch1= nnz(csvdata(:,1));
numMol_ch0=size(csvdata,1)-numMol_ch1;

if ch== 0 || ch== 2
    if numMol_ch0==0
        
        data_ch0=[];
    else
        data_ch0=csvdata(1:numMol_ch0,:);
        
        if nnz(data_ch0(:,1))>0
            error ('Check the Channel 0 values')
        end
    end
end

if ch==1 || ch==2
    if numMol_ch1==0
        data_ch1=[];
    else
        data_ch1=csvdata(numMol_ch0+1:end,:);
        if nnz(~data_ch1(:,1))>0
            error ('Check the Channel 1 values')
        end
    end
end
end

function[PixSizeimage]=ONI_Img_pix(data,ThTable,dotsize,savefilename,disp_ch,C_range,writepixsizeimage)
PixSizeimage=[];
colormap(jet)
xrange=ThTable{2};
yrange=ThTable{3};
if size (xrange,2)==0
    xm=0.5; xp=427.5;
else
    Xrange=str2num(xrange);
    xm=Xrange(1)-.5; xp=Xrange(2)+.5;
end

if size(yrange,2)==0
    ym=0.5; yp=684.5;
else
    Yrange=str2num(yrange);%
    ym=Yrange(1)-.5; yp=Yrange(2)+.5;
end
ImageSize=[yp-ym, xp-xm];

if sum(data(:,5))==0
    figure
    switch disp_ch
        case 0
            figure(1)
            scatter(data(:,8),data(:,9),dotsize,'filled','bo')
            title(savefilename,'FontSize',6, 'Interpreter','none')
            xlabel('X (pix)')
            ylabel('Y (pix)')
        case 1
            figure(2)
            scatter(data(:,8),data(:,9),dotsize,'filled','ro')
            title(savefilename,'FontSize',6, 'Interpreter','none')
            xlabel('X (pix)')
            ylabel('Y (pix)')
            
    end
    
    daspect([1 1 1])
    xlim([xm xp]);
    ylim([ym yp]);
    set(gca,'box','on')
    
else
    if isempty (C_range)
        C_range=[min(data(:,5)), max(data(:,5))];
    end
    
    switch disp_ch
        case 0
            figure
            scatter(data(:,8),data(:,9),dotsize,'filled','CData',data(:,5));
            title(['Z_',savefilename],'FontSize',6, 'Interpreter','none')
            caxis(C_range)
            c= colorbar;
            c.Label.String='Z(nm)';
            set(gca,'box','on')
            daspect([1 1 1])
            xlim([xm xp]);
            ylim([ym yp]);
            xlabel('X (pix)')
            ylabel('Y (pix)')
          
            
        case 1 % for ch1
            
            figure
            scatter(data(:,8),data(:,9),dotsize,'filled','CData',data(:,5));
            title(['Z_',savefilename],'FontSize',6, 'Interpreter','none')
            caxis(C_range)
            c= colorbar;
            c.Label.String='Z(nm)';
            daspect([1 1 1])
            xlim([xm xp]);
            ylim([ym yp]);
            xlabel('X (pix)')
            ylabel('Y (pix)')
            set(gca,'box','on')
    
    end
end

%---- write image
    set(gca,'color','black')
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'box','on'); 
    set(gca,'YColor','k')
    set(gca,'XColor','k')
    pause (1)
    F = getframe;
    Zimage=F.cdata;
    PixSizeimage=imresize(Zimage, ImageSize);
    
    if writepixsizeimage ==1
    save(['STORMpiximg_',savefilename],'PixSizeimage')
    imwrite(PixSizeimage,['STORMpiximg_',savefilename,'.tif'],'Compression','none');
    end
close all
end

function[data_th]=Threshold_ONI(data,savefilename,ThTable,prompt,display)
col= [2,8,9,5:7,13:16,17];
data_th=data;
for pr=1:11
    if isempty(data_th)
        nummol(pr)=0;
        fprintf('%d\n',nummol)
        
    else
        th_range = str2num(ThTable{pr});
        if isempty (th_range)
            nummol(pr)=size(data_th,1);
        else
            
            [data_th]=DataThreshold(data_th,col(pr),th_range);
            if isempty (data_th)
                nummol(pr)= 0;
            else
                nummol(pr)=size(data_th,1);
            end
        end
    end
end

%---- draw a figure
if display==1
    if nummol~=0
        figure (1)
        plot(nummol,'Marker','diamond'),
        ylabel({'number of molecuels'});
        xlabel({'Parameters'});
        xt=(1:11);
        addtick=xt+0.4;
        cnt=0;
        for k=1:10
            cnt=k*2;
            newXtickLabel{cnt-1}=prompt{k};
            newXtickLabel{cnt}=ThTable{k};
            newXtick(cnt-1)=xt(k);
            newXtick(cnt)=addtick(k);
        end
        set(gca,'xtick',newXtick,'XTickLabel',newXtickLabel,'XTickLabelRotation',45)
        set(gca,'TickLength',[0 0]);
        title(savefilename,'FontSize',6, 'Interpreter','none')
        
    end
end
end

function[data_th]=DataThreshold(data,col,th_range)
data_th=[];
[r] = find((data(:,col) >= th_range(1) & data(:,col) <= th_range(2)));
if isempty(r)
    data_th=[];
else
    data_th=data(r,:);
end
end

function [xshift,yshift,xshift_nm,yshift_nm]=shiftcorrect_2DnormXcorr(RefImage,ShiftImage,applyfilteronimage,shiftfindoption,ws,Pixelsize,ApplyfilterCC,filtertype,kernelsizeCC,filtersigma,showpeakchoice,ccresultsave)

if length(size(RefImage))>2
    RefImage=rgb2gray(RefImage);
    RefImage=double(flipud(RefImage));
end

if length(size(ShiftImage))>2
    ShiftImage=rgb2gray(ShiftImage);
    ShiftImage=double(flipud(ShiftImage));
end
pn=(ws-1)/2;
x_range=[1,size(RefImage,2)];
y_range=[1,size(RefImage,1)];

%---- remove edges

RefImage(:,1)=0; % vertical left
RefImage(:,end)=0; % vertical right
RefImage(1,:)=0; % horizontal up
RefImage(end,:)=0; % horizontal down

ShiftImage(:,1)=0; % vertical left
ShiftImage(:,end)=0; % vertical right
ShiftImage(1,:)=0; % horizontal up
ShiftImage(end,:)=0; % horizontal down


%---- apply filter ; gaussian blur
if applyfilteronimage==1
    ShiftImage= imgaussfilt(ShiftImage, 1);
    RefImage = imgaussfilt(RefImage, 1);
end

%--- norxcorr2 - 2D cross correlation

cmean = normxcorr2(RefImage,ShiftImage);

if ApplyfilterCC==1
    cmean=FilterFrameOut_ws(cmean,filtertype,kernelsizeCC,filtersigma);
end

ori_xpeak = size(cmean,2)-x_range(2)+1;
ori_ypeak = size(cmean,1)-y_range(2)+1;

%---- find xpeak,ypeak
[max_cc]= max(cmean(:));
[ypeak, xpeak] =find(cmean==max_cc);

Q=cmean(ypeak-pn:ypeak+pn,xpeak-pn:xpeak+pn);

[max_cmean]= max(Q(:));
[ypeakQ, xpeakQ] =find(Q==max_cmean); %

if showpeakchoice==0
    choice='yes';
    
elseif showpeakchoice==1
    close all
    figure
    imagesc(cmean)
    hold on
    plot(xpeak,ypeak,'+')
    axis image
    title(['xpeak:',num2str(xpeak),' ypeak',num2str(ypeak)])
     choice = questdlg('use this peak?', ...
        'change peak', ...
        'yes','no','manual','yes');
end

switch choice
    case 'yes'
        fprintf('xpeak : %d  ypeak: %d \n',xpeak,ypeak)
    case 'second peak'
        
        [ii] = sort(cmean(:));
        secondmax = cmean(ii([2,end-1]));
        if size(secondmax,2)==1
            
            hold on
            plot(xpeak,ypeak,'bx')
            title('second peak')
            fprintf('xpeak : %d  ypeak: %d \n',xpeak,ypeak)
            
        else
            error('more than two values for secondpeak ')
        end
        fprintf('xpeak : %d  ypeak: %d \n',xpeak,ypeak)
    case 'manual'
        
        prompt = {'Enter new xpeak:','Enter new ypeak:'};
        dlg_title = 'Input';
        num_lines = 1;
        
        defaultans = {num2str(ori_xpeak),num2str(ori_ypeak)};
        newpeak = inputdlg(prompt,dlg_title,num_lines,defaultans);
        xpeak=str2double(newpeak{1});
        ypeak=str2double(newpeak{2});
        
        hold on
        plot(xpeak,ypeak,'gx')
        title('manual peak')
        fprintf('xpeak : %d  ypeak: %d \n',xpeak,ypeak)
end

%---- choose a fitting method
switch shiftfindoption
    case 1 % center of mass
        props = regionprops(true(size(Q)), Q, 'WeightedCentroid');
        WeightedCentroid=props.WeightedCentroid;
        
        x_fit=WeightedCentroid(1);
        y_fit=WeightedCentroid(2);
        
        xshift=-(ori_xpeak-xpeak-(x_fit)+((ws+1)/2));
        yshift=-(ori_ypeak-ypeak-(y_fit)+((ws+1)/2));
        
    case 2 % 2D gaussian fitting
        [bb,ci] = fit2dgaussian(Q,[(ws+1)/2;(ws+1)/2],ws);
        x_fit=bb(2);
        y_fit=bb(3);
        
        xshift=-(ori_xpeak-xpeak-x_fit+((ws+1)/2));
        yshift=-(ori_ypeak-ypeak-y_fit+((ws+1)/2));
end


%---- draw figure

if ccresultsave==1
    close all
    figure('Name','xydrift');
    set(gcf, 'Position', [50    50    1200    500]);
    
    subplot(1,3,1)
    imshowpair(RefImage,ShiftImage);
    set(gca,'YDir','normal');
    title('compare images')
    axis image
    xlim([x_range(1) x_range(2)]);
    ylim([y_range(1) y_range(2)]);
    
    subplot(1,3,2)
    imagesc(cmean)
    hold on
    plot(xpeak,ypeak,'+')
    axis image
    title(['xpeak:',num2str(xpeak),' ypeak',num2str(ypeak)])
    
    subplot(1,3,3)
    imagesc(Q);
    set(gca,'YDir','normal');
    axis image
    title(['max peak'])
    hold on
    plot(xpeakQ,ypeakQ,'b+')
    hold on
    plot(x_fit,y_fit,'r+')
    
    xshift_nm=xshift*Pixelsize;
    yshift_nm=yshift*Pixelsize;
    
    switch shiftfindoption
        case 1
            title(['Weighted Centroid ']);
            
            xlabel({[ 'xshift (pix): ', num2str(xshift,'%0.3f'),' ; yshift (pix):', num2str(yshift, '%0.3f')],...
                [ 'xshift (nm): ', num2str(xshift_nm,'%0.1f'),' ; yshift (nm):', num2str(yshift_nm, '%0.1f')] })
        case 2
            title([' 2D Gaussian Fit'])
            
            xlabel({[ 'xshift (pix): ', num2str(xshift,'%0.3f'),' ; yshift (pix):', num2str(yshift, '%0.3f')],...
                [ 'xshift (nm): ', num2str(xshift_nm,'%0.1f'),' ; yshift (nm):', num2str(yshift_nm, '%0.1f')] })
            
    end
    savefig(gcf,'xcoor2_shift')
    close
else
    xshift_nm=xshift*Pixelsize;
    yshift_nm=yshift*Pixelsize;
end

fprintf('xshift (pix): %.3d  yshift(pix): %.3d \n',xshift,yshift)
fprintf('xshift (nm): %d  yshift(nm): %d \n',round(xshift_nm),round(yshift_nm))

end

function [b,ci] = fit2dgaussian(data,loc,ws)

% Inputs:
%   loc -   a list of central (x,y) pixels.
%   ws -    the size of the region of interest
%   data -  the full matrix in which the 2d gaussian is embeded
%   Note -  this fits the data into a symmetric 2d gaussian, which might
%           not be the case in the cross correlation of the images.

% Outputs:
%   b -     5 parameters: [Intensity x y sigma bkd]
%   ci -    upper and lower limit of the estimate

n=1; % assuming you are only interested in the central peak
[X, Y] = meshgrid(1:ws, 1:ws); xy_data = [X(:) Y(:)];
A = double(data(loc(1,n)-(ws-1)/2:loc(1,n)+(ws-1)/2,loc(2,n)-(ws-1)/2:loc(2,n)+(ws-1)/2));
estimates = fit2dguasswbg(A);
[b,r1,J] = nlinfit(xy_data, reshape(A,[ws^2,1]), @mygauss3, estimates);
ci = nlparci(b,r1,J);
b(2)=loc(1,n)+b(2)-ceil(ws/2);  b(3)=loc(2,n)+b(3)-ceil(ws/2);

end

function [estimates,fval,model,FittedCurve] = fit2dguasswbg(mydata, firstguess)

switch nargin
    case 1
        start_point = [max(mydata(:)),ceil(size(mydata,1)/2),ceil(size(mydata,1)/2),0.703,mean(mydata(:))];
    case 2
        start_point=firstguess;
end

start_point=double(start_point);
[X,Y]=meshgrid(1:size(mydata,1),1:size(mydata,1));
xy_data=[X(:),Y(:)];

model = @guass2dfun;
[estimates,fval] = fminsearch(model, start_point);
[v,FittedCurve] = guass2dfun(estimates);
FittedCurve = reshape(FittedCurve,[size(mydata,1) size(mydata,2)]);

% guass2dfun accepts 2d gaussian parameters as inputs, and outputs sse,
% the sum of squares error for the measured and the FittedCurve.
% FMINSEARCH only needs sse, but we want to plot the FittedCurve at the end.
    function [sse, FittedCurve] = guass2dfun(params)
        A = params(1);
        x0 = params(2);
        y0 = params(3);
        sig = params(4);
        bg = params(5);
        FittedCurve = A .* exp(-((xy_data(:,1)-x0).^2 + (xy_data(:,2)-y0).^2)./(2*sig^2)) + bg;
        ErrorVector = FittedCurve - reshape(mydata,[size(mydata,1)^2,1]);
        sse = sum(ErrorVector .^ 2);
    end
end

function [H] = mygauss3(p,xy_dat)
% Creates a 2D gaussian.
% Inputs -
% X,Y  - x and y values created using the meshgrid function
% p - a list of parameter : amplitude, x-center, y-center, width,
% background

H = p(1)*exp( - ((xy_dat(:,1)-p(2)).^2 + (xy_dat(:,2)-p(3)).^2)./(2*p(4)^2))+p(5) ;
end

function [C] = FilterFrameOut_ws(A,filtertype,filtersize,filtersigma)
B=zeros(size(A));
C=zeros(size(A));
switch filtertype
    
    % edges = imfilter(AA,h,'replicate'); AA=images
    % h = fspecial('log', 14, 0.5); filtertype: log, 14= kernel size (14),
    % gaussian sigma
    %h = fspecial('log', hsize, sigma)  returns a rotationally symmetric Laplacian of Gaussian filter of size hsize with standard deviation sigma (positive). hsize can be a vector specifying the number of rows and columns in h, or it can be a scalar, in which case h is a square matrix. The default value for hsize is [5 5] and 0.5 for sigma.
    % median filter
    
    case 'MedianFilter'
        B = medfilt2(A, [filtersize filtersize]);
        C=A-B;
    case 'MeanFilter'
        h = fspecial('average', [filtersize filtersize]);
        B = imfilter(A,h,'replicate');
        
        C=A-B;
    case 'DiscFilter'
        h = fspecial('disc', filtersize/2);
        B = imfilter(A,h,'replicate');
        C=A-B;
        
    case'LogFilter'
        h = fspecial('log', 14, 0.5);
        B = imfilter(A,h,'replicate');
        C=A-B;
        
    case'GaussFilter'
        h = fspecial('gaussian', [filtersize filtersize],filtersigma) ;
        B = imfilter(A,h,'replicate');
        C=A-B;
        
end
end

function[W,posxyz,posxyz_nm]=renderimage_nm(X,Y,RendSize,range_x,range_y,pixelsize,rescale,savefilename,gaussfiltersigma)

    ROIPix_X=range_y(2)-range_y(1)+1;
    ROIPix_Y=range_x(2)-range_x(1)+1;
    X0nm=((X-range_x(1))*pixelsize); 
    Y0nm=((Y-range_y(1))*pixelsize); 
    

IndYnm=(Y0nm/RendSize);
IndXnm=(X0nm/RendSize);
IndY=floor(Y0nm/RendSize);
IndX=floor(X0nm/RendSize);

ROInmX=floor((ROIPix_X*pixelsize)/RendSize) ; 
ROInmY=floor((ROIPix_Y*pixelsize)/RendSize) ;

IndX(IndX==0)=1;
IndY(IndY==0)=1;
IndX(IndX>ROInmY)=ROInmY; 
IndY(IndY>ROInmX)=ROInmX; 

W=zeros(ROInmX,ROInmY); 
LinearInd=sub2ind(size(W),IndY,IndX);

for i=1:size(LinearInd,1)
    if W(LinearInd(i))>=1
        W(LinearInd(i))= W(LinearInd(i))+1;
    else
        W(LinearInd(i))=1;
    end
end

IndZ=zeros(size(LinearInd));
for k=1:size(LinearInd,1)
    IndZ(k)=W(IndY(k),IndX(k));
end
posxyz=[IndX,IndY,IndZ];
posxyz_nm=[IndXnm,IndYnm,IndZ];

if rescale==1
    W = (W./max(W(:)))*255;
end

if gaussfiltersigma~=0
    
    W = imgaussfilt(W,gaussfiltersigma);
    imshow(W)
    set(gca,'YDir','normal')
    
end
W=W';

W=flipud(W);

imwrite(W,['Rend',num2str(RendSize),'_',savefilename,'.tif'],'Compression','none');
end

function[]=ONI_storm(data1,data2,xshift,yshift,thTable,dotsize,savefilename)
if isempty(thTable{2,1})

    xm=0.5; xp=427.5;
    
else
    xrange=thTable{2,1};
    Xrange=str2num(xrange);
    xm=Xrange(1)-.5; xp=Xrange(2)+.5;
end

if isempty(thTable{3,1})
    ym=0.5; yp=684.5;
else
    yrange=thTable{3,1};
    Yrange=str2num(yrange);
    ym=Yrange(1)-.5; yp=Yrange(2)+.5;
end

ImageSize=[yp-ym, xp-xm];

if sum(data1(:,5))==0
    close all
    figure
    scatter(data1(:,8),data1(:,9),dotsize,'filled','bo')  
    title(['ch0'],'FontSize',10, 'Interpreter','none')
    xlabel('X (pix)')
    ylabel('Y (pix)')
    daspect([1 1 1])
    xlim([xm xp]);
    ylim([ym yp]);  
    set(gca,'box','on')
    savefig(gcf,['Storm_ch0',savefilename]);
    
    figure
     scatter(data2(:,8)-xshift, data2(:,9)-yshift, dotsize,'filled','ro')
     title(['ch1'],'FontSize',10, 'Interpreter','none')
     xlabel('X (pix)')
     ylabel('Y (pix)')
     daspect([1 1 1])
     xlim([xm xp]);
     ylim([ym yp]);  
     set(gca,'box','on')
     savefig(gcf,['Storm_ch1',savefilename]);

else
    close all
    figure
    scatter(data1(:,8),data1(:,9),dotsize,'filled','o','CData',data1(:,5));
    
    title(['ch0'],'FontSize',10, 'Interpreter','none')
    xlabel('X (pix)')
    ylabel('Y (pix)')
    daspect([1 1 1])
    xlim([xm xp]);
    ylim([ym yp]);  
    set(gca,'box','on')
    savefig(gcf,['Storm_ch0',savefilename]);
    
    figure
    scatter(data2(:,8)-xshift,data2(:,9)-yshift,dotsize,'filled','o','CData',data2(:,5));
    title(['ch1'],'FontSize',10, 'Interpreter','none')
    xlabel('X (pix)')
    ylabel('Y (pix)')
    daspect([1 1 1])
    xlim([xm xp]);
    ylim([ym yp]);  
    set(gca,'box','on')
    savefig(gcf,['Storm_ch1',savefilename]);   
end

end

function[T]=statistic_fig_2D_3D(data,dotsize,savefilename)
close all

if sum(data(:,5))==0
    datatype='2D';
else
    datatype='3D';
end

figure
%% channel
T=[];
%% ch
if sum(data(:,1))==0
    ch=0;
elseif mean(data(:,1))==1
    ch=1;
else
    error ('only for single channel')
end

%% 2D or 3D

switch datatype
    case ('2D')
        set(gcf, 'Position', [50    50    1200    700]);
        
        %---- draw image
        subplot(2,4,1)
        if ch==0
            scatter(data(:,8),data(:,9),dotsize,'filled','o','r');
        elseif  ch==1
            scatter(data(:,8),data(:,9),dotsize,'filled','o','r');
        else
            error 'use single channel data'
            
        end
        
        axis square
        daspect([1 1 1])
        if nargin==3
            title(['ch ',num2str(ch),'_',savefilename],'FontSize',6, 'Interpreter','none')
        else
            title(['ch ',num2str(ch)],'FontSize',6, 'Interpreter','none')
        end
        set(gca,'box','on')
        daspect([1 1 1])
        xlabel('X (pix)')
        ylabel('Y (pix)')
        
        %---- plot frame
        subplot(2,4,2)
        mol=(1:size(data,1));
        plot(data(:,2),mol)
        title (['total N: ', num2str(size(data,1))])
        ylabel('# of molecules','FontSize',14');
        xlabel('frame','FontSize',14');
        
        %---- X-Precision and Y Precision
        subplot(2,4,3);
        histogram(data(:,6),'facecolor','m','facealpha',.2,'edgecolor','none');
        pdx = fitdist(data(:,6),'Normal');
        MuX=pdx.mu;%*PixelSize;
        SigmaX=pdx.sigma;%*PixelSize;
        titleline1=['X: \mu=',num2str(MuX,'%.2f'),' ;\sigma=',num2str(SigmaX,'%.2f'),' (nm)' ];
        ylabel('frequency','FontSize',14');
        hold on
        histogram(data(:,7),'facecolor','c','facealpha',.2,'edgecolor','none');
        pdy = fitdist(data(:,7),'Normal');
        MuY=pdy.mu;%*PixelSize;
        SigmaY=pdy.sigma;%*PixelSize;
        FWHMYnm=SigmaY*2.355;
        titleline2=['Y: mean= ',num2str(MuY,'%.2f'),' ;std= ',num2str(SigmaY,'%.2f'),' (nm)'];
        xlabel('Precision (nm)','FontSize',14');
        ylabel('frequency','FontSize',14');
        title({titleline1;titleline2})
        legend('X','Y')
        legend('boxoff')
        
        %---- precision vs photons
        subplot(2,4,4);
        if size (data,2)==16
            data(:,17)=data(:,15)./data(:,16);
        end
        scatter(data(:,8),data(:,9),1,'cdata',data(:,17))
        axis image
        xlabel('X (pix)','FontSize',14');
        ylabel('Y (pix)','FontSize',14');
        colormap(jet)
        c=colorbar ;
        c.Label.String = 'simga X/ simgaY';
        
        %---- photons
        
        subplot(2,4,5);
        histogram(data(:,13),'edgecolor','none');
        AvePho=mean(data(:,13));
        medPho=median(data(:,13));
        minPho=min(data(:,13));
        maxPho=max(data(:,13));
        title(['mean:',num2str(AvePho,'%.0f'),' min: ',num2str(minPho,'%.0f'),' med:',num2str(medPho,'%.0f'), ' max:' num2str(maxPho)])
        xlabel('photons','FontSize',14');
        ylabel('frequency','FontSize',14');
        
        %---- background
        
        subplot(2,4,6);
        histogram(data(:,14),'edgecolor','none');
        meanbg=mean(data(:,14));
        stdbg=std(data(:,14));
        title(['   mean= ',num2str(meanbg,'%.0f'),' std=', num2str(stdbg,'%.0f')])
        xlabel('background','FontSize',14');
        ylabel('frequency','FontSize',14');
        
        %---- sigma S, sigma Y from fitting (nm scale)
        subplot(2,4,7);
        SX=data(:,15);
        SY=data(:,16);
        histogram(SX,'facecolor','m','facealpha',.2,'edgecolor','none');
        hold on
        histogram(SY,'facecolor','c','facealpha',.2,'edgecolor','none');
        MeanSX=mean(SX);
        STDSX=std(SX);
        FWSX=std(SX)*2.355;
        MeanSY=mean(SY);
        STDSY=std(SY);
        FWSY=std(SY)*2.355;
        firstline=['\sigmaX : mean=',num2str(MeanSX,'%.2f'),'; std= ',num2str(STDSX,'%.2f')];
        secondline= ['\sigmaY : mean=',num2str(MeanSY,'%.2f'),'; std= ',num2str(STDSY,'%.2f')];
        title({firstline;secondline})
        xlabel('\sigmaX & \sigmaY (nm)','FontSize',14');
        ylabel('frequency','FontSize',14');
        legend('\sigmaX','\sigmaY')%
        legend('boxoff')
        
        
        %---- sigmax/sigmay
        subplot(2,4,8);
        sigmaratio=SX./SY;
        histogram(sigmaratio,'edgecolor','none');
        xlabel('\sigmaX/\sigmaY','FontSize',14');
        ylabel('frequency','FontSize',14');
        pdr = fitdist(sigmaratio,'Normal');
        MeanRatio=mean(sigmaratio);
        STDRatio=std(sigmaratio);
        title(['mean=',num2str(MeanRatio,'%.2f'),'; std= ',num2str(STDRatio,'%.2f')])%,' ;FWHM=',num2str(FWSX,'%.1f'),' (nm)'];
        
    case '3D'
        set(gcf, 'Position', [50    50    1200    700]);
        
        %----  localization map
        subplot(2,4,1)
        
        if sum(data(:,5))~=0
            scatter(data(:,8),data(:,9),10,'filled','o','CData',data(:,5));
            c= colorbar;
            c.Label.String='Z(nm)';
        end
        axis square
        daspect([1 1 1])
        if nargin ==3
            title(['ch ',num2str(ch),'_',savefilename],'FontSize',6, 'Interpreter','none')
        else
            title(['ch ',num2str(ch)],'FontSize',6, 'Interpreter','none')
        end
        set(gca,'box','on')
        daspect([1 1 1])
        xlabel('X (pix)')
        ylabel('Y (pix)')
        
        %---- plot frame
        subplot(2,4,2)
        mol=(1:size(data,1));
        plot(data(:,2),mol)
        title (['total N: ', num2str(size(data,1))])
        ylabel('# of molecules','FontSize',14');
        xlabel('frame','FontSize',14');
        
        %---- X-Precision and Y Precision
        subplot(2,4,3);
        histogram(data(:,6),'facecolor','m','facealpha',.2,'edgecolor','none');
        pdx = fitdist(data(:,6),'Normal');
        MuX=pdx.mu;%*PixelSize;
        SigmaX=pdx.sigma;%*PixelSize;
        titleline1=['X: \mu=',num2str(MuX,'%.2f'),' ;\sigma=',num2str(SigmaX,'%.2f'),' (nm)' ];
        ylabel('frequency','FontSize',14');
        hold on
        histogram(data(:,7),'facecolor','c','facealpha',.2,'edgecolor','none');
        pdy = fitdist(data(:,7),'Normal');
        MuY=pdy.mu;%*PixelSize;
        SigmaY=pdy.sigma;%*PixelSize;
        FWHMYnm=SigmaY*2.355;
        titleline2=['Y: mean= ',num2str(MuY,'%.2f'),' ;std= ',num2str(SigmaY,'%.2f'),' (nm)'];
        xlabel('Precision (nm)','FontSize',14');
        ylabel('frequency','FontSize',14');
        title({titleline1;titleline2})
        legend('X','Y')
        legend('boxoff')
        
        %---- Z
        if sum(data(:,5))~=0
            subplot(2,4,4);
            histogram(data(:,5),'edgecolor','none')
            pdz = fitdist(data(:,5),'Normal');
            MuZ=pdz.mu;
            SigmaZ=pdz.sigma;
            FWHMnm=SigmaZ *2.355;
            title(['  mean=', num2str(MuZ,'%.1f'), ' std= ',num2str(SigmaZ,'%.1f')])
            xlabel('z (nm)','FontSize',14');
            ylabel('frequency','FontSize',14');
        end
        
        %---- num of photons
        
        subplot(2,4,5);
        histogram(data(:,13),'edgecolor','none');
        AvePho=mean(data(:,13));
        medPho=median(data(:,13));
        minPho=min(data(:,13));
        maxPho=max(data(:,13));
        title(['mean:',num2str(AvePho,'%.0f'),' min: ',num2str(minPho,'%.0f'),' med:',num2str(medPho,'%.0f'), ' max:' num2str(maxPho)])
        xlabel('photons','FontSize',14');
        ylabel('frequency','FontSize',14');
        
        %---- background
        
        subplot(2,4,6);
        histogram(data(:,14),'edgecolor','none');
        meanbg=mean(data(:,14));
        stdbg=std(data(:,14));
        title(['   mean= ',num2str(meanbg,'%.0f'),' std=', num2str(stdbg,'%.0f')])
        xlabel('background','FontSize',14');
        ylabel('frequency','FontSize',14');
        
        
        %---- sigma S, sigma Y from fitting (nm scale)
        subplot(2,4,7);
        
        SX=data(:,15);
        SY=data(:,16);
        histogram(SX,'facecolor','m','facealpha',.2,'edgecolor','none');
        hold on
        histogram(SY,'facecolor','c','facealpha',.2,'edgecolor','none');
        MeanSX=mean(SX);
        STDSX=std(SX);
        FWSX=std(SX)*2.355;
        MeanSY=mean(SY);
        STDSY=std(SY);
        FWSY=std(SY)*2.355;
        firstline=['\sigmaX : mean=',num2str(MeanSX,'%.2f'),'; std= ',num2str(STDSX,'%.2f')];
        secondline= ['\sigmaY : mean=',num2str(MeanSY,'%.2f'),'; std= ',num2str(STDSY,'%.2f')];
        title({firstline;secondline})
        xlabel('\sigmaX & \sigmaY (nm)','FontSize',14');
        ylabel('frequency','FontSize',14');
        legend('\sigmaX','\sigmaY')%
        legend('boxoff')
        
        %---- Z vs (sigmaX/sigmaY), color = photons
        if sum(data(:,5))~=0
            subplot(2,4,8);
            scatter(data(:,5),SX./SY,3,data(:,13));
            c=colorbar;
            c.Label.String='Photons';
            c.Location='NorthOutside';
            xlabel('Z (nm)','FontSize',14');
            ylabel('\sigmaX/\sigmaY','FontSize',14');
        end
        
end

%% multiplot_stat_2D_cdata(data,dotsize,savefilename)
if size(data,2)==16
    data(:,17)=data(:,15)./data(:,16);
end

if nargin==3    
    %% Save Result
    SaveName=['histall_ch',num2str(ch),'_',savefilename];
    savefig(gcf,SaveName);
    print(SaveName,'-dpng')
    head={'statistic','Channel','Frame','X_nm','Y_nm','Z_nm',...
        'X_precision_nm','Y_precision_nm','X_pix','Y_pix','Z_pix',...
        'X_precision_pix','Y_precision_pix','Photons','Background',...
        'PSF_Sigma_X_pix','PSF_Sigma_Y_pix','sigmaRatio_sXoversY'};
    
    result(1,:)=mean(data,1);
    result(2,:)=std(data,1);
    result(3,:)=min(data);
    result(4,:)=median(data);
    result(5,:)=max(data);
    
    A=num2cell(result);
    colhead={'mean';'std';'min';'median';'max'};
    B=[colhead,A];
    T=cell2table(B,'VariableNames',head);
    writetable(T,[SaveName,'.csv']);
end
end

function[]  = stormimage2DColor_merge(locx0,locy0,locx1,locy1,range_x,range_y,dotsize,savefilename)
if isempty(range_x)
    xm=0.5; xp=427.5;
else
    Xrange=range_x;
    xm=Xrange(1)-.5; xp=Xrange(2)+.5;
end

if isempty(range_y)
    ym=0.5; yp=684.5;
else
    Yrange=range_y;
    ym=Yrange(1)-.5; yp=Yrange(2)+.5;
end

figure
scatter(locx0,locy0,dotsize,'bo','filled');
xlim([xm xp]);
ylim([ym yp]);
daspect([1 1 1])
xlabel('X (pix)')
ylabel('Y (pix)')
set(gca,'box','on')
hold on
scatter(locx1,locy1,dotsize,'mo','filled');
xlim([xm xp]);
ylim([ym yp]);
daspect([1 1 1])
xlabel('X (pix)')
ylabel('Y (pix)')
set(gca,'box','on')
legend('ch0','ch1','Location','northeastoutside')
savefig(gcf,[savefilename]);
end

function [c12_3D_tips,set_comb,radii]=tipbytipanalysis_3D_Paircorr(SaveName,Tiploc,data1,data2,tipws,Resolution,pixelsize,randomtest)


rMax=tipws/2;
pw=tipws/2;

for t=1:size(Tiploc,2)

datarange(1)= Tiploc(1,t) - pw; % yloc
datarange(2)= Tiploc(1,t) + pw;
datarange(3)= Tiploc(2,t) - pw; % xloc
datarange(4)= Tiploc(2,t) + pw;


% z Tiploc 
x1=data1.x;
y1=data1.y;

x2=data2.x;
y2=data2.y;

min_x=datarange(1);
max_x=datarange(2);
min_y=datarange(3);
max_y=datarange(4);

z1=data1.z;
z2=data2.z;

ind1=[];
ind1 =(x1>datarange(1) & x1<datarange(2)).* (y1>datarange(3) & y1<datarange(4));
z1_ws=z1(ind1==1);

ind2=[];
ind2 =(x2>datarange(1) & x2<datarange(2)).* (y2>datarange(3) & y2<datarange(4));
z2_ws=z2(ind2==1);

Tiploc(3,t)=mean(z1_ws);

N2_ws = sum(ind2);
x2_ws=x2(ind2==1);
y2_ws=y2(ind2==1);

datarange(5)= Tiploc(3,t) - pw;
datarange(6)= Tiploc(3,t) + pw;

min_z=datarange(5);
max_z=datarange(6);

ws=[(max_x-min_x),(max_y-min_y),(max_z-min_z)];


[c12_3D, radii,N1_ws,N2_ws]= pair_corr_cross_3D(data1, data2, datarange, Resolution, rMax, pixelsize);

c12_3D_tips(t,:)=c12_3D(1:end-1);



%%  randomize localization

if randomtest==1

   rand_data1.x= (max_x-min_x).*rand((N1_ws),1)+min_x;
   rand_data1.y= (max_y-min_y).*rand((N1_ws),1)+min_y; 
   rand_data1.z= (max_z-min_z).*rand((N1_ws),1)+min_z;


   rand_data2.x= (max_x-min_x).*rand((N2_ws),1)+min_x;
   rand_data2.y= (max_y-min_y).*rand((N2_ws),1)+min_y; 
   rand_data2.z= (max_z-min_z).*rand((N2_ws),1)+min_z;

   [c12_3D_Rand,radii]= pair_corr_cross_3D(rand_data1, rand_data2, datarange, Resolution, rMax, pixelsize);

   c12_3D_tips_rand(t,:)=c12_3D_Rand(1:end-1);

end

% 
end

Tiploc_3D=Tiploc;

set_comb.c12_tips= stat_combine_mat(c12_3D_tips,1,0,1);
set_comb.c12_tips_rand= stat_combine_mat(c12_3D_tips_rand,1,0,1);
set_comb.radii=radii;
set_comb.Tiploc_3D=Tiploc_3D;


figure
errorbar(radii(1:end-1),set_comb.c12_tips.mean,set_comb.c12_tips.ste);
hold on 
errorbar(radii(1:end-1),set_comb.c12_tips_rand.mean,set_comb.c12_tips_rand.ste);
legend ('data','random')
title ('pair-cross-correlation')
xlabel('r(nm)')
ylabel('g(r)')
xlim([0,radii(end-1)])

saveas(gcf,['Paircorr_3D_' ,SaveName]);
save(['Paircorr3D_',SaveName,'.mat'],'set_comb','-v7.3');

end

function [g,radii,N1_ws,N2_ws]= pair_corr_cross_3D(data1,data2, myrange, dr, rMax, pixelsize)

x1=data1.x;
y1=data1.y;
z1=data1.z;

x2=data2.x;
y2=data2.y;
z2=data2.z;


%%   RECT=[];

min_x=myrange(1);
max_x=myrange(2);
min_y=myrange(3);
max_y=myrange(4);
min_z=myrange(5);
max_z=myrange(6);

ws=[(max_x-min_x),(max_y-min_y),(max_z-min_z)];

dr =dr/pixelsize;
rad_bins = 0:dr:rMax; % pixelsize
in_radii = rad_bins(1:end-1); 
out_radii = rad_bins(2:end); 
mean_r = (in_radii+out_radii)/2; 
radii=mean_r*pixelsize;
volumns = (4 / 3)*pi*(out_radii.^3 - in_radii.^3); 

ind1 =(x1>myrange(1) & x1<myrange(2)).* (y1>myrange(3) & y1<myrange(4)).* (z1>myrange(5) & z1<myrange(6));

N1_ws = sum(ind1);
x1_ws=x1(ind1==1);
y1_ws=y1(ind1==1);
z1_ws=z1(ind1==1);


ind2 =(x2>myrange(1) & x2<myrange(2)).* (y2>myrange(3) & y2<myrange(4)).* (z2>myrange(5) & z2<myrange(6));

N2_ws = sum(ind2);
x2_ws=x2(ind2==1);
y2_ws=y2(ind2==1);
z2_ws=z2(ind2==1);


 if N1_ws < 1 
     error( " particle must be more than 1 ")
 else
     
     
    g = zeros(N2_ws, length(rad_bins)-1);   
    totalDensity = (N2_ws)/(ws(1)*ws(2)*ws(3));

         cnt=0;
     for p=1:N1_ws
         distance = sqrt((x1_ws(p) - x2_ws).^2 + (y1_ws(p) - y2_ws).^2 + (z1_ws(p) - z2_ws).^2);
         [result] = (hist(distance,rad_bins(1:end-1)));  
         g(p,:) = result/totalDensity;
     end
    

    g = mean(g,1)./volumns;

 end
 end 
