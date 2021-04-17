function [BWout,InBackgroundBW,InEdgeBW,OutBackgroundBW] = segmentation(image,BWimage)

%-----< Remove extracellular background >-----%
BW = imbinarize(image, 'adaptive','Sensitivity',1);
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);

if biggest <= ((size(image,1)*size(image,2))*0.2) % 10% of entire image
    while  1
        se = strel('disk',3);
        BW=imdilate(BW,se);
        BW=imfill(logical(BW),'holes');
        CC = bwconncomp(BW);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [biggest,idx] = max(numPixels);
        if biggest > ((size(image,1)*size(image,2))*0.2)% && secondbiggest < ((size(image,1)*size(image,2))*0.05)
            break
            disp(biggest)
        end
    end
elseif biggest >= (size(image,1)*size(image,2))*0.90
    while  1
        se = strel('disk',10);
        BW=imerode(BW,se);
        BW=imfill(logical(BW),'holes');
        CC = bwconncomp(BW);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [biggest,idx] = max(numPixels);
        if biggest <((size(image,1)*size(image,2))*0.9) && biggest >((size(image,1)*size(image,2))*0.1) %&& secondbiggest < ((size(image,1)*size(image,2))*0.05)
            break
            disp(biggest)
        elseif  biggest <= ((size(image,1)*size(image,2))*0.1)
            while  1
                se = strel('disk',3);
                BW=imdilate(BW,se);
                BW=imfill(logical(BW),'holes');
                CC = bwconncomp(BW);
                numPixels = cellfun(@numel,CC.PixelIdxList);
                [~,idx] = max(numPixels);
            end
        end
    end
else
    BW=imfill(logical(BW),'holes');
    CC = bwconncomp(BW);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    
end

OutBackgroundBW=ones(size(BW));
OutBackgroundBW(CC.PixelIdxList{idx}) = 0;
OutCleanedBW=(~OutBackgroundBW.*BWimage);

%-----< Remove intracellular background >-----%
se = strel('disk',60);
erodedBW = imdilate(OutBackgroundBW , se);
se = strel('disk',90);
erodedBW = imdilate(erodedBW , se);
CleanedBW=OutCleanedBW.*(erodedBW);
se = strel('disk',5);
CleanedBW=imclose(CleanedBW,se);
CC = bwconncomp(CleanedBW);
numPixels = cellfun(@numel,CC.PixelIdxList);

for n=1:length(numPixels)
    if  numPixels(n)<100
        CleanedBW(CC.PixelIdxList{n}) = 0;
    end
end
BW=(~OutBackgroundBW);
se = strel('cube',10);
BW=imerode(BW,se);
OutBackgroundBW=~BW;

%-----< get intracellular background area :InBackgroundBW >-----%
BW=(OutBackgroundBW+CleanedBW);
se = strel('cube',10);
BW=imdilate(BW,se);
BW=~BW;
BW=imerode(BW,se);
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
InBackgroundBW=zeros(size(BW));
InBackgroundBW(CC.PixelIdxList{idx}) = 1;
InBackgroundBW=imfill(InBackgroundBW,'holes');
InBackgroundBW=imdilate(InBackgroundBW,se);
InBackgroundBW=imdilate(InBackgroundBW,se);
SegBWClean=~InBackgroundBW.*~OutBackgroundBW.*BWimage;
CC = bwconncomp(SegBWClean);
numPixels = cellfun(@numel,CC.PixelIdxList);
for n=1:length(numPixels)
    
    if  numPixels(n)<500
        SegBWClean(CC.PixelIdxList{n}) = 0;
    end
end

InBackgroundBW=activecontour(SegBWClean,InBackgroundBW,25,'edge');
BW=InBackgroundBW+SegBWClean;
bw1=imfill(BW,'holes');
BW=(bw1-BW);
BW=BW+InBackgroundBW;
BW=imfill(BW,'holes');
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
InBackgroundBW=zeros(size(BW));
InBackgroundBW(CC.PixelIdxList{idx}) = 1;

se = strel('disk',5);
InBackgroundBW=imdilate(InBackgroundBW,se);
BW=(InBackgroundBW+CleanedBW);
BW=imfill(BW,'holes');
OutBackgroundBW=~BW;

%-----< get the segmented area :BWout>-----%
BWout=~InBackgroundBW.*~OutBackgroundBW.*BWimage;
%-----< get the cell inner edge :InEdgeBW>-----%
InEdgeBW = edge(InBackgroundBW,'sobel',0,'nothinning');
end
