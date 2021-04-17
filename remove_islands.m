function [SegL_new] = remove_islands(segbw,SegL)

seg_bw=~isnan(segbw);
SE = strel('disk',5);
seg_bw=imclose(seg_bw,SE);
CC = bwconncomp(seg_bw);
numPixels = cellfun(@numel,CC.PixelIdxList);

for n=1:length(numPixels)
    
    if  numPixels(n)<=2000
        seg_bw(CC.PixelIdxList{n}) =0;
    end
end

BW=SegL==2;
BW=(BW+seg_bw);
BW=imfill(BW,'holes');
%imagesc(BW)
BW1=BW-seg_bw;
%imagesc(BW1)
BW1=imfill(logical(BW1),'holes');
CC = bwconncomp(BW1);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);

BW2=SegL==1;
SegL_new=zeros(size(BW));
SegL_new(CC.PixelIdxList{idx}) = 2;
SegL_new=(BW2+SegL_new);
SegL_new(SegL_new==3)=2;
SegL_new=SegL_new+(seg_bw.*(SegL==3).*3);
end

