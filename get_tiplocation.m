
function [tiploc]=get_tiplocation(InboundEdge,SegL,MVlengthlimit_min,MVlengthlimit_max,closedistlimit)


MVBW=(SegL==3);
se = strel('disk',10);
BW=imclose(MVBW,se);
se=strel('cube',5);
BW=imerode(BW,se);
CH = bwconvhull(BW);
se=strel('cube',3);
CH=imerode(CH,se);
CH=~CH;
% imshowpair(CH,BW);
%title('Union Convex Hull');
%imshow(CH)
Tips=BW.*CH;
% imshow(Tips)
se = strel('cube',3);
Tips=imdilate(Tips,se);
[L,n]=bwlabel(Tips,8); 
Tips=imerode(Tips,se);
L=L.*Tips;

tiploc=[];
for k=1:n
    [r,c]=find(L==k);
    tiploc(1,k)=mean(r);
    tiploc(2,k)=mean(c);
end


if ~isempty(tiploc)
    %% calculate distance from inbouwndedge
    [r1,c1]=find(InboundEdge==1);
    Y=[c1,r1];
    
    D=[];
    InxEdge=[];
    shortestD=[];
    for m=1:size(tiploc,2)
        X=[tiploc(2,m),tiploc(1,m)];
        [D,I] = pdist2(Y,X,'euclidean','Smallest',1);
        InxEdge(1,m)=Y(I,(2));
        InxEdge(2,m)=Y(I,(1));
        shortestD(m)=D;
    end
    
    %% remove too short or too long MVs
    
    removeind=find((shortestD < MVlengthlimit_min) | (shortestD > MVlengthlimit_max));
    if ~isempty(removeind)
      
        InxEdge(:,removeind)=[];
        tiploc(:,removeind)=[];
        shortestD(removeind)=[];

    end
    clear D I
    removeind=[];
    
    %% remove too close-tips
    
    Y1=[tiploc(2,:);tiploc(1,:)]';
    for m=1:size(tiploc,2)
        X1=[tiploc(2,m),tiploc(1,m)];
        [D,I] = pdist2(Y1,X1,'euclidean','smallest',2);% find the closest tip position
        
        if  length(D)==2 && D(2) <= closedistlimit
            tip1dist=shortestD(m);
            tip2dist=shortestD(I(2));
            
            if tip1dist<=tip2dist
                removeind(m)=m;
            else
                removeind(m)=I(2);
            end
        else
            removeind(m)=0;
            
        end
        D=[]; I=[];
    end
    
    removeind=unique(removeind);
    removeind=removeind(removeind>0);
    
    if ~isempty(removeind) && ~isempty(tiploc)
        tiploc(:,removeind)=[];
    end
    
end
end