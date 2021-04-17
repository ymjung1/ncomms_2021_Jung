function[combCell]=combine_multicells(celldata)
combMat=[];

for qq=1:size(celldata,2)
    for kk=1:size(celldata,1)        
        combMat=[combMat;cell2mat(celldata(kk,qq))];
         n(kk)= size(celldata{kk,qq},1);
    end
    
    if ~isequal(sum(n),size(combMat,1))
        error('check')
    end
    combCell{qq}=  combMat;
    n=[];
    combMat=[];
end


