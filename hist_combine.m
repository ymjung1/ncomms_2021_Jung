function [statresult] = hist_combine(HistCellData,sortingoption,pixelsize,Exfactor)

for p=1:size(HistCellData,2)
    datalength(p) =  size(HistCellData{p},2);
    B=max(datalength);
end

HistDist=nan(p,B);
for f=1:size(HistCellData,1)
for p=1:size(HistCellData,2)
    
    A=HistCellData{f,p};
    if strcmp(sortingoption, 'RISvalue')
    A=A.*pixelsize./Exfactor;
    end
    
    statresult.mean(f,p)=mean(A);
    statresult.std(f,p)=std(A);
    statresult.ste(f,p)=std(A)./sqrt(length(A));
    statresult.median(f,p)=median(A);
    statresult.quatile25(f,p)=quantile(A,[0.25]);
    statresult.quatile50(f,p)=quantile(A,[0.5]);
    statresult.quatile75(f,p)=quantile(A,[0.75]);
    A=[];
   
end
end


    statresult.mean_all=mean(statresult.mean,2,'omitnan');
    statresult.std_all=std(statresult.mean,[],2,'omitnan');
     
    sq=(sqrt(f-length(find(all(isnan(statresult.mean),1)))));
    statresult.ste_all= statresult.std_all./sq ;

end



