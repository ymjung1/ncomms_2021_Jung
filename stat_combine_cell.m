function [statresult] = stat_combine_cell(CellData,dimention)

% dimention :1 vertical cat
% dimention :2 horizantal cat
 
if dimention==1
    k=2;
elseif dimention==2
    k=1;
else
    error('error')
end

for p=1:size(CellData,k)
    
    if p==1 
        A=CellData{p};
    else
        
        A=vertcat(A,CellData{p});
    end

end

    statresult.combined_data=A;
    statresult.n=size(A);
    statresult.mean=mean(A,dimention,'omitnan');
    statresult.std=std(A,[],dimention,'omitnan');
    statresult.ste=std(A)./sqrt(length(A));
    statresult.median=median(A,dimention,'omitnan');
    statresult.quatile25=quantile(A,[0.25],dimention);
    statresult.quatile50=quantile(A,[0.5],dimention);
    statresult.quatile75=quantile(A,[0.75],dimention);
  
end




