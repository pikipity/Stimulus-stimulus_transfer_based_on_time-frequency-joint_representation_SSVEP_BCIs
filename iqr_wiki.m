function [num_outlier,outlier,L,U,Q1,Q2,Q3,IQR]=iqr_wiki(x)
    Q2=median(sort(x));
    small_media_tmp=x(x<Q2);
    large_media_tmp=x(x>Q2);
    Q1=median(sort(small_media_tmp));
    Q3=median(sort(large_media_tmp));
    IQR=Q3-Q1;
    L=Q1-1.5*IQR;
    U=Q3+1.5*IQR;
    outlier_ind=[find(x<=L) find(x>=U)];
    outlier=x(outlier_ind);
    num_outlier=length(outlier_ind);
end