function CI95=cal_CI95(x,varargin)
    if isempty(varargin)
        ci=0.95;
    else
        ci=varargin{1};
    end
    N=size(x,1);
    SEM = std(x) / sqrt(N); 
    CI95 = SEM * tinv(ci, N-1);
end