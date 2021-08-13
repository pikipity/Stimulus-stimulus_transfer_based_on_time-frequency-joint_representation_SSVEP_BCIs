function y=gen_ref_sin(f,fs,L,N,phase)
    t=linspace(0,(L-1)/fs,L);
    y=[];
    for n=1:N
        y=[y;...
            sin(2.*pi.*n.*f.*t+n.*phase);
            cos(2.*pi.*n.*f.*t+n.*phase)];
    end
end