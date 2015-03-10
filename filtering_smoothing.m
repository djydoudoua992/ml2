%compute p(x_1,x_2,...,x_t) or p(x_1,x_2,...,x_T)
function p_x=modelprobability(n,t,Y,E,pi,A,Mode)
a=hmm_alpha(n,t,Y,E,pi,A,Mode);
p_x=a(1)+a(2);
end

%compute p(z_t|x_1,x_2,...,x_t)
function p_zx=filtering(n,t,Y,E,pi,A,Mode)
p_zx=zeros(2,1);
a=hmm_alpha(n,t,Y,E,pi,A,Mode);
p_zx(1)=a(1)/sum(a);
p_zx(2)=a(2)/sum(a);
end

%compute p(z_t|x_1,x_2,...,x_T)
function p_zx=smoothing(n,t,T,Y,E,pi,A,Mode)
b=hmm_beta(n,t,T,Y,E,A,Mode);
a=hmm_alpha(n,t,Y,E,pi,A,Mode);
p_zx=zeros(2,1);
p_zx(1)=a(1)*b(1)/sum(a);
p_zx(2)=a(2)*b(2)/sum(a);
end