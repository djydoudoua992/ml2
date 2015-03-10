%compute alpha function
function p_alpha=hmm_alpha(n,t,Y,E,pi,A,Mode)
p_alpha=zeros(2,1);
if t==1
    p_xz=emission(n,1,Y,E,Mode);
    p_alpha(1)=pi(1)*p_xz(1);
    p_alpha(2)=pi(2)*p_xz(2);

elseif t>1
    p_xz=emission(n,t,Y,E,Mode);
    alpha2=hmm_alpha(n,t-1,Y,E,pi,A,Mode);%recursively calculate alpha
    p_alpha(1)=p_xz(1)*(A(1,1)*alpha2(1)+A(2,1)*alpha2(2));
    p_alpha(2)=p_xz(2)*(A(1,2)*alpha2(1)+A(2,2)*alpha2(2));
end    
end

%compute beta function
function p_beta=hmm_beta(n,t,T,Y,E,A,Mode)
p_beta=zeros(2,1);
if t==T
    p_beta=ones(2,1);

elseif t<T
    p_xz=emission(n,t,Y,E,Mode);
    beta2=hmm_beta(n,t+1,T,Y,E,A,Mode);%recursively calculate beta
    p_beta(1)=p_xz(1)*A(1,1)*beta2(1)+p_xz(2)*A(1,2)*beta2(2);
    p_beta(2)=p_xz(1)*A(2,1)*beta2(1)+p_xz(2)*A(2,2)*beta2(2);
end 
end

%compute p(x_t|z_t)
function p_xz=emission(n,t,Y,E,Mode)
p_xz=zeros(2,1);
if strcmp(Mode, 'discrete')
    p_xz(1)=E(1,Y(n,t));
    p_xz(2)=E(2,Y(n,t));
elseif strcmp(Mode, 'continuous')
    m=E.mu(1);
    si=E.sigma2(1);
    x=Y(n,t);
    p_xz(1)=exp(-(x-m)*(x-m)/(2*si*si))/(si*sqrt(2*3.141592654));
    m=E.mu(2);
    si=E.sigma2(2);
    p_xz(2)=exp(-(x-m)*(x-m)/(2*si*si))/(si*sqrt(2*3.141592654));
end    
end