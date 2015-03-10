%perform EM
%for discrete observation, the results of em_mu and em_sigma are zeros
%for continuous observation, the result of em_bj2 is zero
function [em_pi,em_A,em_bj2,em_mu,em_sigma]=EM(N,T,Y,E,pi,A,Mode)
if nargin < 7
    Mode = 'discrete';
end

%compute E[z_1k] and pi_k
em_pi=zeros(2,1);
e1=zeros(2,1);
for n=1:N
    b=hmm_beta(n,1,T,Y,E,A,Mode);
    a=hmm_alpha(n,1,Y,E,pi,A,Mode);
    p_x=modelprobability(n,T,Y,E,pi,A,Mode);
    e1(1)=e1(1)+b(1)*a(1)/p_x;
    e1(2)=e1(2)+b(2)*a(2)/p_x;
end
em_pi(1)=e1(1)/sum(e1);
em_pi(2)=e1(2)/sum(e1);

%compute E[z_tk] and bj2, eg. j=6
em_bj2=0;
if strcmp(Mode, 'discrete')
et=zeros(2,1);
ex=0;
for n=1:N
    p_x=modelprobability(n,T,Y,E,pi,A,Mode);
    for t=1:T
        b=hmm_beta(n,t,T,Y,E,A,Mode);
        a=hmm_alpha(n,t,Y,E,pi,A,Mode);        
        %et(1)=et(1)+b(1)*a(1)/p_x;%for calculating b_j1
        et(2)=et(2)+b(2)*a(2)/p_x;
        x=Y(n,t);
        if x==6%change j for b_j2 where j=1..6
            ex=ex+et(2)*x;%add E[z_t2]*x_tj
        end
    end
end
em_bj2=ex/et(2);

em_mu=zeros(2,1);
em_sigma=zeros(2,1);
elseif strcmp(Mode, 'continuous')
    et=zeros(2,1);
    ex=zeros(2,1);
    exm=zeros(2,1);
    for n=1:N
        p_x=modelprobability(n,T,Y,E,pi,A,Mode);
        for t=1:T
            b=hmm_beta(n,t,T,Y,E,A,Mode);
            a=hmm_alpha(n,t,Y,E,pi,A,Mode);        
            et(1)=et(1)+b(1)*a(1)/p_x;
            et(2)=et(2)+b(2)*a(2)/p_x;
            x=Y(n,t);         
            ex=ex+et*x;%add E[z_tk]*x_t
            exm(1)=exm(1)+et(1)*(x-E.mu(1))*(x-E.mu(1));
            exm(2)=exm(2)+et(2)*(x-E.mu(2))*(x-E.mu(2));
        end
    end
    em_mu=ex./et;
    em_sigma=exm./et;
end    

%compute E[z_t-1j z_tk] and a_jk
em_A=zeros(2);
%if strcmp(Mode, 'discrete')
    et=zeros(2);
    for n=1:N
        p_x=modelprobability(n,T,Y,E,pi,A,Mode);
        for t=2:T
            b=hmm_beta(n,t,T,Y,E,A,Mode);
            a=hmm_alpha(n,t-1,Y,E,pi,A,Mode);
            p_xz=emission(n,t,Y,E,Mode);
            et(1,1)=et(1,1)+a(1)*p_xz(1)*A(1,1)*b(1)/p_x;
            et(1,2)=et(1,2)+a(1)*p_xz(2)*A(1,2)*b(2)/p_x;
            et(2,1)=et(2,1)+a(2)*p_xz(1)*A(2,1)*b(1)/p_x;
            et(2,2)=et(2,2)+a(2)*p_xz(2)*A(2,2)*b(2)/p_x;
        end
    end
    em_A(1,1)=et(1,1)/(et(1,1)+et(1,2));
    em_A(1,2)=et(1,2)/(et(1,1)+et(1,2));
    em_A(2,1)=et(2,1)/(et(2,1)+et(2,2));
    em_A(2,2)=et(2,2)/(et(2,1)+et(2,2));
%end
end
