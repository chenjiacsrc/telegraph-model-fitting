
% Fitting cross-talk pathway model to telegraph model 
function endtest_crosstalk

%clear
clc
close all


str={'v','la1','la2','q1','gat','bar_v','bar_la','bar_gat','f','H','KL','mean','fano'};
% v:synthesis rate; la1: activation rate of weak pathway; la2: activation rate of strong pathway; q1: selection probability of weak pathway; 
%gat: inactivation rate; bar_v: effective synthesis rate; bar_la: effective activation rate; bar_gat: effective inactivation rate
xlswrite('crosstalk_N100000_10%.xlsx',str,1,'A1');
%%-------------------------------------
clearvars -except zushu

data=xlsread('N100000_crosstalk_10%.xlsx');%% read synthetic distributions under 625 parameter sets of cross-talk pathway model and 100000 cell samples with 10% extinsic noise on synthesis rate.


for zushu=1:1:625;
%%-------------------------------------
  clearvars -except zushu data filename
  
global  xdata tdata tt  mean  fano  

N=100000;
L=70;

xdata=data(8:L,zushu);
real_para=data(1:7,zushu); %distribution data 

xdata=clearnan(xdata);
lam1=real_para(1);
lam2=real_para(2);
gamma=real_para(3);
q1=real_para(4);
q2=1-q1;
v=real_para(5);
m=real_para(6);
fano=real_para(7);

print=[v,lam1,lam2,q1,gamma];

S=length(xdata);

tdata=[];
%%----------------------------
for i=1:1:S;
tdata(i)=i-1;
tt(i)=N.*xdata(i);
end
% minimization step
optims.MaxIter = 10000;
optims.MaxFunEvals = 10000;

f=m*gamma/(v-m); %%% average frequency of gene activation

%%%%%%% initial value for likelyhood method
k=length(xdata);

for i=1:1:length(xdata)
if xdata(k)>0
      break
    else
        k=length(xdata)-i;
end
end
vin=k-1; 
lain=m*(vin+1-fano-m)/(vin*(fano-1)); 
gain=(vin-m)*lain/m;  

%%%%%%%%%%%%%%%%% maximum liklyhood method

b2=[log(lain)  log(gain)  log(vin)];  %%初值  la ga v

[b2min, S2min] = fminsearch(@S11fun_pm1,b2); % 调用搜索函数

rho =exp(b2min(3)); d = 1; sigma1 = exp(b2min(1)); sigma0 = exp(b2min(2));
a = sigma1/d; b = (sigma0+sigma1)/d; u = rho/d;
N = 2*(max(tdata));

% FSP
num = 2*N;
Q = zeros(num);
one = ones(num,1);
for i = 1:N-1
    Q(i+1,i) = i*d;
    Q(i+N+1,i+N) = i*d;
    Q(i+N,i+N+1) = rho;
end
for i = 1:N
    Q(i,i+N) = sigma1;
    Q(i+N,i) = sigma0;
end
temp = Q*one;
for i = 1:num
    Q(i,i) = -temp(i);
end
Qmod = Q;
for i = 1:num
    Qmod(i,num) = 1;
end
vec = zeros(1,num);
vec(num) = 1;
ssd = vec/Qmod;

% tmpoff=0; tmpon=0;
% for t=1:1:N
%     tmpoff=tmpoff+ssd(t);
% end
% 
% for t=N+1:1:2*N
%     tmpon=tmpon+ssd(t);
% end

tmdist = ssd(1:N)+ssd(N+1:2*N);
% tmdistoff=ssd(1:N)./tmpoff;
% tmdiston=ssd(N+1:2*N)./tmpon;
 
for t = 1:1:max(tdata)+1;
 tmx(t)=tmdist(t);
% tmxoff(t)=tmdistoff(t);
%  tmxon(t)=tmdiston(t);
end

% figure(1)
% plot(tdata2,tmx,'k','linewidth',2);
% hold on
% 
% figure(2)
% plot(tdata2,tmxoff,'k','linewidth',2);
% hold on
% 
% figure(3)
% plot(tdata2,tmxon,'k','linewidth',2);
% hold on
% %%%%%%%%%%%%%%%%%%%%%%%%%%5
% f2=tsmpon;
% f22=tmpon;

%%%%%%%%%%%%%%%%%%%%%%%%%HD
dM=length(tdata);

H1=0;
for i=1:1:dM;
    H1=H1+(sqrt(tmx(i))-sqrt(xdata(i))).^2; 
end
H=1/sqrt(2)*sqrt(H1);

% H1off=0;
% for i=1:1:dM;
%     H1off=H1off+(sqrt(tmxoff(i))-sqrt(tsmxoff(i))).^2;
% end
% Hoff=1/sqrt(2)*sqrt(H1off);
% 
% H1on=0;
% for i=1:1:dM;
%     H1on=H1on+(sqrt(tmxon(i))-sqrt(tsmxon(i))).^2;
% end
% Hon=1/sqrt(2)*sqrt(H1on);

KL=0;
for  i=1:1:dM;
KL=KL+xdata(i)*log((xdata(i)/tmx(i)));
end


% KLoff=0;
% for  i=1:1:dM;
% KLoff=KLoff+tsmxoff(i)*log((tsmxoff(i)/tmxoff(i)));
% end
% 
% KLon=0;
% for  i=1:1:dM;
% KLon=KLon+tsmxon(i)*log((tsmxon(i)/tmxon(i)));
% end

%----------------------------------------
print=[print,exp(b2min(3)), exp(b2min(1)), exp(b2min(2)), f,H,KL,m, fano];
xlswrite('crosstalk_N100000_10%.xlsx',print,1,['A',num2str(zushu+1)]);




end

end

