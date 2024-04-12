
%%%%% genrating synthetic distributions under 625 parameter sets of negative feedback model

function negative_N100000

%clear
clc
close all

for zushu=1:1:625;
    zushu
tic
%%-------------------------------------
clearvars -except zushu
  
global  xdata tdata tt  mean  fano 


lat11=[0.3,0.7,1,2,4]; lat22 =[0.05,0.1,0.5,1,1.5]; gatt =[0.3,0.7,1,2,4];  v=[10,15,20,25,30]; %%%% 625 parameter sets

ind=fullfact([length(lat11) length(lat22) length(gatt) length(v)]);
kenen=[reshape(lat11(ind(:,1)),[],1) reshape(lat22(ind(:,2)),[],1) reshape(gatt(ind(:,3)),[],1) reshape(v(ind(:,4)),[],1)];%%所有可能
  
lam=kenen(zushu,1); mu=kenen(zushu,2); gamma=kenen(zushu,3);  v=kenen(zushu,4);

%%-------------------------------------
str1={num2str(lam),num2str(mu),num2str(gamma),num2str(v)}';

if zushu<=26
location=char(65+mod(zushu-1,625));
elseif zushu<=52 
location=char([65+mod(0,625), 65+mod(zushu-27,625)]);
elseif zushu<=78
    location=char([65+mod(1,625), 65+mod(zushu-53,625)]);
elseif  zushu<=104
    location=char([65+mod(2,625), 65+mod(zushu-79,625)]);
elseif  zushu<=130
    location=char([65+mod(3,625), 65+mod(zushu-105,625)]);
    elseif  zushu<=156
    location=char([65+mod(4,625), 65+mod(zushu-131,625)]);
    elseif  zushu<=182
    location=char([65+mod(5,625), 65+mod(zushu-157,625)]);
    elseif  zushu<=208
    location=char([65+mod(6,625), 65+mod(zushu-183,625)]);
    elseif  zushu<=234
    location=char([65+mod(7,625), 65+mod(zushu-209,625)]);
    elseif  zushu<=260
    location=char([65+mod(8,625), 65+mod(zushu-235,625)]);
    elseif  zushu<=286
    location=char([65+mod(9,625), 65+mod(zushu-261,625)]);
    elseif  zushu<=312
    location=char([65+mod(10,625), 65+mod(zushu-287,625)]);
    elseif  zushu<=338
    location=char([65+mod(11,625), 65+mod(zushu-313,625)]);
    elseif  zushu<=364
    location=char([65+mod(12,625), 65+mod(zushu-339,625)]);
    elseif  zushu<=390
    location=char([65+mod(13,625), 65+mod(zushu-365,625)]);
    elseif  zushu<=416
    location=char([65+mod(14,625), 65+mod(zushu-391,625)]);
    elseif  zushu<=442
    location=char([65+mod(15,625), 65+mod(zushu-417,625)]);
    elseif  zushu<=468
    location=char([65+mod(16,625), 65+mod(zushu-443,625)]);
    elseif  zushu<=494
    location=char([65+mod(17,625), 65+mod(zushu-469,625)]);
    elseif  zushu<=520
    location=char([65+mod(18,625), 65+mod(zushu-495,625)]);
    elseif  zushu<=546
    location=char([65+mod(19,625), 65+mod(zushu-521,625)]);
    elseif  zushu<=572
    location=char([65+mod(20,625), 65+mod(zushu-547,625)]);
    elseif  zushu<=598
    location=char([65+mod(21,625), 65+mod(zushu-573,625)]);
    elseif  zushu<=624
    location=char([65+mod(22,625), 65+mod(zushu-599,625)]);
else
    location=char([65+mod(23,625), 65+mod(zushu-625,625)]);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename='N100000_negative_10%.xlsx'; %%%% SSA for 100000 cell samples with 10% extrinsic noise on synthesis rate.
range1=[location,'1:',location,'4'];
xlswrite(filename,str1,1,range1)
%%-------------------------------------


N = 100000;  %% cell samples
delta=1; del=1;

vT=v*delta;
gaT=gamma*delta;
la1=lam*delta;
nu=mu*delta;

%%%%%%%%%%%%%%%%%%%%%%%%%%  Determine time point at steady state
ta=[];

TT=200;  H=200;  eps=0.01;
TP(1)=TT/H;
for i=2:1:H
   TP(i)=TT/H+TP(i-1) ;
end

%%%%%%%%%%%%%%%%%%
M=ceil(3*v+1); 
dt=0.001;   N=TT/dt+1;  % time step; max time;  # of time point
ta=[0:dt:TT]; 

p1=zeros(N,M+1); p2=zeros(N,M+1); 
  p1(1,1)=1; p2(1,1)=0;
for j=2:1:M+1
p2(1,j)=0;  
p1(1,j)=0;
end

for ta=2:N % time
    p1(ta,1)=p1(ta-1,1)+dt*(del*p1(ta-1,2)+gaT*p2(ta-1,1)-la1*p1(ta-1,1));
    p2(ta,1)=p2(ta-1,1)+dt*(del*p2(ta-1,2)+la1*p1(ta-1,1)-(vT+gaT)*p2(ta-1,1));
   
for i=2:M  % number
    p1(ta,i)=p1(ta-1,i)+dt*((gaT+(i-1)*nu)*p2(ta-1,i)-((i-1)*(del+0)+la1)*p1(ta-1,i)+i*del*p1(ta-1,i+1));
    p2(ta,i)=p2(ta-1,i)+dt*((la1+(i-1)*0)*p1(ta-1,i)-(vT+(i-1)*(del+nu)+gaT)*p2(ta-1,i)+i*del*p2(ta-1,i+1)+vT*p2(ta-1,i-1));
  end
   
end
pM=p1+p2; 
% 
% x=[0:1:M]; t=[0:dt:TT];

% plot(x,pM(20./dt+1,:),'-'); 
% hold on

for j=1:1:H
   meandata(j)=0; 
   seconddata(j)=0;
end

for j=1:1:H
    for i=1:1:M+1
meandata(j)=(i-1)*pM(TP(j)./dt+1,i)+meandata(j);
    end
end

for j=1:1:H
    for i=1:1:M+1
seconddata(j)=(i-1)^2*pM(TP(j)./dt+1,i)+seconddata(j);
    end
end

TPP(1)=0; meandataa(1)=0; seconddataa(1)=0;
for j=2:1:H+1
TPP(j)=TP(j-1);
meandataa(j)=meandata(j-1);
seconddataa(j)=seconddata(j-1);
end

k=H; 
for j=1:1:H
  
   if  abs(meandataa(H)-meandataa(H-j))/meandataa(H)<=eps
      k=k-1;
   else
       break
   end
end
meantime=k*(TT/H);

k1=H; 
for j=1:1:H
   if  abs(seconddataa(H)-seconddataa(H-j))/seconddataa(H)<eps
      k1=k1-1;
   else
       break
   end
end
secondtime=k1*(TT/H);

T=3*max(meantime,secondtime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SSA 
M=3*v+1;        
X(1) = 0; t(1) = 0;  % initial number and time
t1= 1*T/5 ;  
t2 =2*T/5 ; 
t3=3*T/5 ; 
t4=4*T/5 ; 
t5=T; %% steady state

S1 = []; S2 = []; S3=[]; S4=[]; S5=[];
        
for i = 1:N
 s(1)=1;
    n = 1; 
    Xall = [X(1)]; Sall = [s(1)]; tall = [t(1)];
    
      v; % mean of synthesis rate 
        sta=0.1*v; % standard deviation being equal 0.1 (presenting 10% extrinsic noise)
        eee = log((v^2)/sqrt(sta+v^2)); 
        sigma = sqrt(log(sta/(v^2)+1)); % log-normal distributed random variable
        v1=lognrnd(eee,sigma); 
    
  while t(n) <= T
        h1 = 1; c1 = v1; a1 = h1*c1; % generate
        h2 = 1; h4 = X(n); c2 = gamma; a2 = h2*c2+h4*nu; % ON --> I1
        h3 = 1; c3 = lam ; a3 = h3*c3;  % I1 --> ON 
        c4 = delta; a4 = h4*c4;  % decay
        n = n+1;
        
        if s(n-1) == 0 % ON
            a0 = a1+a2+a4;
            r1 = rand; r2=rand;
            tau = -log(r1)/a0;  % time interval in which nothing occurs
            if a0*r2 <= a1 % generate occurs
                X(n)=X(n-1)+1;  s(n)=0;
            elseif a0*r2 <=a1+a2  % transition occurs
                X(n)=X(n-1);  s(n)=1;
            else % decay
                X(n)=X(n-1)-1; s(n)=0;
            end
            t(n) = t(n-1) + tau;
        else s(n-1) == 1; % I1
            a0 = a3+a4;
            r1=rand;  r2=rand;
            tau=-log(r1)/a0;
            if a0*r2 <= a3 % transition occurs
                X(n)=X(n-1); s(n)=0;
            else % decy occurs
                X(n)=X(n-1)-1; s(n)=1;
            end
            t(n) = t(n-1) + tau;
        end
        Xall = [Xall X(n)]; Sall = [Sall s(n)]; tall = [tall t(n)];
    end
    if t(n)<T
        Xall=[Xall X(n)]; Sall=[Sall s(n)]; tall=[tall T];
    end
    for j=2:n
        if t(j-1) <=t1 && t(j) >t1
            S1 = [S1 X(j-1)];
        end
        if t(j-1) <=t2 && t(j) >t2
            S2 = [S2 X(j-1)];
        end
        if t(j-1) <=t3 && t(j) >t3
            S3 = [S3 X(j-1)];
        end
        if t(j-1) <=t4 && t(j) >t4
            S4 = [S4 X(j-1)];
        end
        if t(j-1) <=t5 && t(j) >t5
            S5 = [S5 X(j-1)];
        end
    end
   
end



for i = 1:M
    y(i) = sum(S5==(i-1))/N;
end

xdata=y;
tt=N.*xdata;
S=numel(xdata);
for i=1:1:S;
tdata(i)=i-1;
end
xdata;

mean=0;
for i=1:1:S
mean=mean+(i-1).*xdata(i);    
end    
mean;

Var=0;
for i=1:1:S
    Var=Var+((i-1)^2).*xdata(i); 
end
fano=Var/mean-mean;
%%-------------------------------------
str2={num2str(mean),num2str(fano)}';
range2=[location,'5:',location,'6'];
xlswrite(filename,str2,1,range2)
range3=[location,'7'];
 xlswrite(filename,xdata',1,range3)
toc
end
