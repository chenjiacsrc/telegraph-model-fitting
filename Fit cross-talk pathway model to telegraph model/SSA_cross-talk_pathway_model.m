%%%%% genrating synthetic distributions under 625 parameter sets of cross-talk pathway model

function crosstalk_N100000

%clear
clc
close all

for zushu=1:1:625
    zushu
tic
%%-------------------------------------
clearvars -except zushu
  
lat11=0.2; lat22 =[0.5,1,2,4,8]; q1=[0.1,0.3,0.5,0.7,0.9]; gatt =[0.5,1,2,4,8];  v=[10,15,20,25,30]; 

ind=fullfact([length(lat22) length(q1) length(gatt) length(v)]);
kenen=[reshape(lat22(ind(:,1)),[],1) reshape(q1(ind(:,2)),[],1)  reshape(gatt(ind(:,3)),[],1) reshape(v(ind(:,4)),[],1)];%%所有可能
lam1=lat11;
lam2=kenen(zushu,1); 
q1=kenen(zushu,2);
q2=1-q1;
gamma=kenen(zushu,3);  
nu=kenen(zushu,4);

global  xdata tdata tt  mean  fano 


str1={num2str(lam1),num2str(lam2),num2str(gamma),num2str(q1),num2str(nu)}';

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
    
filename='N100000_crosstalk_10%.xlsx';  %%% SSA for 100000 cell samples with 10% extrinsic noise on synthesis rate.
range1=[location,'1:',location,'5'];
xlswrite(filename,str1,1,range1)
%%-------------------------------------


N = 100000;   %% cell samples

delta = 1;  
vT=nu*delta;
gaT=gamma*delta;
la1=lam1*delta;
la2=lam2*delta;

ta=[]; yy=[]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%  Determine time point at steady state

TT=200;  H=200;  eps=0.01;
TP(1)=TT/H;
for i=2:1:H
   TP(i)=TT/H+TP(i-1) ;
end

odefun=@(ta,yy)[
q1*gaT*yy(3)-la1*yy(1);
q2*gaT*yy(3)-la2*yy(2);
la1*yy(1)+la2*yy(2)-gaT*yy(3);
vT*yy(3)-delta*yy(4);
-(gaT+delta)*yy(5)+la1*yy(6)+la2*yy(7)+vT*yy(3);
q1*gaT*yy(5)-(la1+delta)*yy(6);
q2*gaT*yy(5)-(la2+delta)*yy(7);
-2*delta*yy(8)+delta*yy(4)+vT*yy(3)+2*vT*yy(5)
];

y00=[q1,q2,0,0,0,0,0,0];
options=odeset('reltol',1e-6,'abstol',1e-8);
[ta,yy]=ode45(odefun,[0,TT],y00,options);

for j=1:1:H
   for i=1:1:length(ta)
   if ta(i)>=TP(j)
   break
   end
   k(j)=i+1;
   end
end

for j=1:1:H
meandata(j)=yy(k(j),4);
seconddata(j)=yy(k(j),8);
end

k=H; 
for j=1:1:H
  
   if  abs(meandata(H)-meandata(H-j))/meandata(H)<=eps
      k=k-1;
   else
       break
   end
end
meantime=k*(TT/H);

k1=H; 
for j=1:1:H
   if  abs(seconddata(H)-seconddata(H-j))/seconddata(H)<eps
      k1=k1-1;
   else
       break
   end
end
secondtime=k1*(TT/H);

T=6*max(meantime,secondtime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SSA 

 %T = 36; 
M=3*nu+1;
X(1) = 0; t(1) = 0;  % % initial number and time

t1= 0; t2 =T/10; t3=2*T/10; t4=3*T/10; t5=4*T/10;t6=5*T/10; 
t7=6*T/10; t8=7*T/10; t9=8*T/10; t10=9*T/10;
t11=T; %% steady state

S1 = []; S2 = []; S3=[]; S4=[]; S5=[];S6 = []; 
S7 = []; S8=[]; S9=[]; S10=[];S11 = [];


for i = 1:N
    if rand<=q1 
        s(1)=1; 
    else
        s(1)=2;
    end
    n = 1; 
    Xall = [X(1)]; Sall = [s(1)]; tall = [t(1)];
    
     nu;  % mean of synthesis rate 
        sta=0.1*nu; % standard deviation being equal 0.1 (presenting 10% extrinsic noise)
        mu = log((nu^2)/sqrt(sta+nu^2));
        sigma = sqrt(log(sta/(nu^2)+1));  % log-normal distributed random variable 
        nu1=lognrnd(mu,sigma); 
    
    while t(n) <= T
        h1 = 1; c1 = nu1; a1 = h1*c1; % generate
        h2 = 1; c2 = gamma; a2 = h2*c2; % ON --> 
        h3 = 1; c3 = lam1; a3 = h3*c3;  % I1 --> ON
        h4 = 1; c4 = lam2; a4 = h4*c4;  % I2 --> ON 
        h5 = X(n); c5 = delta; a5 = h5*c5;  % decay
        n = n+1;
        if s(n-1) == 0 % ON
            a0 = a1+a2+a5;
            r1 = rand; r2=rand;
            tau = -log(r1)/a0;  % time interval in which nothing occurs
            if a0*r2 <= a1 % generate occurs
                X(n)=X(n-1)+1;  s(n)=0;
            elseif a0*r2 <=a1+a2  % transition occurs
                X(n)=X(n-1);
                if rand <=q1
                    s(n)=1;
                else
                    s(n)=2;
                end
            else % decay
                X(n)=X(n-1)-1; s(n)=0;
            end
            t(n) = t(n-1) + tau;
        elseif s(n-1) == 1 % I1
            a0 = a3+a5;
            r1=rand;  r2=rand;
            tau=-log(r1)/a0;
            if a0*r2 <= a3 % transition occurs
                X(n)=X(n-1); s(n)=0;
            else % decy occurs
                X(n)=X(n-1)-1; s(n)=1;
            end
            t(n) = t(n-1) + tau;
        else  % I2
            a0 = a4+a5;
            r1=rand;  r2=rand;
            tau=-log(r1)/a0;
            if a0*r2 <= a4  % transition occurs
                X(n)=X(n-1); s(n)=0;
            else  % decay occurs
                X(n)=X(n-1)-1; s(n)=2;
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
        if t(j-1) <=t6 && t(j) >t6
            S6 = [S6 X(j-1)];
        end
        if t(j-1) <=t7 && t(j) >t7
            S7 = [S7 X(j-1)];
        end
        if t(j-1) <=t8 && t(j) >t8
            S8 = [S8 X(j-1)];
        end
        if t(j-1) <=t9 && t(j) >t9
            S9 = [S9 X(j-1)];
        end
        if t(j-1) <=t10 && t(j) >t10
            S10 = [S10 X(j-1)];
        end
        if t(j-1) <=t11 && t(j) >t11
            S11 = [S11 X(j-1)];
        end
    end
    
end
    

for i = 1:M
    y(i) = sum(S11==(i-1))./N;
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
range2=[location,'6:',location,'7'];
xlswrite(filename,str2,1,range2)
range3=[location,'8'];
 xlswrite(filename,xdata',1,range3)
toc
end

