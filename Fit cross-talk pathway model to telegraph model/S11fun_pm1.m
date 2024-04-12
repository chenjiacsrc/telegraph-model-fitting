function S1 = S11fun_pm1(b2)
%%%%极大似然


global tdata xdata  tt

rho =exp(b2(3)); d = 1; sigma1 = exp(b2(1)); sigma0 = exp(b2(2));
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
dist = ssd(1:N)+ssd(N+1:2*N);
 
for t = 1:1:max(tdata)+1;
      x(t)=dist(t);
end



% for t = 1:1:max(tdata)+1;
% x(t) = ((b(3)).^(t-1)./factorial(t-1)).*(gamma(b(1)+b(2)).*gamma(b(1)+(t-1)))./(gamma(b(1)).*gamma(b(1)+b(2)+(t-1))).*hypergeom([b(1)+(t-1)],[b(1)+b(2)+(t-1)],-b(3));
% 
% end


%% 数值计算

%% plot result of the integration
figure(2)
plot(tdata,xdata,'r.','MarkerSize',20);
hold on
t = 0:1:max(tdata);
plot(t,x,'k','linewidth',2);
 hold off
drawnow

xpred = interp1(t,x,tdata);

S1 = 0;
for i = 1:length(tdata)
    S1 = S1 - tt(i).*log(xpred(i));
end



end