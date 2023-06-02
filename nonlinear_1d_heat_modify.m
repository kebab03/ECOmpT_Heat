
%https://arxiv.org/ftp/arxiv/papers/1811/1811.06337.pdf

function main
rho=1; Cp=1; kappa0=0.1; chi=0.5;
M=31;N=2% N=41;
tEnd=15; tau=tEnd/(M-1);
a=1; b=3; h=(b-a)/(N-1);
alpha=2; beta=1;
x=zeros(N,1);
u0=zeros(N,1);
for i=1:N
 x(i)=a+(i-1)*h;
 u0(i)=2-(x(i)-1)/2+(x(i)-1)*(x(i)-3);
 
%figure(i);
 plot(i,u0(i),'o')
 
end
t=zeros(M,1);
for n=1:M
 t(n)=(n-1)*tau;
end
u_1=zeros(N,1); u=zeros(N,1);
uNext=zeros(N,1); G=zeros(N,1);
L=zeros(N,N); L(1,1)=1; L(N,N)=1;
U=zeros(N,M); U(:,1)=u0;
for n=2:M
 u=U(:,n-1); u_1=U(:,n-1);
 eps=1;
 while(eps>0.0001)
 G(1)=u(1)-alpha;
 G(N)=u(N)-beta;
 for i=2:N-1
 k=kappa0*exp(chi*u(i));
 Dk=chi*k;
 D2k=chi*Dk;
 v=(u(i+1)-u(i-1))/(2*h);
 A=rho*Cp/tau;
 phi=A*(u(i)-u_1(i))-Dk*v*v;
 f=phi/k;
 q=(-f*Dk+A-D2k*v*v)/k;
 p=-2*Dk*v/k;
 G(i)=u(i+1)-2*u(i)+u(i-1)-h*h*f;
 L(i,i-1)=1+0.5*h*p;
 L(i,i)=-2-h*h*q;
 L(i,i+1)=1-0.5*h*p;
 end
 uNext=u-L\G;
 eps=sqrt(h*(uNext-u)'*(uNext-u));
 u=uNext;
 end
 U(:,n)=u;
 end
 mesh(x,t,U');
end
