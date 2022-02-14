clear;
format long;

k=0;
kk=0;
M=10;
N=500;

temp1=zeros(M,M);
temp2=zeros(M,M);
temp3=zeros(M,M);

temp0=zeros(M,M);
tempR=zeros(M,1);
res_=zeros(M,1);

vec1=zeros(N+1,M);
vec2=zeros(N+1,M);
vec3=zeros(N+1,M);

vecT=zeros(N+1,3);

for k=0:M-1
    for j=0:N
      
      y=j/N;  x=1; r=sqrt(x^2+y^2);  t=atan(y/x);  vecT(j+1,1)=t;
      vec1(j+1,k+1)=cos((k+1/2)*t)*r^(k+1/2)/sqrt(2.5)^(k+2);
      
      x=1-2*j/N;  y=1; r=sqrt(x^2+y^2);  t=pi/2-atan(x/y); vecT(j+1,2)=t;
      vec2(j+1,k+1)=-(k+1/2)*sin((k-1/2)*t)*r^(k-1/2)/sqrt(2.5)^(k+2);
      
      y=1-j/N;  x=-1; r=sqrt(x^2+y^2);  t=pi+atan(y/x); vecT(j+1,3)=t;
      vec3(j+1,k+1)=(k+1/2)*cos((k-1/2)*t)*r^(k-1/2)/sqrt(2.5)^(k+2);
      
    endfor
endfor

for k=0:M-1
  for kk=k:M-1
    temp1(k+1,kk+1)=transpose(vec1(1:N+1,k+1))*vec1(1:N+1,kk+1)*(1/N);
    temp2(k+1,kk+1)=transpose(vec2(1:N+1,k+1))*vec2(1:N+1,kk+1)*(2/N);
    temp3(k+1,kk+1)=transpose(vec3(1:N+1,k+1))*vec3(1:N+1,kk+1)*(1/N);
    
    temp1(kk+1,k+1)=temp1(k+1,kk+1);
    temp2(kk+1,k+1)=temp2(k+1,kk+1);
    temp3(kk+1,k+1)=temp3(k+1,kk+1);
  endfor
endfor

temp0=temp1+1/N*(temp2+temp3);

for k=0:M-1
  tempR(k+1,1)=500*sum(vec1(1:N+1,k+1))*(1/N);
endfor

res_=temp0\tempR;

res2_=res_;

for k=0:M-1
  res2_(k+1,1)=res2_(k+1,1)/sqrt(2.5)^(k+2);
endfor




d_exact=zeros(40,1);

d_exact(1)=401.162453745234;
d_exact(2)=87.6559201950879;
d_exact(3)=17.2379150794467;
d_exact(4)=-8.07121525968247;
d_exact(5)=1.44027271702287;
d_exact(6)=0.331054885920767;
d_exact(7)=0.275437344508727;
d_exact(8)=-0.0869329946589323;
d_exact(9)=0.0336048786193139;
d_exact(10)=0.0153843742840408;
d_exact(11)=0.00730230172821564;
d_exact(12)=-0.00318411377492366;
d_exact(13)=0.00122064638572558;
d_exact(14)=0.000530965296814512;
d_exact(15)=0.000271512155856615;
d_exact(16)=-0.000120045575740674;
d_exact(17)=5.05397300846273E-05;
d_exact(18)=0.000023167016399;
d_exact(19)=1.15353890475415E-05;
d_exact(20)=-5.29397510035331E-06;
d_exact(21)=2.29074482358489E-06;
d_exact(22)=1.06369298555841E-06;
d_exact(23)=5.31498056296845E-07;
d_exact(24)=-2.45403558772047E-07;
d_exact(25)=1.09425820081517E-07;
d_exact(26)=5.20412323411083E-08;
d_exact(27)=2.59454810138248E-08;
d_exact(28)=-1.09175391022716E-08;
d_exact(29)=5.27436083354048E-09;
d_exact(30)=2.73360939552139E-09;
d_exact(31)=1.35407691721933E-09;
d_exact(32)=-1.81685168223046E-10;
d_exact(33)=2.26648038099555E-10;
d_exact(34)=1.54146776907919E-10;
d_exact(35)=7.39636752201555E-11;
d_exact(36)=5.05824030696331E-11;
d_exact(37)=6.1655071915812E-12;
d_exact(38)=6.54434622972211E-12;
d_exact(39)=2.94394253116419E-12;
d_exact(40)=5.12785385295806E-12;


%%%%%%%%%%%%%%%%%%%%%% testing
u_R=zeros(N+1,3);

for j=0:N
  x=1;y=j/N; r=sqrt(x^2+y^2); t=atan(y/x);
  for k=0:40-1
    u_R(j+1,1)=u_R(j+1,1)+d_exact(k+1)*cos((k+1/2)*t)*r^(k+1/2);
  endfor
  for k=0:M-1
    u_R(j+1,2)=u_R(j+1,2)+res2_(k+1)*cos((k+1/2)*t)*r^(k+1/2);
  endfor  
endfor

for j=0:N
  u_R(j+1,3)=abs(u_R(j+1,1)-u_R(j+1,2));
endfor
disp("u_R");
max(u_R(:,3))
mean(u_R(:,3))

%below 0.004% error for u_R


%%%%%%%%%%%%%%%%%%%%%% testing
u_zT=zeros(N+1,3);

for j=0:N
  y=1;x=1-2*j/N; r=sqrt(x^2+y^2); t=pi/2-atan(x/y);
  for k=0:40-1
    u_zT(j+1,1)=u_zT(j+1,1)+d_exact(k+1)*(-1)*(k+1/2)*sin((k-1/2)*t)*r^(k-1/2);
  endfor
  for k=0:M-1
    u_zT(j+1,2)=u_zT(j+1,2)+res2_(k+1)*(-1)*(k+1/2)*sin((k-1/2)*t)*r^(k-1/2);
  endfor  
endfor

for j=0:N
  u_zT(j+1,3)=abs(u_zT(j+1,1)-u_zT(j+1,2));
endfor
disp("u_zT");
max(u_zT(:,3)) %
mean(u_zT(:,3))%

%below 0.08% error for u_zT




%%%%%%%%%%%%%%%%%%%%%% testing
u_xR=zeros(N+1,3);

for j=0:N
  x=-1; y=1-j/N; r=sqrt(x^2+y^2); t=pi-atan(y/x);
  for k=0:40-1
    u_xR(j+1,1)=u_xR(j+1,1)+d_exact(k+1)*(k+1/2)*cos((k-1/2)*t)*r^(k-1/2);
  endfor
  for k=0:M-1
    u_xR(j+1,2)=u_xR(j+1,2)+res2_(k+1)*(k+1/2)*cos((k-1/2)*t)*r^(k-1/2);
  endfor  
endfor

for j=0:N
  u_xR(j+1,3)=abs(u_xR(j+1,1)-u_xR(j+1,2));
endfor
disp("u_xR");
max(u_xR(:,3)) 
mean(u_xR(:,3))

%below 0.2% error for u_xR



%%%%%%%%%%%%%%%%%%%%%% testing
u_domain1=zeros(N+1,N+1);
u_domain2=zeros(N+1,N+1);
u_domain3=zeros(N+1,N+1);

for j=0:N
for jj=0:N  
  x=1-2*j/N;y=jj/N; r=sqrt(x^2+y^2); 
  
  
  t=0  ;  
  #for k=0:0-1
    u_domain1(j+1,jj+1)=t;
  #endfor
  #for k=0:M-1
    u_domain2(j+1,jj+1)=t;
  #endfor  
endfor
endfor
