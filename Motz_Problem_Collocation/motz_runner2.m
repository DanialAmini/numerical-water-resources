clear;  format long;

M=40; N=10000; NN=4*N+3; A=zeros(NN,M); B=zeros(NN,1); V_xyt=zeros(NN,3);

for j=0:N
  V_xyt(j+1,1)=1;   V_xyt(j+1,2)=j/N;   t=atan(j/N);   V_xyt(j+1,3)=t;
  for n=0:M-1
    A(j+1,n+1)=((1/cos(t))^(n+1/2))*cos((n+1/2)*t)/sqrt(2.5)^(n+2);
  endfor
endfor
for j=0:N
  B(j+1)=500;
endfor

for j=0:2*N
  V_xyt(N+1+j+1,1)=(N-j)/N; V_xyt(N+1+j+1,2)=1; t=pi/2-atan((N-j)/N); V_xyt(N+1+j+1,3)=t;
  for n=0:M-1
    A(N+1+j+1,n+1)=-(n+1/2)*((1/sin(t))^(n-1/2))*sin((n-1/2)*t)/sqrt(2.5)^(n+2);
  endfor
endfor
for j=0:2*N
  B(N+1+j+1)=0;
endfor

for j=0:N
  V_xyt(3*N+2+j+1,1)=-1; V_xyt(3*N+2+j+1,2)=1-j/N; t=pi-atan(1-j/N); V_xyt(3*N+2+j+1,3)=t;
  for n=0:M-1
    A(3*N+2+j+1,n+1)=(n+1/2)*((-1/cos(t))^(n-1/2))*cos((n-1/2)*t)/sqrt(2.5)^(n+2);
  endfor  
endfor
for j=0:N
  B(3*N+2+j+1)=0;
endfor

%C=A\B
C=(transpose(A)*A)\(transpose(A)*B)

DD=zeros(M,1);
for n=0:M-1
  DD(n+1)=C(n+1)/sqrt(2.5)^(n+2);
endfor