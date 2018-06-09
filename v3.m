G=load('img1.txt');H=load('img2.txt');%提取原始TXT文件中数据
G1=G(:,2:3);H1=H(:,2:3);
F=[G1,H1];%文件中的坐标顺序应该是x1,y1,x2,y2
[m1,n1]=size(F);%m为点对数目

x1=F(1:end,1);
y1=F(1:end,2);
x2=F(1:end,3);
y2=F(1:end,4);
  
bx=abs(x1(1)-x2(1));
f=153.124;%焦距
% 相对定向元素赋初值
d=[0 0 0 0 0]';
%给相对定向元素改正数附初值
Xd=[1 1 1 1 1]';


num=0;
 
while (abs(Xd(1))>0.00003|| abs(Xd(2))>0.00003|| abs(Xd(3))>0.00003|| abs(Xd(4))>0.00003|| abs(Xd(5))>0.00003)
    
num=num+1;
    

%确定左右旋转矩阵R1和R2,并计算各点的像空间辅助坐标
R1=eye(3);
R2=zeros(3,3);%旋转矩阵R2
 
R2(1,1)=cos(d(1))*cos(d(3))-sin(d(1))*sin(d(2))*sin(d(3));
R2(1,2)=-cos(d(1))*sin(d(3))-sin(d(1))*sin(d(2))*sin(d(3));
R2(1,3)=-sin(d(1))*cos(d(2));
R2(2,1)=cos(d(2))*sin(d(3));
R2(2,2)=cos(d(2))*cos(d(3));
R2(2,3)=-sin(d(2));
R2(3,1)=sin(d(1))*cos(d(3))+cos(d(1))*sin(d(2))*sin(d(3));
R2(3,2)=-sin(d(1))*sin(d(3))+cos(d(1))*sin(d(2))*cos(d(3));
R2(3,3)=cos(d(1))*cos(d(2));
 
X1=x1(1);
Y1=y1(1);
Z1=-f;
X2=R2(1,1)*x2(1)+R2(1,2)*y2(1)-R2(1,3)*f;
Y2=R2(2,1)*x2(1)+R2(2,2)*y2(1)-R2(2,3)*f;
Z2=R2(3,1)*x2(1)+R2(3,2)*y2(1)-R2(3,3)*f;
 
by=bx*d(4);     
bz=bx*d(5);


N1=(bx*Z2-bz*X2)/(X1*Z2-X2*Z1);%N1为初始左投影系数，N2为初始右投影系数
N2=(bx*Z1-bz*X1)/(X1*Z2-X2*Z1);
%计算误差方程的各项系数
a11=-X2*Y2*N2/Z2;
a12=-(Z2+Y2*Y2/Z2)*N2;
a13=X2*N2;
a14=bx;
a15=-Y2*bx/Z2;
 
A=[a11,a12,a13,a14,a15];
L=N1*Y1-N2*Y2-by;
for i=2:m1
    
    X1=x1(i);
    Y1=y1(i);
    Z1=-f;
    X2=R2(1,1)*x2(i)+R2(1,2)*y2(i)-R2(1,3)*f;
    Y2=R2(2,1)*x2(i)+R2(2,2)*y2(i)-R2(2,3)*f;
    Z2=R2(3,1)*x2(i)+R2(3,2)*y2(i)-R2(3,3)*f;
   
    N1=(bx*Z2-bz*X2)/(X1*Z2-X2*Z1);
    N2=(bx*Z1-bz*X1)/(X1*Z2-X2*Z1);
    
    a11=-X2*Y2*N2/Z2;
    a12=-(Z2+Y2*Y2/Z2)*N2;
    a13=X2*N2;
    a14=bx;
    a15=-Y2*bx/Z2;
    
    A1=[a11,a12,a13,a14,a15];
    L1=N1*Y1-N2*Y2-by;
    A=[A;A1];
    L=[L;L1];
end
%求解法方程，求出相对定向元素的改正数
    Xd=(A'*A)\(A'*L);
    d=d+Xd;
    by=bx*d(4);
    bz=bx*d(5);
 
   
    
end



Q=inv(A'*A);%法方程的系数矩阵的逆
V=A*Xd-L;
m0=sqrt((V'*V)/(m1-5));%单位权中误差
B=[Q(1,1),Q(2,2),Q(3,3),Q(4,4),Q(5,5)];
m2=sqrt(B)*m0;%各未知数的中误差

fp=fopen('相对定向结果.txt','wt'); 
fprintf(fp,'迭代次数:%d\n',num); 

fprintf(fp,'\n相对定向元素结果：\n\n');
fprintf(fp,'by=%g\t',by);fprintf(fp,'bz=%g\t',bz);fprintf(fp,'Φ=%g\n',d(1));
fprintf(fp,'ω=%g\t',d(2));fprintf(fp,'κ=%g\t\n',d(3));


fprintf(fp,'\n单位权中误差的值：\n\nmo=%g\n\n',m0);

fprintf(fp,'\n相对定向元素中误差为：\n\n');
fprintf(fp,'mu=%g\t',B(4));fprintf(fp,'mγ=%g\t',B(5));fprintf(fp,'mΦ=%g\n',B(1));
fprintf(fp,'mω=%g\t',B(2));fprintf(fp,'mκ=%g\t',B(3));
 