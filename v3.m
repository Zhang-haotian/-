G=load('img1.txt');H=load('img2.txt');%��ȡԭʼTXT�ļ�������
G1=G(:,2:3);H1=H(:,2:3);
F=[G1,H1];%�ļ��е�����˳��Ӧ����x1,y1,x2,y2
[m1,n1]=size(F);%mΪ�����Ŀ

x1=F(1:end,1);
y1=F(1:end,2);
x2=F(1:end,3);
y2=F(1:end,4);
  
bx=abs(x1(1)-x2(1));
f=153.124;%����
% ��Զ���Ԫ�ظ���ֵ
d=[0 0 0 0 0]';
%����Զ���Ԫ�ظ���������ֵ
Xd=[1 1 1 1 1]';


num=0;
 
while (abs(Xd(1))>0.00003|| abs(Xd(2))>0.00003|| abs(Xd(3))>0.00003|| abs(Xd(4))>0.00003|| abs(Xd(5))>0.00003)
    
num=num+1;
    

%ȷ��������ת����R1��R2,������������ռ丨������
R1=eye(3);
R2=zeros(3,3);%��ת����R2
 
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


N1=(bx*Z2-bz*X2)/(X1*Z2-X2*Z1);%N1Ϊ��ʼ��ͶӰϵ����N2Ϊ��ʼ��ͶӰϵ��
N2=(bx*Z1-bz*X1)/(X1*Z2-X2*Z1);
%�������̵ĸ���ϵ��
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
%��ⷨ���̣������Զ���Ԫ�صĸ�����
    Xd=(A'*A)\(A'*L);
    d=d+Xd;
    by=bx*d(4);
    bz=bx*d(5);
 
   
    
end



Q=inv(A'*A);%�����̵�ϵ���������
V=A*Xd-L;
m0=sqrt((V'*V)/(m1-5));%��λȨ�����
B=[Q(1,1),Q(2,2),Q(3,3),Q(4,4),Q(5,5)];
m2=sqrt(B)*m0;%��δ֪���������

fp=fopen('��Զ�����.txt','wt'); 
fprintf(fp,'��������:%d\n',num); 

fprintf(fp,'\n��Զ���Ԫ�ؽ����\n\n');
fprintf(fp,'by=%g\t',by);fprintf(fp,'bz=%g\t',bz);fprintf(fp,'��=%g\n',d(1));
fprintf(fp,'��=%g\t',d(2));fprintf(fp,'��=%g\t\n',d(3));


fprintf(fp,'\n��λȨ������ֵ��\n\nmo=%g\n\n',m0);

fprintf(fp,'\n��Զ���Ԫ�������Ϊ��\n\n');
fprintf(fp,'mu=%g\t',B(4));fprintf(fp,'m��=%g\t',B(5));fprintf(fp,'m��=%g\n',B(1));
fprintf(fp,'m��=%g\t',B(2));fprintf(fp,'m��=%g\t',B(3));
 