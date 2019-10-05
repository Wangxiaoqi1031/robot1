clc
clear all
syms a1 a2 a3 a4 a5 a6 d1 d2 d3 d4 d5 d6 q1 q2 q3 q4 q5 q6 p1 p2 p3 p4 p5 p6 lx mx nx px ly my ny py lz mz nz pz 
%D-H parameter
a1=0.25;a2=0.95;a3=0.3;a4=0.55;a5=0.16;a6=0.22;d1=0.19;d2=0.37;d3=0.1;d4=1.55;d5=0.21;d6=0.13;
p1=20;p2=30;p3=-45;p4=80;p5=-120;p6=100;
%Randomly select a set of joint variables
q1_i=14.5427;q2_i=-44.5435;q3_i=130.5472;q4_i=-70.2145;q5_i=60.2147;q6_i=160.9745;
%Positive kinematics
c1=cosd(q1_i); c2=cosd(q2_i); c3=cosd(q3_i); c4=cosd(q4_i); c5=cosd(q5_i); c6=cosd(q6_i);
s1=sind(q1_i); s2=sind(q2_i); s3=sind(q3_i); s4=sind(q4_i); s5=sind(q5_i); s6=sind(q6_i);
u1=cosd(p1); u2=cosd(p2); u3=cosd(p3); u4=cosd(p4); u5=cosd(p5); u6=cosd(p6);
y1=sind(p1); y2=sind(p2); y3=sind(p3); y4=sind(p4); y5=sind(p5); y6=sind(p6);
A1=[c1 -s1*y1 s1*u1 a1*c1;s1 c1*y1 -c1*u1 a1*s1; 0 u1 y1 d1;0 0 0 1];
A2=[c2 -s2*y2 s2*u2 a2*c2;s2 c2*y2 -c2*u2 a2*s2; 0 u2 y2 d2;0 0 0 1];
A3=[c3 -s3*y3 s3*u3 a3*c3;s3 c3*y3 -c3*u3 a3*s3; 0 u3 y3 d3;0 0 0 1];
A4=[c4 -s4*y4 s4*u4 a4*c4;s4 c4*y4 -c4*u4 a4*s4; 0 u4 y4 d4;0 0 0 1];
A5=[c5 -s5*y5 s5*u5 a5*c5;s5 c5*y5 -c5*u5 a5*s5; 0 u5 y5 d5;0 0 0 1];
A6=[c6 -s6*y6 s6*u6 a6*c6;s6 c6*y6 -c6*u6 a6*s6; 0 u6 y6 d6;0 0 0 1];
A=A1*A2*A3*A4*A5*A6
lx=A(1,1);mx=A(1,2);nx=A(1,3);px=A(1,4);
ly=A(2,1);my=A(2,2);ny=A(2,3);py=A(2,4);
lz=A(3,1);mz=A(3,2);nz=A(3,3);pz=A(3,4);
%Inverse kinematics
u= mx*u6+nx*y6;
v= my*u6+ny*y6;
w= mz*u6+nz*y6;
p= -lx*a6-(mx*u6+nx*y6)*d6+px;
q= -ly*a6-(my*u6+ny*y6)*d6+py;
y= -lz*a6-(mz*u6+nz*y6)*d6+pz;
Q1=p*p+q*q+(y-d1)*(y-d1)-a1*a1;
Q2=2*(p*u+q*v+(y-d1)*w);
z=p*u1^2*v - q*u*u1^2;
%Omit complex calculations,directly calculate elements of A, B, and C matrices
k11_1=- 2*a5*d2*u2*y4 - (2*a1*a5*p*u*u1*u2*y4)/z - (2*a1*a5*q*u1*u2*v*y4)/z;
k11_2=- 2*a2*a5*y4 - (2*a1*p^2*u1*u2*u5*y3)/z - (2*a1*q^2*u1*u2*u5*y3)/z;
k11_3=- 2*a3*a5*y4 - (2*a1*p^2*u1*u3*u5*y2)/z - (2*a1*q^2*u1*u3*u5*y2)/z;
k12_1= - 2*a2*a5*y3 - (2*a1*p^2*u1*u2*u5*y4)/z - (2*a1*q^2*u1*u2*u5*y4)/z;
k12_2=2*a5*d2*u2*y3 + (2*a1*a5*p*u*u1*u2*y3)/z + (2*a1*a5*q*u1*u2*v*y3)/z;
k12_3=2*a5*d3*u3 + 2*a5*d2*u3*y2 + (2*a1*a5*p*u*u1*u3*y2)/z + (2*a1*a5*q*u1*u3*v*y2)/z;
k13_1= - 2*a2*a5*y3*y4 - (2*a1*p^2*u1*u2*u5)/z - (2*a1*q^2*u1*u2*u5)/z;
k13_2=2*a5*d2*u2*y3*y4 + (2*a1*a5*p*u*u1*u2*y3*y4)/z + (2*a1*a5*q*u1*u2*v*y3*y4)/z;
k13_3=2*a5*d3*u3*y4 + 2*a5*d2*u3*y2*y4 + (2*a1*a5*p*u*u1*u3*y2*y4)/z + (2*a1*a5*q*u1*u3*v*y2*y4)/z;
k14_1= 2*a5*d2*u2 + (2*a1*a5*p*u*u1*u2)/z + (2*a1*a5*q*u1*u2*v)/z;
k14_2=(2*a1*u1*u2*u5*y3*y4*p^2)/z + (2*a1*u1*u2*u5*y3*y4*q^2)/z + 2*a2*a5;
k14_3=(2*a1*u1*u3*u5*y2*y4*p^2)/z + (2*a1*u1*u3*u5*y2*y4*q^2)/z + 2*a3*a5;
k15_1= 2*d2*d5*u2*u4 - 2*a2*a4*y3 - (2*a1*p^2*u1*u2*u4*y5)/z - (2*a1*q^2*u1*u2*u4*y5)/z + (2*a1*d5*p*u*u1*u2*u4)/z + (2*a1*d5*q*u1*u2*u4*v)/z;
k15_2=2*a2*d5*u4 + 2*a4*d2*u2*y3 + (2*a1*a4*p*u*u1*u2*y3)/z + (2*a1*a4*q*u1*u2*v*y3)/z;
k15_3=2*a4*d3*u3 + 2*a3*d5*u4 + 2*a4*d2*u3*y2 + (2*a1*a4*p*u*u1*u3*y2)/z + (2*a1*a4*q*u1*u3*v*y2)/z;
k16_1= 2*a4*d2*u2 + 2*a2*d5*u4*y3 + (2*a1*a4*p*u*u1*u2)/z + (2*a1*a4*q*u1*u2*v)/z;
k16_2=(2*a1*u1*u2*u4*y3*y5*p^2)/z - (2*a1*d5*u*u1*u2*u4*y3*p)/z + (2*a1*u1*u2*u4*y3*y5*q^2)/z - (2*a1*d5*u1*u2*u4*v*y3*q)/z + 2*a2*a4 - 2*d2*d5*u2*u4*y3;
k16_3=(2*a1*u1*u3*u4*y2*y5*p^2)/z - (2*a1*d5*u*u1*u3*u4*y2*p)/z + (2*a1*u1*u3*u4*y2*y5*q^2)/z - (2*a1*d5*u1*u3*u4*v*y2*q)/z + 2*a3*a4 - 2*d3*d5*u3*u4 - 2*d2*d5*u3*u4*y2;
k17_1=2*a2*a5*u3*u4;
k17_2=- 2*a5*d2*u2*u3*u4 - (2*a1*a5*p*u*u1*u2*u3*u4)/z - (2*a1*a5*q*u1*u2*u3*u4*v)/z;
k17_3= 2*a5*d4*u4 + 2*a5*d3*u4*y3 + 2*a5*d2*u4*y2*y3 + (2*a1*a5*p*u*u1*u4*y2*y3)/z + (2*a1*a5*q*u1*u4*v*y2*y3)/z;
k18_1=0;
k18_2=- (2*a1*p^2*u1*u2*u3*u4*u5)/z - (2*a1*q^2*u1*u2*u3*u4*u5)/z;
k18_3=(2*a1*u1*u4*u5*y2*y3*p^2)/z + (2*a1*u1*u4*u5*y2*y3*q^2)/z + 2*a4*a5;
k19_1=2*a3*d2*u2 + 2*a2*d4*u3 + 2*a2*d5*u3*y4 + (2*a1*a3*p*u*u1*u2)/z + (2*a1*a3*q*u1*u2*v)/z;
k19_2=2*a2*a3 - 2*d2*d4*u2*u3 - 2*d2*d5*u2*u3*y4 - (2*a1*d4*p*u*u1*u2*u3)/z - (2*a1*d4*q*u1*u2*u3*v)/z + (2*a1*p^2*u1*u2*u3*y4*y5)/z + (2*a1*q^2*u1*u2*u3*y4*y5)/z - (2*a1*d5*p*u*u1*u2*u3*y4)/z - (2*a1*d5*q*u1*u2*u3*v*y4)/z;
k19_3=2*d1*y - a1^2 + a2^2 + a3^2 + a4^2 + a5^2 - d1^2 + d2^2 + d3^2 + d4^2 + d5^2 - p^2 - q^2 - y^2 + 2*d2*d3*y2 + 2*d3*d4*y3 + 2*d4*d5*y4 + 2*d2*d4*y2*y3 + 2*d3*d5*y3*y4 + 2*d2*d5*y2*y3*y4 + (2*a1*d2*p*u*u1)/z + (2*a1*d2*q*u1*v)/z + (2*a1*p^2*u1*w*y1)/z + (2*a1*q^2*u1*w*y1)/z + (2*a1*d1*p*u*u1*y1)/z + (2*a1*d3*p*u*u1*y2)/z + (2*a1*d1*q*u1*v*y1)/z + (2*a1*d3*q*u1*v*y2)/z - (2*a1*p*u*u1*y*y1)/z - (2*a1*q*u1*v*y*y1)/z + (2*a1*d4*p*u*u1*y2*y3)/z + (2*a1*d4*q*u1*v*y2*y3)/z - (2*a1*p^2*u1*y2*y3*y4*y5)/z - (2*a1*q^2*u1*y2*y3*y4*y5)/z + (2*a1*d5*p*u*u1*y2*y3*y4)/z + (2*a1*d5*q*u1*v*y2*y3*y4)/z;
k21_1=- a2*u5*y3 - (a1*a5*u^2*u1*u2*y4)/z - (a1*a5*u1*u2*v^2*y4)/z;
k21_2= d2*u2*u5*y3 - (a1*p*u*u1*u2*u5*y3)/z - (a1*q*u1*u2*u5*v*y3)/z;
k21_3= d3*u3*u5 + d2*u3*u5*y2 - (a1*p*u*u1*u3*u5*y2)/z - (a1*q*u1*u3*u5*v*y2)/z;
k22_1= d2*u2*u5*y4 - (a1*p*u*u1*u2*u5*y4)/z - (a1*q*u1*u2*u5*v*y4)/z;
k22_2=(a1*a5*u1*u2*y3*u^2)/z + (a1*a5*u1*u2*y3*v^2)/z + a2*u5*y4;
k22_3=(a1*a5*u1*u3*y2*u^2)/z + (a1*a5*u1*u3*y2*v^2)/z + a3*u5*y4;
k23_1= d2*u2*u5 - (a1*p*u*u1*u2*u5)/z - (a1*q*u1*u2*u5*v)/z;
k23_2=(a1*a5*u1*u2*y3*y4*u^2)/z + (a1*a5*u1*u2*y3*y4*v^2)/z + a2*u5;
k23_3= (a1*a5*u1*u3*y2*y4*u^2)/z + (a1*a5*u1*u3*y2*y4*v^2)/z + a3*u5;
k24_1=(a1*a5*u1*u2*u^2)/z + (a1*a5*u1*u2*v^2)/z + a2*u5*y3*y4;
k24_2=(a1*p*u*u1*u2*u5*y3*y4)/z - d2*u2*u5*y3*y4 + (a1*q*u1*u2*u5*v*y3*y4)/z;
k24_3=(a1*p*u*u1*u3*u5*y2*y4)/z - d2*u3*u5*y2*y4 - d3*u3*u5*y4 + (a1*q*u1*u3*u5*v*y2*y4)/z;
k25_1= (a1*d5*u1*u2*u4*u^2)/z - (a1*p*u1*u2*u4*y5*u)/z + (a1*d5*u1*u2*u4*v^2)/z - (a1*q*u1*u2*u4*y5*v)/z + d2*u2*u4*y5;
k25_2=(a1*a4*u1*u2*y3*u^2)/z + (a1*a4*u1*u2*y3*v^2)/z + a2*u4*y5;
k25_3=(a1*a4*u1*u3*y2*u^2)/z + (a1*a4*u1*u3*y2*v^2)/z + a3*u4*y5;
k26_1=(a1*a4*u1*u2*u^2)/z + (a1*a4*u1*u2*v^2)/z + a2*u4*y3*y5;
k26_2= (a1*p*u*u1*u2*u4*y3*y5)/z - (a1*d5*u^2*u1*u2*u4*y3)/z - (a1*d5*u1*u2*u4*v^2*y3)/z - d2*u2*u4*y3*y5 + (a1*q*u1*u2*u4*v*y3*y5)/z;
k26_3= (a1*p*u*u1*u3*u4*y2*y5)/z - d2*u3*u4*y2*y5 - (a1*d5*u^2*u1*u3*u4*y2)/z - (a1*d5*u1*u3*u4*v^2*y2)/z - d3*u3*u4*y5 + (a1*q*u1*u3*u4*v*y2*y5)/z;
k27_1=0;
k27_2=- (a1*a5*u^2*u1*u2*u3*u4)/z - (a1*a5*u1*u2*u3*u4*v^2)/z;
k27_3=(a1*a5*u1*u4*y2*y3*u^2)/z + (a1*a5*u1*u4*y2*y3*v^2)/z + a4*u5;
k28_1=-a2*u3*u4*u5;
k28_2= d2*u2*u3*u4*u5 - (a1*p*u*u1*u2*u3*u4*u5)/z - (a1*q*u1*u2*u3*u4*u5*v)/z;
k28_3=(a1*p*u*u1*u4*u5*y2*y3)/z - d3*u4*u5*y3 - d2*u4*u5*y2*y3 - d4*u4*u5 + (a1*q*u1*u4*u5*v*y2*y3)/z;
k29_1=(a1*a3*u1*u2*u^2)/z + (a1*a3*u1*u2*v^2)/z + a2*u3*y4*y5;
k29_2=(a1*p*u*u1*u2*u3*y4*y5)/z - (a1*d4*u^2*u1*u2*u3)/z - (a1*d4*u1*u2*u3*v^2)/z - (a1*d5*u^2*u1*u2*u3*y4)/z - (a1*d5*u1*u2*u3*v^2*y4)/z - d2*u2*u3*y4*y5 + (a1*q*u1*u2*u3*v*y4*y5)/z;
k29_3=d1*w + d5*y5 - p*u - q*v - w*y + d4*y4*y5 + d3*y3*y4*y5 + (a1*d2*u^2*u1)/z + (a1*d2*u1*v^2)/z + d2*y2*y3*y4*y5 + (a1*d1*u^2*u1*y1)/z + (a1*d3*u^2*u1*y2)/z + (a1*d1*u1*v^2*y1)/z + (a1*d3*u1*v^2*y2)/z - (a1*u^2*u1*y*y1)/z - (a1*u1*v^2*y*y1)/z + (a1*d4*u^2*u1*y2*y3)/z + (a1*d4*u1*v^2*y2*y3)/z + (a1*p*u*u1*w*y1)/z + (a1*q*u1*v*w*y1)/z + (a1*d5*u^2*u1*y2*y3*y4)/z + (a1*d5*u1*v^2*y2*y3*y4)/z - (a1*p*u*u1*y2*y3*y4*y5)/z - (a1*q*u1*v*y2*y3*y4*y5)/z;
k31_1=d3*u2*u5*y3 + d4*u2*u5+ d5*u2*u5*y4 - (a5*d1*u^2*u1^2*u2*y4)/z - (a5*d1*u1^2*u2*v^2*y4)/z + (a5*u^2*u1^2*u2*y*y4)/z + (a5*u1^2*u2*v^2*y*y4)/z - (a5*p*u*u1^2*u2*w*y4)/z - (a5*q*u1^2*u2*v*w*y4)/z;
k31_2=a3*u2*u3*u5 - a2*u5*y2*y3 - a5*u2*y3*y5 - (p^2*u1^2*u2*u5*w*y3)/z - (q^2*u1^2*u2*u5*w*y3)/z - (d1*p*u*u1^2*u2*u5*y3)/z - (d1*q*u1^2*u2*u5*v*y3)/z + (p*u*u1^2*u2*u5*y*y3)/z + (q*u1^2*u2*u5*v*y*y3)/z + (a1*p*u1*u2*u5*v*y1*y3)/z - (a1*q*u*u1*u2*u5*y1*y3)/z;
k31_3=a2*u2*u3*u5 - a3*u5*y2*y3 - a5*u3*y2*y5 - (p^2*u1^2*u3*u5*w*y2)/z - (q^2*u1^2*u3*u5*w*y2)/z - (d1*p*u*u1^2*u3*u5*y2)/z - (d1*q*u1^2*u3*u5*v*y2)/z + (p*u*u1^2*u3*u5*y*y2)/z + (q*u1^2*u3*u5*v*y*y2)/z + (a1*p*u1*u3*u5*v*y1*y2)/z - (a1*q*u*u1*u3*u5*y1*y2)/z;
k32_1=a4*u2*u4*u5 - a2*u5*y2*y4 - a5*u2*y4*y5 - (p^2*u1^2*u2*u5*w*y4)/z - (q^2*u1^2*u2*u5*w*y4)/z - (d1*p*u*u1^2*u2*u5*y4)/z - (d1*q*u1^2*u2*u5*v*y4)/z + (p*u*u1^2*u2*u5*y*y4)/z + (q*u1^2*u2*u5*v*y*y4)/z + (a1*p*u1*u2*u5*v*y1*y4)/z - (a1*q*u*u1*u2*u5*y1*y4)/z;
k32_2=(a5*d1*u^2*u1^2*u2*y3)/z - d4*u2*u5*y3*y4 - d5*u2*u5*y3 - d3*u2*u5*y4 + (a5*d1*u1^2*u2*v^2*y3)/z - (a5*u^2*u1^2*u2*y*y3)/z - (a5*u1^2*u2*v^2*y*y3)/z + (a5*p*u*u1^2*u2*w*y3)/z + (a5*q*u1^2*u2*v*w*y3)/z;
k32_3=(a5*d1*u^2*u1^2*u3*y2)/z - d4*u3*u5*y2*y4 - d5*u3*u5*y2 + (a5*d1*u1^2*u3*v^2*y2)/z - (a5*u^2*u1^2*u3*y*y2)/z - (a5*u1^2*u3*v^2*y*y2)/z + (a5*p*u*u1^2*u3*w*y2)/z + (a5*q*u1^2*u3*v*w*y2)/z;
k33_1= (p*u*u1^2*u2*u5*y)/z - a5*u2*y5 - (p^2*u1^2*u2*u5*w)/z - (q^2*u1^2*u2*u5*w)/z - (d1*p*u*u1^2*u2*u5)/z - (d1*q*u1^2*u2*u5*v)/z - a2*u5*y2 + (q*u1^2*u2*u5*v*y)/z + (a1*p*u1*u2*u5*v*y1)/z - (a1*q*u*u1*u2*u5*y1)/z;
k33_2=(a5*d1*u^2*u1^2*u2*y3*y4)/z - d4*u2*u5*y3 - d5*u2*u5*y3*y4 - d3*u2*u5 + (a5*d1*u1^2*u2*v^2*y3*y4)/z - (a5*u^2*u1^2*u2*y*y3*y4)/z - (a5*u1^2*u2*v^2*y*y3*y4)/z + (a5*p*u*u1^2*u2*w*y3*y4)/z + (a5*q*u1^2*u2*v*w*y3*y4)/z;
k33_3=(a5*d1*u^2*u1^2*u3*y2*y4)/z - d5*u3*u5*y2*y4 - d4*u3*u5*y2 + (a5*d1*u1^2*u3*v^2*y2*y4)/z - (a5*u^2*u1^2*u3*y*y2*y4)/z - (a5*u1^2*u3*v^2*y*y2*y4)/z + (a5*p*u*u1^2*u3*w*y2*y4)/z + (a5*q*u1^2*u3*v*w*y2*y4)/z;
k34_1=(a5*d1*u^2*u1^2*u2)/z - d5*u2*u5 - d3*u2*u5*y3*y4 - d4*u2*u5*y4 + (a5*d1*u1^2*u2*v^2)/z - (a5*u^2*u1^2*u2*y)/z - (a5*u1^2*u2*v^2*y)/z + (a5*p*u*u1^2*u2*w)/z + (a5*q*u1^2*u2*v*w)/z;
k34_2=a2*u5*y2*y3*y4 - a4*u2*u4*u5*y3 - a3*u2*u3*u5*y4 + a5*u2*y3*y4*y5 + (p^2*u1^2*u2*u5*w*y3*y4)/z + (q^2*u1^2*u2*u5*w*y3*y4)/z + (d1*p*u*u1^2*u2*u5*y3*y4)/z + (d1*q*u1^2*u2*u5*v*y3*y4)/z - (p*u*u1^2*u2*u5*y*y3*y4)/z - (q*u1^2*u2*u5*v*y*y3*y4)/z - (a1*p*u1*u2*u5*v*y1*y3*y4)/z + (a1*q*u*u1*u2*u5*y1*y3*y4)/z;
k34_3=a3*u5*y2*y3*y4 - a4*u3*u4*u5*y2 - a2*u2*u3*u5*y4 + a5*u3*y2*y4*y5 + (p^2*u1^2*u3*u5*w*y2*y4)/z + (q^2*u1^2*u3*u5*w*y2*y4)/z + (d1*p*u*u1^2*u3*u5*y2*y4)/z + (d1*q*u1^2*u3*u5*v*y2*y4)/z - (p*u*u1^2*u3*u5*y*y2*y4)/z - (q*u1^2*u3*u5*v*y*y2*y4)/z - (a1*p*u1*u3*u5*v*y1*y2*y4)/z + (a1*q*u*u1*u3*u5*y1*y2*y4)/z;
k35_1=a5*u2*u4*u5 - a2*u4*y2*y5 - a4*u2*y4*y5 + (d1*d5*u^2*u1^2*u2*u4)/z + (d1*d5*u1^2*u2*u4*v^2)/z - (d5*u^2*u1^2*u2*u4*y)/z - (d5*u1^2*u2*u4*v^2*y)/z - (p^2*u1^2*u2*u4*w*y5)/z - (q^2*u1^2*u2*u4*w*y5)/z + (d5*p*u*u1^2*u2*u4*w)/z + (d5*q*u1^2*u2*u4*v*w)/z - (d1*p*u*u1^2*u2*u4*y5)/z - (d1*q*u1^2*u2*u4*v*y5)/z + (p*u*u1^2*u2*u4*y*y5)/z + (q*u1^2*u2*u4*v*y*y5)/z + (a1*p*u1*u2*u4*v*y1*y5)/z - (a1*q*u*u1*u2*u4*y1*y5)/z;
k35_2= (a4*d1*u^2*u1^2*u2*y3)/z - d4*u2*u4*y3*y5 - d3*u2*u4*y5 + (a4*d1*u1^2*u2*v^2*y3)/z - (a4*u^2*u1^2*u2*y*y3)/z - (a4*u1^2*u2*v^2*y*y3)/z + (a4*p*u*u1^2*u2*w*y3)/z + (a4*q*u1^2*u2*v*w*y3)/z;
k35_3=(a4*d1*u^2*u1^2*u3*y2)/z - d4*u3*u4*y2*y5 + (a4*d1*u1^2*u3*v^2*y2)/z - (a4*u^2*u1^2*u3*y*y2)/z - (a4*u1^2*u3*v^2*y*y2)/z + (a4*p*u*u1^2*u3*w*y2)/z + (a4*q*u1^2*u3*v*w*y2)/z;
k36_1= (a4*d1*u^2*u1^2*u2)/z - d4*u2*u4*y5 - d3*u2*u4*y3*y5 + (a4*d1*u1^2*u2*v^2)/z - (a4*u^2*u1^2*u2*y)/z - (a4*u1^2*u2*v^2*y)/z + (a4*p*u*u1^2*u2*w)/z + (a4*q*u1^2*u2*v*w)/z;
k36_2=a2*u4*y2*y3*y5 - a5*u2*u4*u5*y3 - a3*u2*u3*u4*y5 + a4*u2*y3*y4*y5 - (d1*d5*u^2*u1^2*u2*u4*y3)/z - (d1*d5*u1^2*u2*u4*v^2*y3)/z + (d5*u^2*u1^2*u2*u4*y*y3)/z + (d5*u1^2*u2*u4*v^2*y*y3)/z + (p^2*u1^2*u2*u4*w*y3*y5)/z + (q^2*u1^2*u2*u4*w*y3*y5)/z - (d5*p*u*u1^2*u2*u4*w*y3)/z - (d5*q*u1^2*u2*u4*v*w*y3)/z + (d1*p*u*u1^2*u2*u4*y3*y5)/z + (d1*q*u1^2*u2*u4*v*y3*y5)/z - (p*u*u1^2*u2*u4*y*y3*y5)/z - (q*u1^2*u2*u4*v*y*y3*y5)/z - (a1*p*u1*u2*u4*v*y1*y3*y5)/z + (a1*q*u*u1*u2*u4*y1*y3*y5)/z;
k36_3=a3*u4*y2*y3*y5 - a5*u3*u4*u5*y2 - a2*u2*u3*u4*y5 + a4*u3*y2*y4*y5 - (d1*d5*u^2*u1^2*u3*u4*y2)/z - (d1*d5*u1^2*u3*u4*v^2*y2)/z + (d5*u^2*u1^2*u3*u4*y*y2)/z + (d5*u1^2*u3*u4*v^2*y*y2)/z + (p^2*u1^2*u3*u4*w*y2*y5)/z + (q^2*u1^2*u3*u4*w*y2*y5)/z - (d5*p*u*u1^2*u3*u4*w*y2)/z - (d5*q*u1^2*u3*u4*v*w*y2)/z + (d1*p*u*u1^2*u3*u4*y2*y5)/z + (d1*q*u1^2*u3*u4*v*y2*y5)/z - (p*u*u1^2*u3*u4*y*y2*y5)/z - (q*u1^2*u3*u4*v*y*y2*y5)/z - (a1*p*u1*u3*u4*v*y1*y2*y5)/z + (a1*q*u*u1*u3*u4*y1*y2*y5)/z;
k37_1=0;
k37_2= d5*u2*u3*u4*u5 - (a5*d1*u^2*u1^2*u2*u3*u4)/z - (a5*d1*u1^2*u2*u3*u4*v^2)/z + (a5*u^2*u1^2*u2*u3*u4*y)/z + (a5*u1^2*u2*u3*u4*v^2*y)/z - (a5*p*u*u1^2*u2*u3*u4*w)/z - (a5*q*u1^2*u2*u3*u4*v*w)/z;
k37_3= (a5*d1*u^2*u1^2*u4*y2*y3)/z - d5*u4*u5*y2*y3 + (a5*d1*u1^2*u4*v^2*y2*y3)/z - (a5*u^2*u1^2*u4*y*y2*y3)/z - (a5*u1^2*u4*v^2*y*y2*y3)/z + (a5*p*u*u1^2*u4*w*y2*y3)/z + (a5*q*u1^2*u4*v*w*y2*y3)/z;
k38_1= d3*u2*u3*u4*u5;
k38_2= (p*u*u1^2*u2*u3*u4*u5*y)/z - a3*u2*u4*u5*y3 - a4*u2*u3*u5*y4 - a5*u2*u3*u4*y5 - (p^2*u1^2*u2*u3*u4*u5*w)/z - (q^2*u1^2*u2*u3*u4*u5*w)/z - (d1*p*u*u1^2*u2*u3*u4*u5)/z - (d1*q*u1^2*u2*u3*u4*u5*v)/z - a2*u3*u4*u5*y2 + (q*u1^2*u2*u3*u4*u5*v*y)/z + (a1*p*u1*u2*u3*u4*u5*v*y1)/z - (a1*q*u*u1*u2*u3*u4*u5*y1)/z;
k38_3=a4*u5*y2*y3*y4 - a3*u3*u4*u5*y2 - a2*u2*u4*u5*y3 + a5*u4*y2*y3*y5 + (p^2*u1^2*u4*u5*w*y2*y3)/z + (q^2*u1^2*u4*u5*w*y2*y3)/z + (d1*p*u*u1^2*u4*u5*y2*y3)/z + (d1*q*u1^2*u4*u5*v*y2*y3)/z - (p*u*u1^2*u4*u5*y*y2*y3)/z - (q*u1^2*u4*u5*v*y*y2*y3)/z - (a1*p*u1*u4*u5*v*y1*y2*y3)/z + (a1*q*u*u1*u4*u5*y1*y2*y3)/z;
k39_1= (a3*d1*u^2*u1^2*u2)/z - d3*u2*u3*y4*y5 + (a3*d1*u1^2*u2*v^2)/z - (a3*u^2*u1^2*u2*y)/z - (a3*u1^2*u2*v^2*y)/z + (a3*p*u*u1^2*u2*w)/z + (a3*q*u1^2*u2*v*w)/z;
k39_2=a2*u3*y2*y4*y5 - a5*u2*u3*u5*y4 - a4*u2*u3*u4*y5 + a3*u2*y3*y4*y5 - (d1*d4*u^2*u1^2*u2*u3)/z - (d1*d4*u1^2*u2*u3*v^2)/z + (d4*u^2*u1^2*u2*u3*y)/z + (d4*u1^2*u2*u3*v^2*y)/z - (d1*d5*u^2*u1^2*u2*u3*y4)/z - (d1*d5*u1^2*u2*u3*v^2*y4)/z + (d5*u^2*u1^2*u2*u3*y*y4)/z + (d5*u1^2*u2*u3*v^2*y*y4)/z + (p^2*u1^2*u2*u3*w*y4*y5)/z + (q^2*u1^2*u2*u3*w*y4*y5)/z - (d4*p*u*u1^2*u2*u3*w)/z - (d4*q*u1^2*u2*u3*v*w)/z - (d5*p*u*u1^2*u2*u3*w*y4)/z - (d5*q*u1^2*u2*u3*v*w*y4)/z + (d1*p*u*u1^2*u2*u3*y4*y5)/z + (d1*q*u1^2*u2*u3*v*y4*y5)/z - (p*u*u1^2*u2*u3*y*y4*y5)/z - (q*u1^2*u2*u3*v*y*y4*y5)/z - (a1*p*u1*u2*u3*v*y1*y4*y5)/z + (a1*q*u*u1*u2*u3*y1*y4*y5)/z;
k39_3=p*v*y1 - a1*u1*w - q*u*y1 + (d1^2*u^2*u1^2*y1)/z + (d1^2*u1^2*v^2*y1)/z + (p^2*u1^2*w^2*y1)/z + (q^2*u1^2*w^2*y1)/z + (u^2*u1^2*y^2*y1)/z + (u1^2*v^2*y^2*y1)/z + a2*u2*y3*y4*y5 + a3*u3*y2*y4*y5 + a4*u4*y2*y3*y5 + a5*u5*y2*y3*y4 + (d1*d2*u^2*u1^2)/z + (d1*d2*u1^2*v^2)/z - (d2*u^2*u1^2*y)/z - (d2*u1^2*v^2*y)/z + (d1*d3*u^2*u1^2*y2)/z + (d1*d3*u1^2*v^2*y2)/z - (2*d1*u^2*u1^2*y*y1)/z - (d3*u^2*u1^2*y*y2)/z - (2*d1*u1^2*v^2*y*y1)/z - (d3*u1^2*v^2*y*y2)/z + (d2*p*u*u1^2*w)/z + (d2*q*u1^2*v*w)/z - (a1*p*u1*v*w*y1^2)/z + (a1*q*u*u1*w*y1^2)/z + (2*d1*p*u*u1^2*w*y1)/z + (d3*p*u*u1^2*w*y2)/z + (2*d1*q*u1^2*v*w*y1)/z + (d3*q*u1^2*v*w*y2)/z - (2*p*u*u1^2*w*y*y1)/z - (2*q*u1^2*v*w*y*y1)/z + (d1*d4*u^2*u1^2*y2*y3)/z + (d1*d4*u1^2*v^2*y2*y3)/z - (d4*u^2*u1^2*y*y2*y3)/z - (d4*u1^2*v^2*y*y2*y3)/z + (d1*d5*u^2*u1^2*y2*y3*y4)/z + (d1*d5*u1^2*v^2*y2*y3*y4)/z - (d5*u^2*u1^2*y*y2*y3*y4)/z - (d5*u1^2*v^2*y*y2*y3*y4)/z - (p^2*u1^2*w*y2*y3*y4*y5)/z - (q^2*u1^2*w*y2*y3*y4*y5)/z + (d4*p*u*u1^2*w*y2*y3)/z + (d4*q*u1^2*v*w*y2*y3)/z + (d5*p*u*u1^2*w*y2*y3*y4)/z + (d5*q*u1^2*v*w*y2*y3*y4)/z - (d1*p*u*u1^2*y2*y3*y4*y5)/z - (d1*q*u1^2*v*y2*y3*y4*y5)/z + (p*u*u1^2*y*y2*y3*y4*y5)/z + (q*u1^2*v*y*y2*y3*y4*y5)/z + (a1*p*u1*v*y1*y2*y3*y4*y5)/z - (a1*q*u*u1*y1*y2*y3*y4*y5)/z;
k41_1=2*a2*d2*u5*y3 - 2*a3*d3*u2*u3*u5 - 2*a4*d5*u2*u4*u5 + 2*a2*d3*u5*y2*y3 + 2*a5*d5*u2*y4*y5 + 2*a2*d4*u5*y2 + 2*a5*d4*u2*y5 + 2*a2*d5*u5*y2*y4 + 2*a5*d3*u2*y3*y5 + (2*a5*p*q*u^2*u1^2*u2*y4)/z - (2*a5*p*q*u1^2*u2*v^2*y4)/z - (2*a5*p^2*u*u1^2*u2*v*y4)/z + (2*a5*q^2*u*u1^2*u2*v*y4)/z - (2*a1*a5*d1*u^2*u1*u2*y1*y4)/z - (2*a1*a5*d1*u1*u2*v^2*y1*y4)/z + (2*a5*d1*p*u1^2*u2*v*w*y4)/z - (2*a5*d1*q*u*u1^2*u2*w*y4)/z + (2*a1*a5*u^2*u1*u2*y*y1*y4)/z + (2*a1*a5*u1*u2*v^2*y*y1*y4)/z - (2*a5*p*u1^2*u2*v*w*y*y4)/z + (2*a5*q*u*u1^2*u2*w*y*y4)/z - (2*a1*a5*p*u*u1*u2*w*y1*y4)/z - (2*a1*a5*q*u1*u2*v*w*y1*y4)/z;
k41_2=- (u2*u5*v*a1^2*p*u1^2*y3)/z + (u*u2*u5*a1^2*q*u1^2*y3)/z - (2*u*u2*u5*y1*a1*d1*p*u1*y3)/z - (2*u2*u5*v*y1*a1*d1*q*u1*y3)/z - (2*u2*u5*w*y1*a1*p^2*u1*y3)/z + (2*u*u2*u5*y1*a1*p*u1*y*y3)/z - (2*u2*u5*w*y1*a1*q^2*u1*y3)/z + (2*u2*u5*v*y1*a1*q*u1*y*y3)/z + u2*u5*a2^2*y3 + 2*u5*y2*a2*a3*u3 - 2*y2*y5*a2*a5*y3 + u2*u5*a3^2*y3 + 2*u2*y5*a3*a5*u3 - u2*u5*a4^2*y3 + u2*u5*a5^2*y3 - (u2*u5*v*d1^2*p*u1^2*y3)/z + (u*u2*u5*d1^2*q*u1^2*y3)/z + (2*u2*u5*v*d1*p*u1^2*y*y3)/z - (2*u*u2*u5*d1*q*u1^2*y*y3)/z - u2*u5*d2^2*y3 + u2*u5*d3^2*y3 + 2*u2*u5*d3*d4 + 2*u2*u5*d3*d5*y4 + u2*u5*d4^2*y3 + 2*u2*u5*d4*d5*y3*y4 + u2*u5*d5^2*y3 - (u2*u5*v*p^3*u1^2*y3)/z + (u*u2*u5*p^2*q*u1^2*y3)/z - (u2*u5*v*p*q^2*u1^2*y3)/z - (u2*u5*v*p*u1^2*y^2*y3)/z + (u*u2*u5*q^3*u1^2*y3)/z + (u*u2*u5*q*u1^2*y^2*y3)/z;
k41_3=- (u3*u5*v*y2*a1^2*p*u1^2)/z + (u*u3*u5*y2*a1^2*q*u1^2)/z - (2*u*u3*u5*y1*y2*a1*d1*p*u1)/z - (2*u3*u5*v*y1*y2*a1*d1*q*u1)/z - (2*u3*u5*w*y1*y2*a1*p^2*u1)/z + (2*u*u3*u5*y1*y2*a1*p*u1*y)/z - (2*u3*u5*w*y1*y2*a1*q^2*u1)/z + (2*u3*u5*v*y1*y2*a1*q*u1*y)/z + u3*u5*y2*a2^2 + 2*u2*u5*y3*a2*a3 + 2*u2*u3*y5*a2*a5 + u3*u5*y2*a3^2 - 2*y2*y3*y5*a3*a5 - u3*u5*y2*a4^2 + u3*u5*y2*a5^2 - (u3*u5*v*y2*d1^2*p*u1^2)/z + (u*u3*u5*y2*d1^2*q*u1^2)/z + (2*u3*u5*v*y2*d1*p*u1^2*y)/z - (2*u*u3*u5*y2*d1*q*u1^2*y)/z - u3*u5*y2*d2^2 - 2*u3*u5*d2*d3 - u3*u5*y2*d3^2 + u3*u5*y2*d4^2 + 2*u3*u5*y2*d4*d5*y4 + u3*u5*y2*d5^2 - (u3*u5*v*y2*p^3*u1^2)/z + (u*u3*u5*y2*p^2*q*u1^2)/z - (u3*u5*v*y2*p*q^2*u1^2)/z - (u3*u5*v*y2*p*u1^2*y^2)/z + (u*u3*u5*y2*q^3*u1^2)/z + (u*u3*u5*y2*q*u1^2*y^2)/z;
k42_1= - (u2*u5*v*a1^2*p*u1^2*y4)/z + (u*u2*u5*a1^2*q*u1^2*y4)/z - (2*u*u2*u5*y1*a1*d1*p*u1*y4)/z - (2*u2*u5*v*y1*a1*d1*q*u1*y4)/z - (2*u2*u5*w*y1*a1*p^2*u1*y4)/z + (2*u*u2*u5*y1*a1*p*u1*y*y4)/z - (2*u2*u5*w*y1*a1*q^2*u1*y4)/z + (2*u2*u5*v*y1*a1*q*u1*y*y4)/z + u2*u5*a2^2*y4 + 2*u5*y2*a2*a4*u4 - 2*y2*y5*a2*a5*y4 - u2*u5*a3^2*y4 + u2*u5*a4^2*y4 + 2*u2*y5*a4*a5*u4 + u2*u5*a5^2*y4 - (u2*u5*v*d1^2*p*u1^2*y4)/z + (u*u2*u5*d1^2*q*u1^2*y4)/z + (2*u2*u5*v*d1*p*u1^2*y*y4)/z - (2*u*u2*u5*d1*q*u1^2*y*y4)/z - u2*u5*d2^2*y4 + u2*u5*d3^2*y4 + 2*u2*u5*d3*d4*y3*y4 + 2*u2*u5*d3*d5*y3 + u2*u5*d4^2*y4 + 2*u2*u5*d4*d5 + u2*u5*d5^2*y4 - (u2*u5*v*p^3*u1^2*y4)/z + (u*u2*u5*p^2*q*u1^2*y4)/z - (u2*u5*v*p*q^2*u1^2*y4)/z - (u2*u5*v*p*u1^2*y^2*y4)/z + (u*u2*u5*q^3*u1^2*y4)/z + (u*u2*u5*q*u1^2*y^2*y4)/z;
k42_2=2*a4*d3*u2*u4*u5 - 2*a2*d3*u5*y2*y4 - 2*a5*d5*u2*y3*y5 - 2*a2*d4*u5*y2*y3*y4 - 2*a5*d4*u2*y3*y4*y5 - 2*a2*d2*u5*y4 + 2*a3*d5*u2*u3*u5 - 2*a2*d5*u5*y2*y3 - 2*a5*d3*u2*y4*y5 + 2*a3*d4*u2*u3*u5*y4 + 2*a4*d4*u2*u4*u5*y3 - (2*a5*p*q*u^2*u1^2*u2*y3)/z + (2*a5*p*q*u1^2*u2*v^2*y3)/z + (2*a5*p^2*u*u1^2*u2*v*y3)/z - (2*a5*q^2*u*u1^2*u2*v*y3)/z + (2*a1*a5*d1*u^2*u1*u2*y1*y3)/z + (2*a1*a5*d1*u1*u2*v^2*y1*y3)/z - (2*a5*d1*p*u1^2*u2*v*w*y3)/z + (2*a5*d1*q*u*u1^2*u2*w*y3)/z - (2*a1*a5*u^2*u1*u2*y*y1*y3)/z - (2*a1*a5*u1*u2*v^2*y*y1*y3)/z + (2*a5*p*u1^2*u2*v*w*y*y3)/z - (2*a5*q*u*u1^2*u2*w*y*y3)/z + (2*a1*a5*p*u*u1*u2*w*y1*y3)/z + (2*a1*a5*q*u1*u2*v*w*y1*y3)/z;
k42_3=2*a2*d5*u2*u3*u5 - 2*a3*d2*u5*y4 - 2*a3*d3*u5*y2*y4 - 2*a5*d5*u3*y2*y5 - 2*a3*d4*u5*y2*y3*y4 - 2*a5*d4*u3*y2*y4*y5 - 2*a3*d5*u5*y2*y3 + 2*a2*d4*u2*u3*u5*y4 + 2*a4*d4*u3*u4*u5*y2 - (2*a5*p*q*u^2*u1^2*u3*y2)/z + (2*a5*p*q*u1^2*u3*v^2*y2)/z + (2*a5*p^2*u*u1^2*u3*v*y2)/z - (2*a5*q^2*u*u1^2*u3*v*y2)/z + (2*a1*a5*d1*u^2*u1*u3*y1*y2)/z + (2*a1*a5*d1*u1*u3*v^2*y1*y2)/z - (2*a5*d1*p*u1^2*u3*v*w*y2)/z + (2*a5*d1*q*u*u1^2*u3*w*y2)/z - (2*a1*a5*u^2*u1*u3*y*y1*y2)/z - (2*a1*a5*u1*u3*v^2*y*y1*y2)/z + (2*a5*p*u1^2*u3*v*w*y*y2)/z - (2*a5*q*u*u1^2*u3*w*y*y2)/z + (2*a1*a5*p*u*u1*u3*w*y1*y2)/z + (2*a1*a5*q*u1*u3*v*w*y1*y2)/z;
k43_1=- (u2*u5*v*a1^2*p*u1^2)/z + (u*u2*u5*a1^2*q*u1^2)/z - (2*u*u2*u5*y1*a1*d1*p*u1)/z - (2*u2*u5*v*y1*a1*d1*q*u1)/z - (2*u2*u5*w*y1*a1*p^2*u1)/z + (2*u*u2*u5*y1*a1*p*u1*y)/z - (2*u2*u5*w*y1*a1*q^2*u1)/z + (2*u2*u5*v*y1*a1*q*u1*y)/z + u2*u5*a2^2 - 2*y2*y5*a2*a5 - u2*u5*a3^2 - u2*u5*a4^2 + u2*u5*a5^2 - (u2*u5*v*d1^2*p*u1^2)/z + (u*u2*u5*d1^2*q*u1^2)/z + (2*u2*u5*v*d1*p*u1^2*y)/z - (2*u*u2*u5*d1*q*u1^2*y)/z - u2*u5*d2^2 + u2*u5*d3^2 + 2*u2*u5*d3*d4*y3 + 2*u2*u5*d3*d5*y3*y4 + u2*u5*d4^2 + 2*u2*u5*d4*d5*y4 + u2*u5*d5^2 - (u2*u5*v*p^3*u1^2)/z + (u*u2*u5*p^2*q*u1^2)/z - (u2*u5*v*p*q^2*u1^2)/z - (u2*u5*v*p*u1^2*y^2)/z + (u*u2*u5*q^3*u1^2)/z + (u*u2*u5*q*u1^2*y^2)/z;
k43_2=2*a3*d4*u2*u3*u5 - 2*a2*d3*u5*y2 - 2*a2*d2*u5 - 2*a2*d4*u5*y2*y3 - 2*a2*d5*u5*y2*y3*y4 - 2*a5*d5*u2*y3*y4*y5 - 2*a5*d4*u2*y3*y5 - 2*a5*d3*u2*y5 + 2*a3*d5*u2*u3*u5*y4 + 2*a4*d5*u2*u4*u5*y3 - (2*a5*p*q*u^2*u1^2*u2*y3*y4)/z + (2*a5*p*q*u1^2*u2*v^2*y3*y4)/z + (2*a5*p^2*u*u1^2*u2*v*y3*y4)/z - (2*a5*q^2*u*u1^2*u2*v*y3*y4)/z + (2*a1*a5*d1*u^2*u1*u2*y1*y3*y4)/z + (2*a1*a5*d1*u1*u2*v^2*y1*y3*y4)/z - (2*a5*d1*p*u1^2*u2*v*w*y3*y4)/z + (2*a5*d1*q*u*u1^2*u2*w*y3*y4)/z - (2*a1*a5*u^2*u1*u2*y*y1*y3*y4)/z - (2*a1*a5*u1*u2*v^2*y*y1*y3*y4)/z + (2*a5*p*u1^2*u2*v*w*y*y3*y4)/z - (2*a5*q*u*u1^2*u2*w*y*y3*y4)/z + (2*a1*a5*p*u*u1*u2*w*y1*y3*y4)/z + (2*a1*a5*q*u1*u2*v*w*y1*y3*y4)/z;
k43_3=2*a2*d4*u2*u3*u5 - 2*a3*d3*u5*y2 - 2*a3*d2*u5 - 2*a3*d4*u5*y2*y3 - 2*a3*d5*u5*y2*y3*y4 - 2*a5*d5*u3*y2*y4*y5 - 2*a5*d4*u3*y2*y5 + 2*a2*d5*u2*u3*u5*y4 + 2*a4*d5*u3*u4*u5*y2 - (2*a5*p*q*u^2*u1^2*u3*y2*y4)/z + (2*a5*p*q*u1^2*u3*v^2*y2*y4)/z + (2*a5*p^2*u*u1^2*u3*v*y2*y4)/z - (2*a5*q^2*u*u1^2*u3*v*y2*y4)/z + (2*a1*a5*d1*u^2*u1*u3*y1*y2*y4)/z + (2*a1*a5*d1*u1*u3*v^2*y1*y2*y4)/z - (2*a5*d1*p*u1^2*u3*v*w*y2*y4)/z + (2*a5*d1*q*u*u1^2*u3*w*y2*y4)/z - (2*a1*a5*u^2*u1*u3*y*y1*y2*y4)/z - (2*a1*a5*u1*u3*v^2*y*y1*y2*y4)/z + (2*a5*p*u1^2*u3*v*w*y*y2*y4)/z - (2*a5*q*u*u1^2*u3*w*y*y2*y4)/z + (2*a1*a5*p*u*u1*u3*w*y1*y2*y4)/z + (2*a1*a5*q*u1*u3*v*w*y1*y2*y4)/z;
k44_1= 2*a4*d4*u2*u4*u5 - 2*a5*d5*u2*y5 - 2*a2*d2*u5*y3*y4 - 2*a5*d4*u2*y4*y5 - 2*a2*d3*u5*y2*y3*y4 - 2*a5*d3*u2*y3*y4*y5 - 2*a2*d4*u5*y2*y4 - 2*a2*d5*u5*y2 + 2*a3*d3*u2*u3*u5*y4 + 2*a4*d3*u2*u4*u5*y3 - (2*a5*p*q*u^2*u1^2*u2)/z + (2*a5*p*q*u1^2*u2*v^2)/z + (2*a5*p^2*u*u1^2*u2*v)/z - (2*a5*q^2*u*u1^2*u2*v)/z + (2*a1*a5*d1*u^2*u1*u2*y1)/z + (2*a1*a5*d1*u1*u2*v^2*y1)/z - (2*a5*d1*p*u1^2*u2*v*w)/z + (2*a5*d1*q*u*u1^2*u2*w)/z - (2*a1*a5*u^2*u1*u2*y*y1)/z - (2*a1*a5*u1*u2*v^2*y*y1)/z + (2*a5*p*u1^2*u2*v*w*y)/z - (2*a5*q*u*u1^2*u2*w*y)/z + (2*a1*a5*p*u*u1*u2*w*y1)/z + (2*a1*a5*q*u1*u2*v*w*y1)/z;
k44_2=(u2*u5*v*a1^2*p*u1^2*y3*y4)/z - (u*u2*u5*a1^2*q*u1^2*y3*y4)/z + (2*u*u2*u5*y1*a1*d1*p*u1*y3*y4)/z + (2*u2*u5*v*y1*a1*d1*q*u1*y3*y4)/z + (2*u2*u5*w*y1*a1*p^2*u1*y3*y4)/z - (2*u*u2*u5*y1*a1*p*u1*y*y3*y4)/z + (2*u2*u5*w*y1*a1*q^2*u1*y3*y4)/z - (2*u2*u5*v*y1*a1*q*u1*y*y3*y4)/z - u2*u5*a2^2*y3*y4 - 2*u5*y2*a2*a3*u3*y4 - 2*u5*y2*a2*a4*u4*y3 + 2*y2*y5*a2*a5*y3*y4 - u2*u5*a3^2*y3*y4 + 2*u2*u5*a3*a4*u3*u4 - 2*u2*y5*a3*a5*u3*y4 - u2*u5*a4^2*y3*y4 - 2*u2*y5*a4*a5*u4*y3 - u2*u5*a5^2*y3*y4 + (u2*u5*v*d1^2*p*u1^2*y3*y4)/z - (u*u2*u5*d1^2*q*u1^2*y3*y4)/z - (2*u2*u5*v*d1*p*u1^2*y*y3*y4)/z + (2*u*u2*u5*d1*q*u1^2*y*y3*y4)/z + u2*u5*d2^2*y3*y4 - u2*u5*d3^2*y3*y4 - 2*u2*u5*d3*d4*y4 - 2*u2*u5*d3*d5 - u2*u5*d4^2*y3*y4 - 2*u2*u5*d4*d5*y3 - u2*u5*d5^2*y3*y4 + (u2*u5*v*p^3*u1^2*y3*y4)/z - (u*u2*u5*p^2*q*u1^2*y3*y4)/z + (u2*u5*v*p*q^2*u1^2*y3*y4)/z + (u2*u5*v*p*u1^2*y^2*y3*y4)/z - (u*u2*u5*q^3*u1^2*y3*y4)/z - (u*u2*u5*q*u1^2*y^2*y3*y4)/z;
k44_3= (u3*u5*v*y2*a1^2*p*u1^2*y4)/z - (u*u3*u5*y2*a1^2*q*u1^2*y4)/z + (2*u*u3*u5*y1*y2*a1*d1*p*u1*y4)/z + (2*u3*u5*v*y1*y2*a1*d1*q*u1*y4)/z + (2*u3*u5*w*y1*y2*a1*p^2*u1*y4)/z - (2*u*u3*u5*y1*y2*a1*p*u1*y*y4)/z + (2*u3*u5*w*y1*y2*a1*q^2*u1*y4)/z - (2*u3*u5*v*y1*y2*a1*q*u1*y*y4)/z - u3*u5*y2*a2^2*y4 - 2*u2*u5*y3*a2*a3*y4 + 2*u2*u3*u5*a2*a4*u4 - 2*u2*u3*y5*a2*a5*y4 - u3*u5*y2*a3^2*y4 - 2*u5*y2*y3*a3*a4*u4 + 2*y2*y3*y5*a3*a5*y4 - u3*u5*y2*a4^2*y4 - 2*u3*y2*y5*a4*a5*u4 - u3*u5*y2*a5^2*y4 + (u3*u5*v*y2*d1^2*p*u1^2*y4)/z - (u*u3*u5*y2*d1^2*q*u1^2*y4)/z - (2*u3*u5*v*y2*d1*p*u1^2*y*y4)/z + (2*u*u3*u5*y2*d1*q*u1^2*y*y4)/z + u3*u5*y2*d2^2*y4 + 2*u3*u5*d2*d3*y4 + u3*u5*y2*d3^2*y4 - u3*u5*y2*d4^2*y4 - 2*u3*u5*y2*d4*d5 - u3*u5*y2*d5^2*y4 + (u3*u5*v*y2*p^3*u1^2*y4)/z - (u*u3*u5*y2*p^2*q*u1^2*y4)/z + (u3*u5*v*y2*p*q^2*u1^2*y4)/z + (u3*u5*v*y2*p*u1^2*y^2*y4)/z - (u*u3*u5*y2*q^3*u1^2*y4)/z - (u*u3*u5*y2*q*u1^2*y^2*y4)/z;
k45_1= - (u2*u4*y5*a1^2*p*u1^2*v)/z + (u2*u4*y5*a1^2*q*u*u1^2)/z + (2*u2*u4*y1*a1*d1*d5*u^2*u1)/z + (2*u2*u4*y1*a1*d1*d5*u1*v^2)/z - (2*u2*u4*y1*y5*a1*d1*p*u*u1)/z - (2*u2*u4*y1*y5*a1*d1*q*u1*v)/z + (2*u2*u4*w*y1*a1*d5*p*u*u1)/z + (2*u2*u4*w*y1*a1*d5*q*u1*v)/z - (2*u2*u4*y1*a1*d5*u^2*u1*y)/z - (2*u2*u4*y1*a1*d5*u1*v^2*y)/z - (2*u2*u4*w*y1*y5*a1*p^2*u1)/z + (2*u2*u4*y1*y5*a1*p*u*u1*y)/z - (2*u2*u4*w*y1*y5*a1*q^2*u1)/z + (2*u2*u4*y1*y5*a1*q*u1*v*y)/z + u2*u4*y5*a2^2 - 2*y2*y4*y5*a2*a4 + 2*u4*u5*y2*a2*a5 - u2*u4*y5*a3^2 + u2*u4*y5*a4^2 + 2*u2*u5*y4*a4*a5 + u2*u4*y5*a5^2 - (u2*u4*y5*d1^2*p*u1^2*v)/z + (u2*u4*y5*d1^2*q*u*u1^2)/z - (2*u2*u4*w*d1*d5*p*u1^2*v)/z + (2*u2*u4*w*d1*d5*q*u*u1^2)/z + (2*u2*u4*y5*d1*p*u1^2*v*y)/z - (2*u2*u4*y5*d1*q*u*u1^2*y)/z - u2*u4*y5*d2^2 + u2*u4*y5*d3^2 + 2*u2*u4*y5*d3*d4*y3 + u2*u4*y5*d4^2 - u2*u4*y5*d5^2 + (2*u2*u4*d5*p^2*u*u1^2*v)/z - (2*u2*u4*d5*p*q*u^2*u1^2)/z + (2*u2*u4*d5*p*q*u1^2*v^2)/z + (2*u2*u4*w*d5*p*u1^2*v*y)/z - (2*u2*u4*d5*q^2*u*u1^2*v)/z - (2*u2*u4*w*d5*q*u*u1^2*y)/z - (u2*u4*y5*p^3*u1^2*v)/z + (u2*u4*y5*p^2*q*u*u1^2)/z - (u2*u4*y5*p*q^2*u1^2*v)/z - (u2*u4*y5*p*u1^2*v*y^2)/z + (u2*u4*y5*q^3*u*u1^2)/z + (u2*u4*y5*q*u*u1^2*y^2)/z;
k45_2=2*a5*d3*u2*u4*u5 - 2*a2*d2*u4*y5 - 2*a2*d3*u4*y2*y5 - 2*a4*d5*u2*y3*y5 - 2*a2*d4*u4*y2*y3*y5 - 2*a4*d4*u2*y3*y4*y5 - 2*a4*d3*u2*y4*y5 + 2*a3*d4*u2*u3*u4*y5 + 2*a5*d4*u2*u4*u5*y3 - (2*a4*p*q*u^2*u1^2*u2*y3)/z + (2*a4*p*q*u1^2*u2*v^2*y3)/z + (2*a4*p^2*u*u1^2*u2*v*y3)/z - (2*a4*q^2*u*u1^2*u2*v*y3)/z + (2*a1*a4*d1*u^2*u1*u2*y1*y3)/z + (2*a1*a4*d1*u1*u2*v^2*y1*y3)/z - (2*a4*d1*p*u1^2*u2*v*w*y3)/z + (2*a4*d1*q*u*u1^2*u2*w*y3)/z - (2*a1*a4*u^2*u1*u2*y*y1*y3)/z - (2*a1*a4*u1*u2*v^2*y*y1*y3)/z + (2*a4*p*u1^2*u2*v*w*y*y3)/z - (2*a4*q*u*u1^2*u2*w*y*y3)/z + (2*a1*a4*p*u*u1*u2*w*y1*y3)/z + (2*a1*a4*q*u1*u2*v*w*y1*y3)/z;
k45_3= 2*a2*d4*u2*u3*u4*y5 - 2*a3*d3*u4*y2*y5 - 2*a4*d5*u3*y2*y5 - 2*a3*d4*u4*y2*y3*y5 - 2*a4*d4*u3*y2*y4*y5 - 2*a3*d2*u4*y5 + 2*a5*d4*u3*u4*u5*y2 - (2*a4*p*q*u^2*u1^2*u3*y2)/z + (2*a4*p*q*u1^2*u3*v^2*y2)/z + (2*a4*p^2*u*u1^2*u3*v*y2)/z - (2*a4*q^2*u*u1^2*u3*v*y2)/z + (2*a1*a4*d1*u^2*u1*u3*y1*y2)/z + (2*a1*a4*d1*u1*u3*v^2*y1*y2)/z - (2*a4*d1*p*u1^2*u3*v*w*y2)/z + (2*a4*d1*q*u*u1^2*u3*w*y2)/z - (2*a1*a4*u^2*u1*u3*y*y1*y2)/z - (2*a1*a4*u1*u3*v^2*y*y1*y2)/z + (2*a4*p*u1^2*u3*v*w*y*y2)/z - (2*a4*q*u*u1^2*u3*w*y*y2)/z + (2*a1*a4*p*u*u1*u3*w*y1*y2)/z + (2*a1*a4*q*u1*u3*v*w*y1*y2)/z;
k46_1=2*a5*d4*u2*u4*u5 - 2*a4*d5*u2*y5 - 2*a2*d2*u4*y3*y5 - 2*a4*d4*u2*y4*y5 - 2*a2*d3*u4*y2*y3*y5 - 2*a4*d3*u2*y3*y4*y5 - 2*a2*d4*u4*y2*y5 + 2*a3*d3*u2*u3*u4*y5 + 2*a5*d3*u2*u4*u5*y3 - (2*a4*p*q*u^2*u1^2*u2)/z + (2*a4*p*q*u1^2*u2*v^2)/z + (2*a4*p^2*u*u1^2*u2*v)/z - (2*a4*q^2*u*u1^2*u2*v)/z + (2*a1*a4*d1*u^2*u1*u2*y1)/z + (2*a1*a4*d1*u1*u2*v^2*y1)/z - (2*a4*d1*p*u1^2*u2*v*w)/z + (2*a4*d1*q*u*u1^2*u2*w)/z - (2*a1*a4*u^2*u1*u2*y*y1)/z - (2*a1*a4*u1*u2*v^2*y*y1)/z + (2*a4*p*u1^2*u2*v*w*y)/z - (2*a4*q*u*u1^2*u2*w*y)/z + (2*a1*a4*p*u*u1*u2*w*y1)/z + (2*a1*a4*q*u1*u2*v*w*y1)/z;
k46_2= (u2*u4*y5*a1^2*p*u1^2*v*y3)/z - (u2*u4*y5*a1^2*q*u*u1^2*y3)/z - (2*u2*u4*y1*a1*d1*d5*u^2*u1*y3)/z - (2*u2*u4*y1*a1*d1*d5*u1*v^2*y3)/z + (2*u2*u4*y1*y5*a1*d1*p*u*u1*y3)/z + (2*u2*u4*y1*y5*a1*d1*q*u1*v*y3)/z - (2*u2*u4*w*y1*a1*d5*p*u*u1*y3)/z - (2*u2*u4*w*y1*a1*d5*q*u1*v*y3)/z + (2*u2*u4*y1*a1*d5*u^2*u1*y*y3)/z + (2*u2*u4*y1*a1*d5*u1*v^2*y*y3)/z + (2*u2*u4*w*y1*y5*a1*p^2*u1*y3)/z - (2*u2*u4*y1*y5*a1*p*u*u1*y*y3)/z + (2*u2*u4*w*y1*y5*a1*q^2*u1*y3)/z - (2*u2*u4*y1*y5*a1*q*u1*v*y*y3)/z - u2*u4*y5*a2^2*y3 - 2*u4*y2*y5*a2*a3*u3 + 2*y2*y4*y5*a2*a4*y3 - 2*u4*u5*y2*a2*a5*y3 - u2*u4*y5*a3^2*y3 - 2*u2*y4*y5*a3*a4*u3 + 2*u2*u4*u5*a3*a5*u3 - u2*u4*y5*a4^2*y3 - 2*u2*u5*y4*a4*a5*y3 - u2*u4*y5*a5^2*y3 + (u2*u4*y5*d1^2*p*u1^2*v*y3)/z - (u2*u4*y5*d1^2*q*u*u1^2*y3)/z + (2*u2*u4*w*d1*d5*p*u1^2*v*y3)/z - (2*u2*u4*w*d1*d5*q*u*u1^2*y3)/z - (2*u2*u4*y5*d1*p*u1^2*v*y*y3)/z + (2*u2*u4*y5*d1*q*u*u1^2*y*y3)/z + u2*u4*y5*d2^2*y3 - u2*u4*y5*d3^2*y3 - 2*u2*u4*y5*d3*d4 - u2*u4*y5*d4^2*y3 + u2*u4*y5*d5^2*y3 - (2*u2*u4*d5*p^2*u*u1^2*v*y3)/z + (2*u2*u4*d5*p*q*u^2*u1^2*y3)/z - (2*u2*u4*d5*p*q*u1^2*v^2*y3)/z - (2*u2*u4*w*d5*p*u1^2*v*y*y3)/z + (2*u2*u4*d5*q^2*u*u1^2*v*y3)/z + (2*u2*u4*w*d5*q*u*u1^2*y*y3)/z + (u2*u4*y5*p^3*u1^2*v*y3)/z - (u2*u4*y5*p^2*q*u*u1^2*y3)/z + (u2*u4*y5*p*q^2*u1^2*v*y3)/z + (u2*u4*y5*p*u1^2*v*y^2*y3)/z - (u2*u4*y5*q^3*u*u1^2*y3)/z - (u2*u4*y5*q*u*u1^2*y^2*y3)/z;
k46_3= (u3*u4*y2*y5*a1^2*p*u1^2*v)/z - (u3*u4*y2*y5*a1^2*q*u*u1^2)/z - (2*u3*u4*y1*y2*a1*d1*d5*u^2*u1)/z - (2*u3*u4*y1*y2*a1*d1*d5*u1*v^2)/z + (2*u3*u4*y1*y2*y5*a1*d1*p*u*u1)/z + (2*u3*u4*y1*y2*y5*a1*d1*q*u1*v)/z - (2*u3*u4*w*y1*y2*a1*d5*p*u*u1)/z - (2*u3*u4*w*y1*y2*a1*d5*q*u1*v)/z + (2*u3*u4*y1*y2*a1*d5*u^2*u1*y)/z + (2*u3*u4*y1*y2*a1*d5*u1*v^2*y)/z + (2*u3*u4*w*y1*y2*y5*a1*p^2*u1)/z - (2*u3*u4*y1*y2*y5*a1*p*u*u1*y)/z + (2*u3*u4*w*y1*y2*y5*a1*q^2*u1)/z - (2*u3*u4*y1*y2*y5*a1*q*u1*v*y)/z - u3*u4*y2*y5*a2^2 - 2*u2*u4*y3*y5*a2*a3 - 2*u2*u3*y4*y5*a2*a4 + 2*u2*u3*u4*u5*a2*a5 - u3*u4*y2*y5*a3^2 + 2*y2*y3*y4*y5*a3*a4 - 2*u4*u5*y2*y3*a3*a5 - u3*u4*y2*y5*a4^2 - 2*u3*u5*y2*y4*a4*a5 - u3*u4*y2*y5*a5^2 + (u3*u4*y2*y5*d1^2*p*u1^2*v)/z - (u3*u4*y2*y5*d1^2*q*u*u1^2)/z + (2*u3*u4*w*y2*d1*d5*p*u1^2*v)/z - (2*u3*u4*w*y2*d1*d5*q*u*u1^2)/z - (2*u3*u4*y2*y5*d1*p*u1^2*v*y)/z + (2*u3*u4*y2*y5*d1*q*u*u1^2*y)/z + u3*u4*y2*y5*d2^2 + 2*u3*u4*y5*d2*d3 + u3*u4*y2*y5*d3^2 - u3*u4*y2*y5*d4^2 + u3*u4*y2*y5*d5^2 - (2*u3*u4*y2*d5*p^2*u*u1^2*v)/z + (2*u3*u4*y2*d5*p*q*u^2*u1^2)/z - (2*u3*u4*y2*d5*p*q*u1^2*v^2)/z - (2*u3*u4*w*y2*d5*p*u1^2*v*y)/z + (2*u3*u4*y2*d5*q^2*u*u1^2*v)/z + (2*u3*u4*w*y2*d5*q*u*u1^2*y)/z + (u3*u4*y2*y5*p^3*u1^2*v)/z - (u3*u4*y2*y5*p^2*q*u*u1^2)/z + (u3*u4*y2*y5*p*q^2*u1^2*v)/z + (u3*u4*y2*y5*p*u1^2*v*y^2)/z - (u3*u4*y2*y5*q^3*u*u1^2)/z - (u3*u4*y2*y5*q*u*u1^2*y^2)/z;
k47_1=- 2*a3*a4*u2*u5 - 2*d3*d5*u2*u3*u4*u5;
k47_2=2*a4*d4*u2*u3*u5 + 2*a2*d5*u3*u4*u5*y2 + 2*a3*d5*u2*u4*u5*y3 + 2*a4*d5*u2*u3*u5*y4 + 2*a5*d5*u2*u3*u4*y5 + (2*a5*p*q*u^2*u1^2*u2*u3*u4)/z - (2*a5*p*q*u1^2*u2*u3*u4*v^2)/z - (2*a5*p^2*u*u1^2*u2*u3*u4*v)/z + (2*a5*q^2*u*u1^2*u2*u3*u4*v)/z - (2*a1*a5*d1*u^2*u1*u2*u3*u4*y1)/z - (2*a1*a5*d1*u1*u2*u3*u4*v^2*y1)/z + (2*a5*d1*p*u1^2*u2*u3*u4*v*w)/z - (2*a5*d1*q*u*u1^2*u2*u3*u4*w)/z + (2*a1*a5*u^2*u1*u2*u3*u4*y*y1)/z + (2*a1*a5*u1*u2*u3*u4*v^2*y*y1)/z - (2*a5*p*u1^2*u2*u3*u4*v*w*y)/z + (2*a5*q*u*u1^2*u2*u3*u4*w*y)/z - (2*a1*a5*p*u*u1*u2*u3*u4*w*y1)/z - (2*a1*a5*q*u1*u2*u3*u4*v*w*y1)/z;
k47_3=2*a2*d5*u2*u4*u5*y3 - 2*a4*d3*u5*y2 - 2*a4*d4*u5*y2*y3 - 2*a4*d5*u5*y2*y3*y4 - 2*a5*d5*u4*y2*y3*y5 - 2*a4*d2*u5 + 2*a3*d5*u3*u4*u5*y2 - (2*a5*p*q*u^2*u1^2*u4*y2*y3)/z + (2*a5*p*q*u1^2*u4*v^2*y2*y3)/z + (2*a5*p^2*u*u1^2*u4*v*y2*y3)/z - (2*a5*q^2*u*u1^2*u4*v*y2*y3)/z + (2*a1*a5*d1*u^2*u1*u4*y1*y2*y3)/z + (2*a1*a5*d1*u1*u4*v^2*y1*y2*y3)/z - (2*a5*d1*p*u1^2*u4*v*w*y2*y3)/z + (2*a5*d1*q*u*u1^2*u4*w*y2*y3)/z - (2*a1*a5*u^2*u1*u4*y*y1*y2*y3)/z - (2*a1*a5*u1*u4*v^2*y*y1*y2*y3)/z + (2*a5*p*u1^2*u4*v*w*y*y2*y3)/z - (2*a5*q*u*u1^2*u4*w*y*y2*y3)/z + (2*a1*a5*p*u*u1*u4*w*y1*y2*y3)/z + (2*a1*a5*q*u1*u4*v*w*y1*y2*y3)/z;
k48_1=2*a2*d2*u3*u4*u5 + 2*a3*d4*u2*u4*u5 + 2*a2*d3*u3*u4*u5*y2 + 2*a3*d3*u2*u4*u5*y3 + 2*a4*d3*u2*u3*u5*y4 + 2*a5*d3*u2*u3*u4*y5;
k48_2= - (u2*u3*u4*u5*v*a1^2*p*u1^2)/z + (u*u2*u3*u4*u5*a1^2*q*u1^2)/z - (2*u*u2*u3*u4*u5*y1*a1*d1*p*u1)/z - (2*u2*u3*u4*u5*v*y1*a1*d1*q*u1)/z - (2*u2*u3*u4*u5*w*y1*a1*p^2*u1)/z + (2*u*u2*u3*u4*u5*y1*a1*p*u1*y)/z - (2*u2*u3*u4*u5*w*y1*a1*q^2*u1)/z + (2*u2*u3*u4*u5*v*y1*a1*q*u1*y)/z + u2*u3*u4*u5*a2^2 - 2*u4*u5*y2*y3*a2*a3 - 2*u3*u5*y2*y4*a2*a4 - 2*u3*u4*y2*y5*a2*a5 + u2*u3*u4*u5*a3^2 - 2*u2*u5*y3*y4*a3*a4 - 2*u2*u4*y3*y5*a3*a5 + u2*u3*u4*u5*a4^2 - 2*u2*u3*y4*y5*a4*a5 + u2*u3*u4*u5*a5^2 - (u2*u3*u4*u5*v*d1^2*p*u1^2)/z + (u*u2*u3*u4*u5*d1^2*q*u1^2)/z + (2*u2*u3*u4*u5*v*d1*p*u1^2*y)/z - (2*u*u2*u3*u4*u5*d1*q*u1^2*y)/z - u2*u3*u4*u5*d2^2 + u2*u3*u4*u5*d3^2 - u2*u3*u4*u5*d4^2 + u2*u3*u4*u5*d5^2 - (u2*u3*u4*u5*v*p^3*u1^2)/z + (u*u2*u3*u4*u5*p^2*q*u1^2)/z - (u2*u3*u4*u5*v*p*q^2*u1^2)/z - (u2*u3*u4*u5*v*p*u1^2*y^2)/z + (u*u2*u3*u4*u5*q^3*u1^2)/z + (u*u2*u3*u4*u5*q*u1^2*y^2)/z;
k48_3=(u4*u5*v*y2*y3*a1^2*p*u1^2)/z - (u*u4*u5*y2*y3*a1^2*q*u1^2)/z + (2*u*u4*u5*y1*y2*y3*a1*d1*p*u1)/z + (2*u4*u5*v*y1*y2*y3*a1*d1*q*u1)/z + (2*u4*u5*w*y1*y2*y3*a1*p^2*u1)/z - (2*u*u4*u5*y1*y2*y3*a1*p*u1*y)/z + (2*u4*u5*w*y1*y2*y3*a1*q^2*u1)/z - (2*u4*u5*v*y1*y2*y3*a1*q*u1*y)/z - u4*u5*y2*y3*a2^2 + 2*u2*u3*u4*u5*a2*a3 - 2*u2*u5*y3*y4*a2*a4 - 2*u2*u4*y3*y5*a2*a5 - u4*u5*y2*y3*a3^2 - 2*u3*u5*y2*y4*a3*a4 - 2*u3*u4*y2*y5*a3*a5 - u4*u5*y2*y3*a4^2 + 2*y2*y3*y4*y5*a4*a5 - u4*u5*y2*y3*a5^2 + (u4*u5*v*y2*y3*d1^2*p*u1^2)/z - (u*u4*u5*y2*y3*d1^2*q*u1^2)/z - (2*u4*u5*v*y2*y3*d1*p*u1^2*y)/z + (2*u*u4*u5*y2*y3*d1*q*u1^2*y)/z + u4*u5*y2*y3*d2^2 + 2*u4*u5*y3*d2*d3 + 2*u4*u5*d2*d4 + u4*u5*y2*y3*d3^2 + 2*u4*u5*y2*d3*d4 + u4*u5*y2*y3*d4^2 - u4*u5*y2*y3*d5^2 + (u4*u5*v*y2*y3*p^3*u1^2)/z - (u*u4*u5*y2*y3*p^2*q*u1^2)/z + (u4*u5*v*y2*y3*p*q^2*u1^2)/z + (u4*u5*v*y2*y3*p*u1^2*y^2)/z - (u*u4*u5*y2*y3*q^3*u1^2)/z - (u*u4*u5*y2*y3*q*u1^2*y^2)/z;
k49_1=2*a4*d3*u2*u3*u4*y5 - 2*a2*d2*u3*y4*y5 - 2*a3*d4*u2*y4*y5 - 2*a2*d3*u3*y2*y4*y5 - 2*a3*d3*u2*y3*y4*y5 - 2*a3*d5*u2*y5 + 2*a5*d3*u2*u3*u5*y4 - (2*a3*p*q*u^2*u1^2*u2)/z + (2*a3*p*q*u1^2*u2*v^2)/z + (2*a3*p^2*u*u1^2*u2*v)/z - (2*a3*q^2*u*u1^2*u2*v)/z + (2*a1*a3*d1*u^2*u1*u2*y1)/z + (2*a1*a3*d1*u1*u2*v^2*y1)/z - (2*a3*d1*p*u1^2*u2*v*w)/z + (2*a3*d1*q*u*u1^2*u2*w)/z - (2*a1*a3*u^2*u1*u2*y*y1)/z - (2*a1*a3*u1*u2*v^2*y*y1)/z + (2*a3*p*u1^2*u2*v*w*y)/z - (2*a3*q*u*u1^2*u2*w*y)/z + (2*a1*a3*p*u*u1*u2*w*y1)/z + (2*a1*a3*q*u1*u2*v*w*y1)/z;
k49_2=(u2*u3*y4*y5*a1^2*p*u1^2*v)/z - (u2*u3*y4*y5*a1^2*q*u*u1^2)/z - (2*u2*u3*y1*a1*d1*d4*u^2*u1)/z - (2*u2*u3*y1*a1*d1*d4*u1*v^2)/z - (2*u2*u3*y1*y4*a1*d1*d5*u^2*u1)/z - (2*u2*u3*y1*y4*a1*d1*d5*u1*v^2)/z + (2*u2*u3*y1*y4*y5*a1*d1*p*u*u1)/z + (2*u2*u3*y1*y4*y5*a1*d1*q*u1*v)/z - (2*u2*u3*w*y1*a1*d4*p*u*u1)/z - (2*u2*u3*w*y1*a1*d4*q*u1*v)/z + (2*u2*u3*y1*a1*d4*u^2*u1*y)/z + (2*u2*u3*y1*a1*d4*u1*v^2*y)/z - (2*u2*u3*w*y1*y4*a1*d5*p*u*u1)/z - (2*u2*u3*w*y1*y4*a1*d5*q*u1*v)/z + (2*u2*u3*y1*y4*a1*d5*u^2*u1*y)/z + (2*u2*u3*y1*y4*a1*d5*u1*v^2*y)/z + (2*u2*u3*w*y1*y4*y5*a1*p^2*u1)/z - (2*u2*u3*y1*y4*y5*a1*p*u*u1*y)/z + (2*u2*u3*w*y1*y4*y5*a1*q^2*u1)/z - (2*u2*u3*y1*y4*y5*a1*q*u1*v*y)/z - u2*u3*y4*y5*a2^2 + 2*y2*y3*y4*y5*a2*a3 - 2*u3*u4*y2*y5*a2*a4 - 2*u3*u5*y2*y4*a2*a5 - u2*u3*y4*y5*a3^2 - 2*u2*u4*y3*y5*a3*a4 - 2*u2*u5*y3*y4*a3*a5 - u2*u3*y4*y5*a4^2 + 2*u2*u3*u4*u5*a4*a5 - u2*u3*y4*y5*a5^2 + (u2*u3*y4*y5*d1^2*p*u1^2*v)/z - (u2*u3*y4*y5*d1^2*q*u*u1^2)/z + (2*u2*u3*w*d1*d4*p*u1^2*v)/z - (2*u2*u3*w*d1*d4*q*u*u1^2)/z + (2*u2*u3*w*y4*d1*d5*p*u1^2*v)/z - (2*u2*u3*w*y4*d1*d5*q*u*u1^2)/z - (2*u2*u3*y4*y5*d1*p*u1^2*v*y)/z + (2*u2*u3*y4*y5*d1*q*u*u1^2*y)/z + u2*u3*y4*y5*d2^2 - u2*u3*y4*y5*d3^2 + u2*u3*y4*y5*d4^2 + 2*u2*u3*y5*d4*d5 - (2*u2*u3*d4*p^2*u*u1^2*v)/z + (2*u2*u3*d4*p*q*u^2*u1^2)/z - (2*u2*u3*d4*p*q*u1^2*v^2)/z - (2*u2*u3*w*d4*p*u1^2*v*y)/z + (2*u2*u3*d4*q^2*u*u1^2*v)/z + (2*u2*u3*w*d4*q*u*u1^2*y)/z + u2*u3*y4*y5*d5^2 - (2*u2*u3*y4*d5*p^2*u*u1^2*v)/z + (2*u2*u3*y4*d5*p*q*u^2*u1^2)/z - (2*u2*u3*y4*d5*p*q*u1^2*v^2)/z - (2*u2*u3*w*y4*d5*p*u1^2*v*y)/z + (2*u2*u3*y4*d5*q^2*u*u1^2*v)/z + (2*u2*u3*w*y4*d5*q*u*u1^2*y)/z + (u2*u3*y4*y5*p^3*u1^2*v)/z - (u2*u3*y4*y5*p^2*q*u*u1^2)/z + (u2*u3*y4*y5*p*q^2*u1^2*v)/z + (u2*u3*y4*y5*p*u1^2*v*y^2)/z - (u2*u3*y4*y5*q^3*u*u1^2)/z - (u2*u3*y4*y5*q*u*u1^2*y^2)/z;
k49_3=(a1^2*p*u1^2*v*w*y1)/z - (y2*y3*y4*y5*a1^2*p*u1^2*v)/z - (a1^2*q*u*u1^2*w*y1)/z + (y2*y3*y4*y5*a1^2*q*u*u1^2)/z - a1^2*w*y1 + (2*a1*d1^2*u^2*u1*y1^2)/z + (2*a1*d1^2*u1*v^2*y1^2)/z + (2*a1*d1*d2*u^2*u1*y1)/z + (2*a1*d1*d2*u1*v^2*y1)/z + (2*y2*a1*d1*d3*u^2*u1*y1)/z + (2*y2*a1*d1*d3*u1*v^2*y1)/z + (2*y2*y3*a1*d1*d4*u^2*u1*y1)/z + (2*y2*y3*a1*d1*d4*u1*v^2*y1)/z + (2*y2*y3*y4*a1*d1*d5*u^2*u1*y1)/z + (2*y2*y3*y4*a1*d1*d5*u1*v^2*y1)/z + (4*a1*d1*p*u*u1*w*y1^2)/z - (2*y2*y3*y4*y5*a1*d1*p*u*u1*y1)/z + (4*a1*d1*q*u1*v*w*y1^2)/z - (2*y2*y3*y4*y5*a1*d1*q*u1*v*y1)/z - (4*a1*d1*u^2*u1*y*y1^2)/z - (4*a1*d1*u1*v^2*y*y1^2)/z + (2*a1*d2*p*u*u1*w*y1)/z + (2*a1*d2*q*u1*v*w*y1)/z - (2*a1*d2*u^2*u1*y*y1)/z - (2*a1*d2*u1*v^2*y*y1)/z + (2*y2*a1*d3*p*u*u1*w*y1)/z + (2*y2*a1*d3*q*u1*v*w*y1)/z - (2*y2*a1*d3*u^2*u1*y*y1)/z - (2*y2*a1*d3*u1*v^2*y*y1)/z + (2*y2*y3*a1*d4*p*u*u1*w*y1)/z + (2*y2*y3*a1*d4*q*u1*v*w*y1)/z - (2*y2*y3*a1*d4*u^2*u1*y*y1)/z - (2*y2*y3*a1*d4*u1*v^2*y*y1)/z + (2*y2*y3*y4*a1*d5*p*u*u1*w*y1)/z + (2*y2*y3*y4*a1*d5*q*u1*v*w*y1)/z - (2*y2*y3*y4*a1*d5*u^2*u1*y*y1)/z - (2*y2*y3*y4*a1*d5*u1*v^2*y*y1)/z + (2*a1*p^2*u1*w^2*y1^2)/z - (2*y2*y3*y4*y5*a1*p^2*u1*w*y1)/z - (4*a1*p*u*u1*w*y*y1^2)/z + (2*y2*y3*y4*y5*a1*p*u*u1*y*y1)/z - 2*a1*p*u1*v + (2*a1*q^2*u1*w^2*y1^2)/z - (2*y2*y3*y4*y5*a1*q^2*u1*w*y1)/z + 2*a1*q*u*u1 - (4*a1*q*u1*v*w*y*y1^2)/z + (2*y2*y3*y4*y5*a1*q*u1*v*y*y1)/z + (2*a1*u^2*u1*y^2*y1^2)/z + (2*a1*u1*v^2*y^2*y1^2)/z + y2*y3*y4*y5*a2^2 - 2*u2*u3*y4*y5*a2*a3 - 2*u2*u4*y3*y5*a2*a4 - 2*u2*u5*y3*y4*a2*a5 + y2*y3*y4*y5*a3^2 - 2*u3*u4*y2*y5*a3*a4 - 2*u3*u5*y2*y4*a3*a5 + y2*y3*y4*y5*a4^2 - 2*u4*u5*y2*y3*a4*a5 + y2*y3*y4*y5*a5^2 - (d1^2*p*u1^2*v*w*y1)/z - (y2*y3*y4*y5*d1^2*p*u1^2*v)/z + (d1^2*q*u*u1^2*w*y1)/z + (y2*y3*y4*y5*d1^2*q*u*u1^2)/z + d1^2*w*y1 - (2*d1*d2*p*u1^2*v*w)/z + (2*d1*d2*q*u*u1^2*w)/z - (2*y2*d1*d3*p*u1^2*v*w)/z + (2*y2*d1*d3*q*u*u1^2*w)/z - (2*y2*y3*d1*d4*p*u1^2*v*w)/z + (2*y2*y3*d1*d4*q*u*u1^2*w)/z - (2*y2*y3*y4*d1*d5*p*u1^2*v*w)/z + (2*y2*y3*y4*d1*d5*q*u*u1^2*w)/z + (2*d1*p^2*u*u1^2*v*y1)/z - (2*d1*p*q*u^2*u1^2*y1)/z + (2*d1*p*q*u1^2*v^2*y1)/z - 2*d1*p*u*y1 + (2*d1*p*u1^2*v*w*y*y1)/z + (2*y2*y3*y4*y5*d1*p*u1^2*v*y)/z - (2*d1*q^2*u*u1^2*v*y1)/z - (2*d1*q*u*u1^2*w*y*y1)/z - (2*y2*y3*y4*y5*d1*q*u*u1^2*y)/z - 2*d1*q*v*y1 - 2*d1*w*y*y1 - y2*y3*y4*y5*d2^2 - 2*y3*y4*y5*d2*d3 - 2*y4*y5*d2*d4 - 2*y5*d2*d5 + (2*d2*p^2*u*u1^2*v)/z - (2*d2*p*q*u^2*u1^2)/z + (2*d2*p*q*u1^2*v^2)/z + (2*d2*p*u1^2*v*w*y)/z - (2*d2*q^2*u*u1^2*v)/z - (2*d2*q*u*u1^2*w*y)/z - y2*y3*y4*y5*d3^2 - 2*y2*y4*y5*d3*d4 - 2*y2*y5*d3*d5 + (2*y2*d3*p^2*u*u1^2*v)/z - (2*y2*d3*p*q*u^2*u1^2)/z + (2*y2*d3*p*q*u1^2*v^2)/z + (2*y2*d3*p*u1^2*v*w*y)/z - (2*y2*d3*q^2*u*u1^2*v)/z - (2*y2*d3*q*u*u1^2*w*y)/z - y2*y3*y4*y5*d4^2 - 2*y2*y3*y5*d4*d5 + (2*y2*y3*d4*p^2*u*u1^2*v)/z - (2*y2*y3*d4*p*q*u^2*u1^2)/z + (2*y2*y3*d4*p*q*u1^2*v^2)/z + (2*y2*y3*d4*p*u1^2*v*w*y)/z - (2*y2*y3*d4*q^2*u*u1^2*v)/z - (2*y2*y3*d4*q*u*u1^2*w*y)/z - y2*y3*y4*y5*d5^2 + (2*y2*y3*y4*d5*p^2*u*u1^2*v)/z - (2*y2*y3*y4*d5*p*q*u^2*u1^2)/z + (2*y2*y3*y4*d5*p*q*u1^2*v^2)/z + (2*y2*y3*y4*d5*p*u1^2*v*w*y)/z - (2*y2*y3*y4*d5*q^2*u*u1^2*v)/z - (2*y2*y3*y4*d5*q*u*u1^2*w*y)/z + (p^3*u1^2*v*w*y1)/z - (y2*y3*y4*y5*p^3*u1^2*v)/z - (p^2*q*u*u1^2*w*y1)/z + (y2*y3*y4*y5*p^2*q*u*u1^2)/z - (2*p^2*u*u1^2*v*y*y1)/z - p^2*w*y1 + (p*q^2*u1^2*v*w*y1)/z - (y2*y3*y4*y5*p*q^2*u1^2*v)/z + (2*p*q*u^2*u1^2*y*y1)/z - (2*p*q*u1^2*v^2*y*y1)/z + 2*p*u*y*y1 - (p*u1^2*v*w*y^2*y1)/z - (y2*y3*y4*y5*p*u1^2*v*y^2)/z - (q^3*u*u1^2*w*y1)/z + (y2*y3*y4*y5*q^3*u*u1^2)/z + (2*q^2*u*u1^2*v*y*y1)/z - q^2*w*y1 + (q*u*u1^2*w*y^2*y1)/z + (y2*y3*y4*y5*q*u*u1^2*y^2)/z + 2*q*v*y*y1 + w*y^2*y1;
k51_1= (Q1*u1^2*u5*y3)/(2*a1) + a5*u1*w*y2*y4 + (a2^2*u1^2*u5*y3)/(2*a1) - (a3^2*u1^2*u5*y3)/(2*a1) + (a4^2*u1^2*u5*y3)/(2*a1) - (a5^2*u1^2*u5*y3)/(2*a1) - (d2^2*u1^2*u5*y3)/(2*a1) - (d3^2*u1^2*u5*y3)/(2*a1) - (d4^2*u1^2*u5*y3)/(2*a1) - (d5^2*u1^2*u5*y3)/(2*a1) - (d3*d4*u1^2*u5)/a1 + a3*u1*u3*u5*y1 - (d2*d3*u1^2*u5*y2*y3)/a1 - (d4*d5*u1^2*u5*y3*y4)/a1 - a5*u1* y1*y3*y5 - (a3*a5*u1^2*u3*y5)/a1 - (d2*d4*u1^2*u5*y2)/a1 - (d3*d5*u1^2*u5*y4)/a1 - (d2*d5*u1^2*u5*y2*y4)/a1;
k51_2= u1*u5*y*y2*y3 - d1*u1*u5*y2*y3 - (Q2*a5*u1^2*y4)/(2*a1) + (a5*d4*u1^2*y5)/a1 - d2*u1*u5*y1*y2*y3 - d3*u1*u5*y1*y3 - (a3*d3*u1^2*u3*u5)/a1 - (a4*d5*u1^2*u4*u5)/a1 + (a5*d5*u1^2*y4*y5)/a1 - d4*u1*u5*y1 - (a2*d2*u1^2*u2*u5*y3)/a1 - (a3*d2*u1^2*u3*u5*y2)/a1 - d5*u1*u5*y1*y4 + (a5*d3*u1^2* y3*y5)/a1 + (a5*d2*u1^2* y2*y3*y5)/a1;
k51_3= d1*u1*u2*u3*u5 - u1*u2*u3*u5*y + d2*u1*u2*u3*u5*y1 - (a2*d3*u1^2*u3*u5)/a1 - (a2*d2*u1^2*u3*u5*y2)/a1 - (a3*d2*u1^2*u2*u5*y3)/a1 - (a5*d2*u1^2*u2*u3*y5)/a1;
k52_1= u1*u5*y*y2*y4 - d1*u1*u5*y2*y4 - (Q2*a5*u1^2*y3)/(2*a1) - d2*u1*u5*y1*y2*y4 - d3*u1*u5*y1*y4 + (a5*d5*u1^2*y3*y5)/a1 - (a2*d2*u1^2*u2*u5*y4)/a1 - (a3*d4*u1^2*u3*u5*y4)/a1 - (a4*d4*u1^2*u4*u5*y3)/a1 + (a5*d4*u1^2*y3*y4*y5)/a1 - d5*u1*u5*y1*y3 - (a4*d3*u1^2*u4*u5)/a1 - (a3*d5*u1^2*u3*u5)/a1 + (a5*d3*u1^2*y4*y5)/a1 - d4*u1*u5*y1*y3*y4 - (a4*d2*u1^2*u4*u5*y2)/a1 + (a5*d2*u1^2*y2*y4*y5)/a1;
k52_2= (a4^2*u1^2*u5*y4)/(2*a1) - a5*u1*w*y2*y3 - (a2^2*u1^2*u5*y4)/(2*a1) - (a3^2*u1^2*u5*y4)/(2*a1) - (Q1*u1^2*u5*y4)/(2*a1) + (a5^2*u1^2*u5*y4)/(2*a1) + (d2^2*u1^2*u5*y4)/(2*a1) + (d3^2*u1^2*u5*y4)/(2*a1) + (d4^2*u1^2*u5*y4)/(2*a1) + (d5^2*u1^2*u5*y4)/(2*a1) + (d4*d5*u1^2*u5)/a1 + (a4*a5*u1^2*u4*y5)/a1 + (d2*d3*u1^2*u5*y2*y4)/a1 + (d3*d4*u1^2*u5*y3*y4)/a1 - a4*u1* u4*u5*y1 + a5*u1* y1*y4*y5 + (d3*d5*u1^2*u5*y3)/a1 + (d2*d4*u1^2*u5*y2*y3*y4)/a1 + (d2*d5*u1^2*u5*y2*y3)/a1;
k52_3= a5*u1*u2*u3*w - (a2*a3*u1^2*u5*y4)/a1 - (d2*d5*u1^2*u2*u3*u5)/a1 - (d2*d4*u1^2*u2*u3*u5*y4)/a1;
k53_1= u1*u5*y*y2 - d1*u1*u5*y2 - d2*u1*u5*y1*y2 - d3*u1*u5*y1 - d4*u1*u5*y1*y3 - (Q2*a5*u1^2*y3*y4)/(2*a1) - (a2*d2*u1^2*u2*u5)/a1 - (a3*d4*u1^2*u3*u5)/a1 - (a3*d5*u1^2*u3*u5*y4)/a1 - (a4*d5*u1^2*u4*u5*y3)/a1 + (a5*d5*u1^2*y3*y4*y5)/a1 + (a5*d4*u1^2*y3*y5)/a1 + (a5*d3*u1^2*y5)/a1 - d5*u1*u5*y1*y3*y4 + (a5*d2*u1^2*y2*y5)/a1;
k53_2=(a5^2*u1^2*u5)/(2*a1) - (a3^2*u1^2*u5)/(2*a1) - (a4^2*u1^2*u5)/(2*a1) - (a2^2*u1^2*u5)/(2*a1) + (d2^2*u1^2*u5)/(2*a1) + (d3^2*u1^2*u5)/(2*a1) + (d4^2*u1^2*u5)/(2*a1) + (d5^2*u1^2*u5)/(2*a1) - (Q1*u1^2*u5)/(2*a1) - a5*u1*w*y2*y3*y4 + (d2*d3*u1^2*u5*y2)/a1 + (d3*d4*u1^2*u5*y3)/a1 + (d4*d5*u1^2*u5*y4)/a1 + (d2*d4*u1^2*u5*y2*y3)/a1 + (d3*d5*u1^2*u5*y3*y4)/a1 + a5*u1*y1*y5 + (d2*d5*u1^2*u5*y2*y3*y4)/a1;
k53_3= a5*u1*u2*u3*w*y4 - (a2*a3*u1^2*u5)/a1 - (d2*d4*u1^2*u2*u3*u5)/a1 - (d2*d5*u1^2*u2*u3*u5*y4)/a1;
k54_1= (a3^2*u1^2*u5*y3*y4)/(2*a1) - (a2^2*u1^2*u5*y3*y4)/(2*a1) - a5*u1*w*y2 + (a4^2*u1^2*u5*y3*y4)/(2*a1) + (a5^2*u1^2*u5*y3*y4)/(2*a1) + (d2^2*u1^2*u5*y3*y4)/(2*a1) + (d3^2*u1^2*u5*y3*y4)/(2*a1) + (d4^2*u1^2*u5*y3*y4)/(2*a1) + (d5^2*u1^2*u5*y3*y4)/(2*a1) - (Q1*u1^2*u5*y3*y4)/(2*a1) - (a3*a4*u1^2*u3*u4*u5)/a1 + (a3*a5*u1^2*u3*y4*y5)/a1 + (a4*a5*u1^2*u4*y3*y5)/a1 + (d3*d4*u1^2*u5*y4)/a1 + (d4*d5*u1^2*u5*y3)/a1 + (d3*d5*u1^2*u5)/a1 - a3*u1*u3*u5*y1*y4 - a4*u1*u4*u5*y1*y3 + a5*u1* y1*y3*y4*y5 + (d2*d5*u1^2*u5*y2)/a1 + (d2*d3*u1^2*u5*y2*y3*y4)/a1 + (d2*d4*u1^2* u5*y2*y4)/a1;
k54_2=(Q2*a5*u1^2)/(2*a1) - (a5*d5*u1^2*y5)/a1 + d1*u1*u5*y2*y3*y4 - u1*u5*y*y2*y3*y4 + (a4*d4*u1^2*u4*u5)/a1 - (a5*d4*u1^2*y4*y5)/a1 + (a3*d3*u1^2*u3*u5*y4)/a1 + (a4*d3*u1^2*u4*u5*y3)/a1 - (a5*d3*u1^2*y3*y4*y5)/a1 + d4*u1*u5*y1*y4 + d2*u1*u5*y1*y2*y3*y4 + d5*u1*u5*y1 + d3*u1*u5*y1*y3*y4 + (a2*d2*u1^2*u2*u5*y3*y4)/a1 + (a3*d2*u1^2*u3*u5*y2*y4)/a1 + (a4*d2*u1^2*u4*u5*y2*y3)/a1 - (a5*d2*u1^2*y2*y3*y4*y5)/a1;
k54_3= u1*u2*u3*u5*y*y4 - d1*u1*u2*u3*u5*y4 + (a2*d3*u1^2*u3*u5*y4)/a1 - d2*u1*u2*u3*u5*y1*y4 - (a4*d2*u1^2*u2*u3*u4*u5)/a1 + (a2*d2*u1^2*u3*u5*y2*y4)/a1 + (a3*d2*u1^2*u2*u5*y3*y4)/a1 + (a5*d2*u1^2*u2*u3*y4*y5)/a1;
k55_1= u1*u4*y*y2*y5 - d5*u1*u4*w*y2 - d1*u1*u4*y2*y5 - (Q2*a4*u1^2*y3)/(2*a1) - d2*u1*u4*y1*y2*y5 - d3*u1*u4*y1*y5 - (a5*d3*u1^2*u4*u5)/a1 + (a4*d5*u1^2*y3*y5)/a1 - (a2*d2*u1^2*u2*u4*y5)/a1 - (a5*d2*u1^2*u4*u5*y2)/a1 - (a3*d4*u1^2*u3*u4*y5)/a1 - (a5*d4*u1^2*u4*u5*y3)/a1 + (a4*d4*u1^2*y3*y4*y5)/a1 + (a4*d3*u1^2*y4*y5)/a1 - d4*u1*u4*y1*y3*y5 + (a4*d2*u1^2*y2*y4*y5)/a1;
k55_2= (Q2*d5*u1^2*u4)/(2*a1) - (Q1*u1^2*u4*y5)/(2*a1) - a5*u1*u4*u5*y1 - a4*u1*w*y2*y3 - (a2^2*u1^2*u4*y5)/(2*a1) - (a3^2*u1^2*u4*y5)/(2*a1) + (a4^2*u1^2*u4*y5)/(2*a1) + (a5^2*u1^2*u4*y5)/(2*a1) + (d2^2*u1^2*u4*y5)/(2*a1) + (d3^2*u1^2*u4*y5)/(2*a1) + (d4^2*u1^2*u4*y5)/(2*a1) - (d5^2*u1^2*u4*y5)/(2*a1) + (a4*a5*u1^2*u5*y4)/a1 + (d2*d3*u1^2*u4*y2*y5)/a1 + (d3*d4*u1^2*u4*y3*y5)/a1 + a4*u1*y1*y4*y5 + (d2*d4*u1^2*u4*y2*y3*y5)/a1;
k55_3= a4*u1*u2*u3*w - (a2*a3*u1^2*u4*y5)/a1 - (d2*d4*u1^2*u2*u3*u4*y5)/a1;
k56_1= (a3^2*u1^2*u4*y3*y5)/(2*a1) - a5*u1*u4*u5*y1*y3 - (a2^2*u1^2*u4*y3*y5)/(2*a1) - a4*u1*w*y2 + (a4^2*u1^2*u4*y3*y5)/(2*a1) + (a5^2*u1^2*u4*y3*y5)/(2*a1) + (d2^2*u1^2*u4*y3*y5)/(2*a1) + (d3^2*u1^2*u4*y3*y5)/(2*a1) + (d4^2*u1^2*u4*y3*y5)/(2*a1) - (d5^2*u1^2*u4*y3*y5)/(2*a1) + (Q2*d5*u1^2*u4*y3)/(2*a1) - (Q1*u1^2*u4*y3*y5)/(2*a1) - (a3*a5*u1^2*u3*u4*u5)/a1 + (a3*a4*u1^2*u3*y4*y5)/a1 + (a4*a5*u1^2*u5*y3*y4)/a1 + (d3*d4*u1^2*u4*y5)/a1 - a3*u1*u3*u4*y1*y5 + a4*u1*y1*y3*y4*y5 + (d2*d3*u1^2*u4*y2*y3*y5)/a1 + (d2*d4*u1^2* u4*y2*y5)/a1;
k56_2= (Q2*a4*u1^2)/(2*a1) - (a4*d5*u1^2*y5)/a1 + d5*u1*u4*w*y2*y3 + d1*u1*u4*y2*y3*y5 - u1*u4*y*y2*y3*y5 + (a5*d4*u1^2*u4*u5)/a1 - (a4*d4*u1^2*y4*y5)/a1 + (a3*d3*u1^2*u3*u4*y5)/a1 + (a5*d3*u1^2*u4*u5*y3)/a1 - (a4*d3*u1^2*y3*y4*y5)/a1 + d4*u1*u4*y1*y5 + d2*u1*u4*y1*y2*y3*y5 + d3*u1* u4*y1*y3*y5 + (a2*d2*u1^2*u2*u4*y3*y5)/a1 + (a3*d2*u1^2*u3*u4*y2*y5)/a1 + (a5*d2*u1^2*u4*u5*y2*y3)/a1 - (a4*d2*u1^2*y2*y3*y4*y5)/a1;
k56_3= u1*u2*u3*u4*y*y5 - d1*u1*u2*u3*u4*y5 - d5*u1*u2*u3*u4*w + (a2*d3*u1^2*u3*u4*y5)/a1 - d2*u1*u2*u3*u4*y1*y5 - (a5*d2*u1^2*u2*u3*u4*u5)/a1 + (a2*d2*u1^2*u3*u4*y2*y5)/a1 + (a3*d2*u1^2*u2*u4*y3*y5)/a1 + (a4*d2*u1^2*u2*u3*y4*y5)/a1;
k57_1= d5*u1*u3*u4*u5*y1 + (Q2*a5*u1^2*u3*u4)/(2*a1) - (a4*d4*u1^2*u3*u5)/a1 - (a3*d5*u1^2*u4*u5*y3)/a1 - (a4*d5*u1^2*u3*u5*y4)/a1 - (a5*d5*u1^2*u3*u4*y5)/a1;
k57_2= a5*u1*u3*u4*w*y2 - (a3*a4*u1^2*u5)/a1 - (d3*d5*u1^2*u3*u4*u5)/a1 - (d2*d5*u1^2*u3*u4*u5*y2)/a1;
k57_3= a5*u1*u2*u4*w*y3 - (a2*a4*u1^2*u5)/a1 - (d2*d5*u1^2*u2*u4*u5*y3)/a1;
k58_1= (a2^2*u1^2*u3*u4*u5)/(2*a1) - (a3^2*u1^2*u3*u4*u5)/(2*a1) - (a4^2*u1^2*u3*u4*u5)/(2*a1) - (a5^2*u1^2*u3*u4*u5)/(2*a1) - (d2^2*u1^2*u3*u4*u5)/(2*a1) - (d3^2*u1^2*u3*u4*u5)/(2*a1) + (d4^2*u1^2*u3*u4*u5)/(2*a1) - (d5^2*u1^2*u3*u4*u5)/(2*a1) + (Q1*u1^2*u3*u4*u5)/(2*a1) + (a3*a4*u1^2*u5*y3*y4)/a1 + (a3*a5*u1^2*u4*y3*y5)/a1 + (a4*a5*u1^2*u3*y4*y5)/a1 - a3*u1*u4*u5*y1*y3 - a5*u1*u3*u4*y1*y5 - a4*u1*u3*u5*y1*y2^2*y4 - (d2*d3*u1^2*u3*u4*u5*y2)/a1;
k58_2= u1*u3*u4*u5*y*y2 - d1*u1*u3*u4*u5*y2 + (a3*d4*u1^2*u4*u5)/a1 + (a3*d3*u1^2*u4*u5*y3)/a1 + (a4*d3*u1^2*u3*u5*y4)/a1 + (a5*d3*u1^2*u3*u4*y5)/a1 - d2*u1*u3*u4*u5*y1*y2 - d3*u1*u3*u4*u5*y1 - (a2*d2*u1^2*u2*u3*u4*u5)/a1 + (a3*d2*u1^2*u4*u5*y2*y3)/a1 + (a4*d2*u1^2*u3*u5*y2*y4)/a1 + (a5*d2*u1^2*u3*u4*y2*y5)/a1;
k58_3= u1*u2*u4*u5*y*y3 - d1*u1*u2*u4*u5*y3 + (a2*d4*u1^2*u4*u5)/a1 + (a2*d3*u1^2*u4*u5*y3)/a1 - d2*u1*u2*u4*u5*y1*y3 - (a3*d2*u1^2*u2*u3*u4*u5)/a1 + (a2*d2*u1^2*u4*u5*y2*y3)/a1 + (a4*d2*u1^2*u2*u5*y3*y4)/a1 + (a5*d2*u1^2*u2*u4*y3*y5)/a1;
k59_1= (Q2*d4*u1^2*u3)/(2*a1) - a3*u1*w*y2 - a4*u1*u3*u4*y1*y5 - a5*u1*u3*u5*y1*y4 - (a2^2*u1^2*u3*y4*y5)/(2*a1) + (a3^2*u1^2*u3*y4*y5)/(2*a1) + (a4^2*u1^2*u3*y4*y5)/(2*a1) + (a5^2*u1^2*u3*y4*y5)/(2*a1) + (d2^2*u1^2*u3*y4*y5)/(2*a1) + (d3^2*u1^2*u3*y4*y5)/(2*a1) - (d4^2*u1^2*u3*y4*y5)/(2*a1) - (d5^2*u1^2*u3*y4*y5)/(2*a1) + (Q2*d5*u1^2*u3*y4)/(2*a1) - (d4*d5*u1^2*u3*y5)/a1 - (Q1*u1^2*u3*y4*y5)/(2*a1) - (a4*a5*u1^2*u3*u4*u5)/a1 + (a3*a4*u1^2*u4*y3*y5)/a1 + (a3*a5*u1^2*u5*y3*y4)/a1 + a3*u1*y1*y3*y4*y5 - a4*c5*u1*u2^2*u3*u5*y1*y4 + (d2*d3*u1^2*u3*y2*y4*y5)/a1;
k59_2= (Q2*a3*u1^2)/(2*a1) - (a3*d5*u1^2*y5)/a1 + d4*u1*u3*w*y2 + d5*u1*u3*w*y2*y4 + d1*u1*u3*y2*y4*y5 - u1*u3*y*y2*y4*y5 - (a3*d4*u1^2*y4*y5)/a1 + (a4*d3*u1^2*u3*u4*y5)/a1 + (a5*d3*u1^2*u3*u5*y4)/a1 - (a3*d3*u1^2*y3*y4*y5)/a1 + d2*u1*u3*y1*y2*y4*y5 + d3*u1*u3*y1*y4*y5 + (a2*d2*u1^2*u2*u3*y4*y5)/a1 + (a4*d2*u1^2*u3*u4*y2*y5)/a1 + (a5*d2*u1^2*u3*u5*y2*y4)/a1 - (a3*d2*u1^2*y2*y3*y4*y5)/a1;
k59_3= d3*u1*u2*w + (Q2*a2*u1^2)/(2*a1) - (a2*d5*u1^2*y5)/a1 + d4*u1*u2*w*y3 + d5*u1*u2*w*y3*y4 + d1*u1*u2*y3*y4*y5 - u1*u2*y*y3*y4*y5 - (a2*d4*u1^2*y4*y5)/a1 - (a2*d3*u1^2*y3*y4*y5)/a1 + d2*u1*u2*y1*y3*y4*y5 + (a3*d2*u1^2*u2*u3*y4*y5)/a1 + (a4*d2*u1^2*u2*u4*y3*y5)/a1 + (a5*d2*u1^2*u2*u5*y3*y4)/a1 - (a2*d2*u1^2*y2*y3*y4*y5)/a1;
k61_1= d1*u1*u5*y3 - u1*u5*y*y3 + d2*u1*u5*y1*y3 + d3*u1*u5*y1*y2*y3 + d4*u1*u5*y1*y2 + (Q2*a5*u1^2*y2*y4)/(2*a1) + (a2*d3*u1^2*u2*u5*y3)/a1 + (a3*d3*u1^2*u3*u5*y2)/a1 + (a4*d5*u1^2*u4*u5*y2)/a1 - (a5*d5*u1^2*y2*y4*y5)/a1 + (a3*d2*u1^2*u3*u5)/a1 + (a2*d4*u1^2*u2*u5)/a1 - (a5*d4*u1^2*y2*y5)/a1 + d5*u1*u5*y1*y2*y4 - (a5*d2*u1^2*y3*y5)/a1 + (a2*d5*u1^2*u2*u5*y4)/a1 - (a5*d3*u1^2*y2*y3*y5)/a1;
k61_2= a5*u1*w*y4 + a2*u1*u2*u5*y1*y3 + a3*u1*u3*u5*y1*y2 - (a2^2*u1^2*u5*y2*y3)/(2*a1) - (a3^2*u1^2*u5*y2*y3)/(2*a1) + (a4^2*u1^2*u5*y2*y3)/(2*a1) - (a5^2*u1^2*u5*y2*y3)/(2*a1) - (d2^2*u1^2*u5*y2*y3)/(2*a1) - (d3^2*u1^2*u5*y2*y3)/(2*a1) - (d4^2*u1^2*u5*y2*y3)/(2*a1) - (d5^2*u1^2*u5*y2*y3)/(2*a1) + (Q1*u1^2*u5*y2*y3)/(2*a1) + (a2*a3*u1^2*u2*u3*u5)/a1 - (d2*d3*u1^2*u5*y3)/a1 - (d3*d4*u1^2*u5*y2)/a1 - (d2*d4*u1^2*u5)/a1 - a5*u1*y1*y2*y3*y5 - (d2*d5*u1^2*u5*y4)/a1 - (d4*d5*u1^2*u5*y2*y3*y4)/a1 - (a2*a5*u1^2*u2*y3*y5)/a1 - (a3*a5*u1^2*u3*y2*y5)/a1 - (d3*d5*u1^2*u5*y2*y4)/a1;
k61_3= a2*u1*u3*u5*y1*y2 + a3*u1*u2*u5*y1*y3 + a5*u1*u2*u3*y1*y5 + (a2^2*u1^2*u2*u3*u5)/(2*a1) + (a3^2*u1^2*u2*u3*u5)/(2*a1) - (a4^2*u1^2*u2*u3*u5)/(2*a1) + (a5^2*u1^2*u2*u3*u5)/(2*a1) + (d2^2*u1^2*u2*u3*u5)/(2*a1) - (d3^2*u1^2*u2*u3*u5)/(2*a1) + (d4^2*u1^2*u2*u3*u5)/(2*a1) + (d5^2*u1^2*u2*u3*u5)/(2*a1) - (Q1*u1^2*u2*u3*u5)/(2*a1) - (a2*a3*u1^2*u5*y2*y3)/a1 - (a2*a5*u1^2*u3*y2*y5)/a1 + (d4*d5*u1^2*u2*u3*u5*y4)/a1 - (a3*a5*u1^2*u2*y3*y5)/a1;
k62_1= a5*u1*w*y3 + a2*u1*u2*u5*y1*y4 - (a2^2*u1^2*u5*y2*y4)/(2*a1) + (a3^2*u1^2*u5*y2*y4)/(2*a1) - (a4^2*u1^2*u5*y2*y4)/(2*a1) - (a5^2*u1^2*u5*y2*y4)/(2*a1) - (d2^2*u1^2*u5*y2*y4)/(2*a1) - (d3^2*u1^2*u5*y2*y4)/(2*a1) - (d4^2*u1^2*u5*y2*y4)/(2*a1) - (d5^2*u1^2*u5*y2*y4)/(2*a1) + (Q1*u1^2*u5*y2*y4)/(2*a1) - (a4*a5*u1^2*u4*y2*y5)/a1 - (d2*d3*u1^2*u5*y4)/a1 - (d4*d5*u1^2*u5*y2)/a1 + a4*u1*u4*u5*y1*y2 - a5*u1*y1*y2*y4*y5 - (d2*d5*u1^2*u5*y3)/a1 - (d3*d4*u1^2*u5*y2*y3*y4)/a1 + (a2*a4*u1^2*u2*u4*u5)/a1 - (a2*a5*u1^2*u2*y4*y5)/a1 - (d2*d4*u1^2*u5*y3*y4)/a1 - (d3*d5*u1^2*u5*y2*y3)/a1;
k62_2= u1*u5*y*y4 - d1*u1*u5*y4 - d2*u1*u5*y1*y4 - d3*u1*u5*y1*y2*y4 - (Q2*a5*u1^2*y2*y3)/(2*a1) - (a2*d3*u1^2*u2*u5*y4)/a1 + (a5*d5*u1^2*y2*y3*y5)/a1 - d4*u1*u5*y1*y2*y3*y4 - d5*u1*u5*y1*y2*y3 - (a4*d2*u1^2*u4*u5)/a1 + (a5*d2*u1^2*y4*y5)/a1 - (a2*d4*u1^2*u2*u5*y3*y4)/a1 - (a3*d4*u1^2*u3*u5*y2*y4)/a1 - (a4*d4*u1^2*u4*u5*y2*y3)/a1 + (a5*d4*u1^2*y2*y3*y4*y5)/a1 - (a2*d5*u1^2*u2*u5*y3)/a1 - (a4*d3*u1^2*u4*u5*y2)/a1 - (a3*d5*u1^2*u3*u5*y2)/a1 + (a5*d3*u1^2*y2*y4*y5)/a1;
k62_3= d5*u1*u2*u3*u5*y1 + (Q2*a5*u1^2*u2*u3)/(2*a1) - (a2*d5*u1^2*u3*u5*y2)/a1 - (a3*d3*u1^2*u2*u5*y4)/a1 - (a5*d5*u1^2*u2*u3*y5)/a1 + d4*u1*u2*u3*u5*y1*y4 + (a4*d4*u1^2*u2*u3*u4*u5)/a1 - (a2*d4*u1^2*u3*u5*y2*y4)/a1 - (a3*d4*u1^2*u2*u5*y3*y4)/a1 - (a5*d4*u1^2*u2*u3*y4*y5)/a1 - (a3*d5*u1^2*u2*u5*y3)/a1;
k63_1= (Q1*u1^2*u5*y2)/(2*a1) + a2*u1*u2*u5*y1 + a5*u1*w*y3*y4 - (a2^2*u1^2*u5*y2)/(2*a1) + (a3^2*u1^2*u5*y2)/(2*a1) + (a4^2*u1^2*u5*y2)/(2*a1) - (a5^2*u1^2*u5*y2)/(2*a1) - (d2^2*u1^2*u5*y2)/(2*a1) - (d3^2*u1^2*u5*y2)/(2*a1) - (d4^2*u1^2*u5*y2)/(2*a1) - (d5^2*u1^2*u5*y2)/(2*a1) - (d2*d3*u1^2*u5)/a1 - (d3*d4*u1^2*u5*y2*y3)/a1 - (d4*d5*u1^2*u5*y2*y4)/a1 - a5*u1* y1*y2*y5 - (d2*d4*u1^2*u5*y3)/a1 - (a2*a5*u1^2*u2*y5)/a1 - (d3*d5*u1^2*u5*y2*y3*y4)/a1 - (d2*d5*u1^2*u5*y3*y4)/a1;
k63_2= u1*u5*y - d1*u1*u5 - d2*u1*u5*y1 - d3*u1*u5*y1*y2 - d4*u1*u5*y1*y2*y3 - (a2*d3*u1^2*u2*u5)/a1 - (Q2*a5*u1^2*y2*y3*y4)/(2*a1) - (a2*d4*u1^2*u2*u5*y3)/a1 - (a3*d4*u1^2*u3*u5*y2)/a1 - d5*u1*u5*y1*y2*y3*y4 + (a5*d3*u1^2*y2* y5)/a1 - (a2*d5*u1^2*u2*u5*y3*y4)/a1 - (a3*d5*u1^2*u3*u5*y2*y4)/a1 - (a4*d5*u1^2*u4*u5*y2*y3)/a1 + (a5*d5*u1^2*y2*y3*y4*y5)/a1 + (a5*d2*u1^2*y5)/a1 + (a5*d4*u1^2*y2*y3*y5)/a1;
k63_3= d4*u1*u2*u3*u5*y1 - (a3*d3*u1^2*u2*u5)/a1 + (Q2*a5*u1^2*u2*u3*y4)/(2*a1) - (a2*d4*u1^2*u3*u5*y2)/a1 - (a3*d4*u1^2*u2*u5*y3)/a1 + d5*u1*u2*u3*u5*y1*y4 + (a4*d5*u1^2*u2*u3*u4*u5)/a1 - (a2*d5*u1^2*u3*u5*y2*y4)/a1 - (a3*d5*u1^2*u2*u5*y3*y4)/a1 - (a5*d5*u1^2*u2*u3*y4*y5)/a1 - (a5*d4*u1^2*u2*u3*y5)/a1;
k64_1= u1*u5*y*y3*y4 - d1*u1*u5*y3*y4 - (Q2*a5*u1^2*y2)/(2*a1) - d2*u1*u5*y1*y3*y4 + (a5*d5*u1^2*y2*y5)/a1 - (a4*d4*u1^2*u4*u5*y2)/a1 + (a5*d4*u1^2*y2*y4*y5)/a1 - d5*u1*u5*y1*y2 - d3*u1*u5*y1*y2*y3*y4 - d4*u1* u5*y1*y2*y4 - (a2*d5*u1^2*u2*u5)/a1 - (a2*d3*u1^2*u2*u5*y3*y4)/a1 - (a3*d3*u1^2*u3*u5*y2*y4)/a1 - (a4*d3*u1^2*u4*u5*y2*y3)/a1 + (a5*d3*u1^2*y2*y3*y4*y5)/a1 - (a3*d2*u1^2*u3*u5*y4)/a1 - (a2*d4*u1^2*u2*u5*y4)/a1 - (a4*d2*u1^2*u4*u5*y3)/a1 + (a5*d2*u1^2*y3*y4*y5)/a1;
k64_2=(u5*a2^2*u1^2*y2*y3*y4)/(2*a1) - (u5*a2*a3*u1^2*u2*u3*y4)/a1 - (u5*a2*a4*u1^2*u2*u4*y3)/a1 + (y5*a2*a5*u1^2*u2*y3*y4)/a1 - u5*y1*a2*u1*u2*y3*y4 + (u5*a3^2*u1^2*y2*y3*y4)/(2*a1) - (u5*a3*a4*u1^2*u3*u4*y2)/a1 + (y5*a3*a5*u1^2*u3*y2*y4)/a1 - u5*y1*a3*u1*u3*y2*y4 + (u5*a4^2*u1^2*y2*y3*y4)/(2*a1) + (y5*a4*a5*u1^2*u4*y2*y3)/a1 - u5*y1*a4*u1*u4*y2*y3 + (u5*a5^2*u1^2*y2*y3*y4)/(2*a1) + y1*y5*a5*u1*y2*y3*y4 - w*a5*u1 + (u5*d2^2*u1^2*y2*y3*y4)/(2*a1) + (u5*d2*d3*u1^2*y3*y4)/a1 + (u5*d2*d4*u1^2*y4)/a1 + (u5*d2*d5*u1^2)/a1 + (u5*d3^2*u1^2*y2*y3*y4)/(2*a1) + (u5*d3*d4*u1^2*y2*y4)/a1 + (u5*d3*d5*u1^2*y2)/a1 + (u5*d4^2*u1^2*y2*y3*y4)/(2*a1) + (u5*d4*d5*u1^2*y2*y3)/a1 + (u5*d5^2*u1^2*y2*y3*y4)/(2*a1) - (Q1*u5*u1^2*y2*y3*y4)/(2*a1);
k64_3=(Q1*u1^2*u2*u3*u5*y4)/(2*a1) + a4*u1*u2*u3*u4*u5*y1 - a2*u1*u3*u5*y1*y2*y4 - a3*u1*u2*u5*y1*y3*y4 - a5*u1*u2*u3*y1*y4*y5 - (a2^2*u1^2*u2*u3*u5*y4)/(2*a1) - (a3^2*u1^2*u2*u3*u5*y4)/(2*a1) - (a4^2*u1^2*u2*u3*u5*y4)/(2*a1) - (a5^2*u1^2*u2*u3*u5*y4)/(2*a1) - (d2^2*u1^2*u2*u3*u5*y4)/(2*a1) + (d3^2*u1^2*u2*u3*u5*y4)/(2*a1) - (d4^2*u1^2*u2*u3*u5*y4)/(2*a1) - (d5^2*u1^2*u2*u3*u5*y4)/(2*a1) - (a2*a4*u1^2*u3*u4*u5*y2)/a1 - (a3*a4*u1^2*u2*u4*u5*y3)/a1 - (a4*a5*u1^2*u2*u3*u4*y5)/a1 + (a2*a3*u1^2*u5*y2*y3*y4)/a1 + (a2*a5*u1^2*u3*y2*y4*y5)/a1 + (a3*a5*u1^2*u2*y3*y4*y5)/a1 - (d4*d5*u1^2*u2*u3*u5)/a1;
k65_1= a4*u1*w*y3 + a2*u1*u2*u4*y1*y5 + a5*u1*u4*u5*y1*y2 - (a2^2*u1^2*u4*y2*y5)/(2*a1) + (a3^2*u1^2*u4*y2*y5)/(2*a1) - (a4^2*u1^2*u4*y2*y5)/(2*a1) - (a5^2*u1^2*u4*y2*y5)/(2*a1) - (d2^2*u1^2*u4*y2*y5)/(2*a1) - (d3^2*u1^2*u4*y2*y5)/(2*a1) - (d4^2*u1^2*u4*y2*y5)/(2*a1) + (d5^2*u1^2*u4*y2*y5)/(2*a1) - (Q2*d5*u1^2*u4*y2)/(2*a1) + (Q1*u1^2*u4*y2*y5)/(2*a1) + (a2*a5*u1^2*u2*u4*u5)/a1 - (a4*a5*u1^2*u5*y2*y4)/a1 - (d2*d3*u1^2*u4*y5)/a1 - a4*u1*y1*y2*y4*y5 - (d3*d4*u1^2*u4*y2*y3*y5)/a1 - (a2*a4*u1^2*u2*y4*y5)/a1 - (d2*d4*u1^2*u4*y3*y5)/a1;
k65_2= u1*u4*y*y5 - d1*u1*u4*y5 - d5*u1*u4*w - d2*u1*u4*y1*y5 - d3*u1*u4*y1*y2*y5 - (Q2*a4*u1^2*y2*y3)/(2*a1) - (a5*d2*u1^2*u4*u5)/a1 - (a2*d3*u1^2*u2*u4*y5)/a1 - (a5*d3*u1^2*u4*u5*y2)/a1 + (a4*d5*u1^2*y2*y3*y5)/a1 - d4*u1*u4*y1*y2*y3*y5 + (a4*d2*u1^2*y4*y5)/a1 - (a2*d4*u1^2*u2*u4*y3*y5)/a1 - (a3*d4*u1^2*u3*u4*y2*y5)/a1 - (a5*d4*u1^2*u4*u5*y2*y3)/a1 + (a4*d4*u1^2*y2*y3*y4*y5)/a1 + (a4*d3*u1^2*y2*y4*y5)/a1;
k65_3=(Q2*a4*u1^2*u2*u3)/(2*a1) - (a3*d3*u1^2*u2*u4*y5)/a1 - (a4*d5*u1^2*u2*u3*y5)/a1 + d4*u1*u2*u3*u4*y1*y5 + (a5*d4*u1^2*u2*u3*u4*u5)/a1 - (a2*d4*u1^2*u3*u4*y2*y5)/a1 - (a3*d4*u1^2*u2*u4*y3*y5)/a1 - (a4*d4*u1^2*u2*u3*y4*y5)/a1;
k66_1= u1*u4*y*y3*y5 - d5*u1*u4*w*y3 - d1*u1*u4*y3*y5 - (Q2*a4*u1^2*y2)/(2*a1) - d2*u1*u4*y1*y3*y5 + (a4*d5*u1^2*y2*y5)/a1 - (a5*d2*u1^2*u4*u5*y3)/a1 - (a5*d4*u1^2*u4*u5*y2)/a1 + (a4*d4*u1^2*y2*y4*y5)/a1 - d3*u1*u4*y1*y2*y3*y5 - d4*u1*u4*y1*y2*y5 - (a2*d3*u1^2*u2*u4*y3*y5)/a1 - (a3*d3*u1^2*u3*u4*y2*y5)/a1 - (a5*d3*u1^2*u4*u5*y2*y3)/a1 + (a4*d3*u1^2*y2*y3*y4*y5)/a1 - (a3*d2*u1^2*u3*u4*y5)/a1 - (a2*d4*u1^2*u2*u4*y5)/a1 + (a4*d2*u1^2*y3*y4*y5)/a1;
k66_2=(Q2*d5*u1^2*u4*y2*y3)/(2*a1) - a4*u1*w - (Q1*u1^2*u4*y2*y3*y5)/(2*a1) - a2*u1*u2*u4*y1*y3*y5 - a3*u1*u3*u4*y1*y2*y5 - a5*u1*u4*u5*y1*y2*y3 + a4*u1*y1*y2*y3*y4*y5 + (a2^2*u1^2*u4*y2*y3*y5)/(2*a1) + (a3^2*u1^2*u4*y2*y3*y5)/(2*a1) + (a4^2*u1^2*u4*y2*y3*y5)/(2*a1) + (a5^2*u1^2*u4*y2*y3*y5)/(2*a1) + (d2^2*u1^2*u4*y2*y3*y5)/(2*a1) + (d3^2*u1^2*u4*y2*y3*y5)/(2*a1) + (d4^2*u1^2*u4*y2*y3*y5)/(2*a1) - (d5^2*u1^2*u4*y2*y3*y5)/(2*a1) + (d2*d4*u1^2*u4*y5)/a1 - (a2*a3*u1^2*u2*u3*u4*y5)/a1 - (a2*a5*u1^2*u2*u4*u5*y3)/a1 - (a3*a5*u1^2*u3*u4*u5*y2)/a1 + (a2*a4*u1^2*u2*y3*y4*y5)/a1 + (a3*a4*u1^2*u3*y2*y4*y5)/a1 + (a4*a5*u1^2*u5*y2*y3*y4)/a1 + (d2*d3*u1^2*u4*y3*y5)/a1 + (d3*d4*u1^2*u4*y2*y5)/a1;
k66_3= (Q1*u1^2*u2*u3*u4*y5)/(2*a1) - (Q2*d5*u1^2*u2*u3*u4)/(2*a1) + a5*u1*u2*u3*u4*u5*y1 - a2*u1*u3*u4*y1*y2*y5 - a3*u1*u2*u4*y1*y3*y5 - a4*u1*u2*u3*y1*y4*y5 - (a2^2*u1^2*u2*u3*u4*y5)/(2*a1) - (a3^2*u1^2*u2*u3*u4*y5)/(2*a1) - (a4^2*u1^2*u2*u3*u4*y5)/(2*a1) - (a5^2*u1^2*u2*u3*u4*y5)/(2*a1) - (d2^2*u1^2*u2*u3*u4*y5)/(2*a1) + (d3^2*u1^2*u2*u3*u4*y5)/(2*a1) - (d4^2*u1^2*u2*u3*u4*y5)/(2*a1) + (d5^2*u1^2*u2*u3*u4*y5)/(2*a1) - (a2*a5*u1^2*u3*u4*u5*y2)/a1 - (a3*a5*u1^2*u2*u4*u5*y3)/a1 - (a4*a5*u1^2*u2*u3*u5*y4)/a1 + (a2*a3*u1^2*u4*y2*y3*y5)/a1 + (a2*a4*u1^2*u3*y2*y4*y5)/a1 + (a3*a4*u1^2*u2*y3*y4*y5)/a1;
k67_1= (a3*a4*u1^2*u5*y2)/a1 - a5*u1*u3*u4*w + (d2*d5*u1^2*u3*u4*u5)/a1 + (d3*d5*u1^2*u3*u4*u5*y2)/a1;
k67_2=(Q2*a5*u1^2*u3*u4*y2)/(2*a1) - (a4*d4*u1^2*u3*u5*y2)/a1 + d5*u1*u3*u4*u5*y1*y2 + (a2*d5*u1^2*u2*u3*u4*u5)/a1 - (a3*d5*u1^2*u4*u5*y2*y3)/a1 - (a4*d5*u1^2*u3*u5*y2*y4)/a1 - (a5*d5*u1^2*u3*u4*y2*y5)/a1;
k67_3= (Q2*a5*u1^2*u2*u4*y3)/(2*a1) - (a4*d3*u1^2*u2*u5)/a1 - (a4*d4*u1^2*u2*u5*y3)/a1 + d5*u1*u2*u4*u5*y1*y3 + (a3*d5*u1^2*u2*u3*u4*u5)/a1 - (a2*d5*u1^2*u4*u5*y2*y3)/a1 - (a4*d5*u1^2*u2*u5*y3*y4)/a1 - (a5*d5*u1^2*u2*u4*y3*y5)/a1;
k68_1= d1*u1*u3*u4*u5 - u1*u3*u4*u5*y + d2*u1*u3*u4*u5*y1 - (a3*d4*u1^2*u4*u5*y2)/a1 - (a4*d2*u1^2*u3*u5*y4)/a1 - (a5*d2*u1^2*u3*u4*y5)/a1 + d3*u1*u3*u4*u5*y1*y2 + (a2*d3*u1^2*u2*u3*u4*u5)/a1 - (a3*d3*u1^2*u4*u5*y2*y3)/a1 - (a4*d3*u1^2*u3*u5*y2*y4)/a1 - (a5*d3*u1^2*u3*u4*y2*y5)/a1 - (a3*d2*u1^2*u4*u5*y3)/a1;
k68_2=(Q1*u1^2*u3*u4*u5*y2)/(2*a1) + a2*u1*u2*u3*u4*u5*y1 - a3*u1*u4*u5*y1*y2*y3 - a4*u1*u3*u5*y1*y2*y4 - a5*u1*u3*u4*y1*y2*y5 - (a2^2*u1^2*u3*u4*u5*y2)/(2*a1) - (a3^2*u1^2*u3*u4*u5*y2)/(2*a1) - (a4^2*u1^2*u3*u4*u5*y2)/(2*a1) - (a5^2*u1^2*u3*u4*u5*y2)/(2*a1) - (d2^2*u1^2*u3*u4*u5*y2)/(2*a1) - (d3^2*u1^2*u3*u4*u5*y2)/(2*a1) + (d4^2*u1^2*u3*u4*u5*y2)/(2*a1) - (d5^2*u1^2*u3*u4*u5*y2)/(2*a1) - (a2*a3*u1^2*u2*u4*u5*y3)/a1 - (a2*a4*u1^2*u2*u3*u5*y4)/a1 - (a2*a5*u1^2*u2*u3*u4*y5)/a1 + (a3*a4*u1^2*u5*y2*y3*y4)/a1 + (a3*a5*u1^2*u4*y2*y3*y5)/a1 + (a4*a5*u1^2*u3*y2*y4*y5)/a1 - (d2*d3*u1^2*u3*u4*u5)/a1;
k68_3=(d3*d4*u1^2*u2*u4*u5)/a1 + (Q1*u1^2*u2*u4*u5*y3)/(2*a1) + a3*u1*u2*u3*u4*u5*y1 - a2*u1*u4*u5*y1*y2*y3 - a4*u1*u2*u5*y1*y3*y4 - a5*u1*u2*u4*y1*y3*y5 - (a2^2*u1^2*u2*u4*u5*y3)/(2*a1) - (a3^2*u1^2*u2*u4*u5*y3)/(2*a1) - (a4^2*u1^2*u2*u4*u5*y3)/(2*a1) - (a5^2*u1^2*u2*u4*u5*y3)/(2*a1) - (d2^2*u1^2*u2*u4*u5*y3)/(2*a1) + (d3^2*u1^2*u2*u4*u5*y3)/(2*a1) + (d4^2*u1^2*u2*u4*u5*y3)/(2*a1) - (d5^2*u1^2*u2*u4*u5*y3)/(2*a1) - (a2*a3*u1^2*u3*u4*u5*y2)/a1 - (a3*a4*u1^2*u2*u3*u5*y4)/a1 - (a3*a5*u1^2*u2*u3*u4*y5)/a1 + (a2*a4*u1^2*u5*y2*y3*y4)/a1 + (a2*a5*u1^2*u4*y2*y3*y5)/a1 + (a4*a5*u1^2*u2*y3*y4*y5)/a1;
k69_1= u1*u3*y*y4*y5 - (Q2*a3*u1^2*y2)/(2*a1) - d5*u1*u3*w*y4 - d1*u1*u3*y4*y5 - d4*u1*u3*w - d2*u1*u3*y1*y4*y5 + (a3*d5*u1^2*y2*y5)/a1 - (a4*d2*u1^2*u3*u4*y5)/a1 - (a5*d2*u1^2*u3*u5*y4)/a1 + (a3*d4*u1^2*y2*y4*y5)/a1 - d3*u1*u3*y1*y2*y4*y5 - (a2*d3*u1^2*u2*u3*y4*y5)/a1 - (a4*d3*u1^2*u3*u4*y2*y5)/a1 - (a5*d3*u1^2*u3*u5*y2*y4)/a1 + (a3*d3*u1^2*y2*y3*y4*y5)/a1 + (a3*d2*u1^2*y3*y4*y5)/a1;
k69_2= (Q2*d4*u1^2*u3*y2)/(2*a1) - a3*u1*w + (Q2*d5*u1^2*u3*y2*y4)/(2*a1) - (d4*d5*u1^2*u3*y2*y5)/a1 - (Q1*u1^2*u3*y2*y4*y5)/(2*a1) - a2*u1*u2*u3*y1*y4*y5 - a4*u1*u3*u4*y1*y2*y5 - a5*u1*u3*u5*y1*y2*y4 + a3*u1*y1*y2*y3*y4*y5 + (a2^2*u1^2*u3*y2*y4*y5)/(2*a1) + (a3^2*u1^2*u3*y2*y4*y5)/(2*a1) + (a4^2*u1^2*u3*y2*y4*y5)/(2*a1) + (a5^2*u1^2*u3*y2*y4*y5)/(2*a1) + (d2^2*u1^2*u3*y2*y4*y5)/(2*a1) + (d3^2*u1^2*u3*y2*y4*y5)/(2*a1) - (d4^2*u1^2*u3*y2*y4*y5)/(2*a1) - (d5^2*u1^2*u3*y2*y4*y5)/(2*a1) - (a2*a4*u1^2*u2*u3*u4*y5)/a1 - (a2*a5*u1^2*u2*u3*u5*y4)/a1 - (a4*a5*u1^2*u3*u4*u5*y2)/a1 + (a2*a3*u1^2*u2*y3*y4*y5)/a1 + (a3*a4*u1^2*u4*y2*y3*y5)/a1 + (a3*a5*u1^2*u5*y2*y3*y4)/a1 + (d2*d3*u1^2*u3*y4*y5)/a1;
k69_3= (Q2*d3*u1^2*u2)/(2*a1) - a2*u1*w + (Q2*d4*u1^2*u2*y3)/(2*a1) - (d3*d5*u1^2*u2*y5)/a1 + (Q2*d5*u1^2*u2*y3*y4)/(2*a1) - (d3*d4*u1^2*u2*y4*y5)/a1 - (d4*d5*u1^2*u2*y3*y5)/a1 - (Q1*u1^2*u2*y3*y4*y5)/(2*a1) - a3*u1*u2*u3*y1*y4*y5 - a4*u1*u2*u4*y1*y3*y5 - a5*u1*u2*u5*y1*y3*y4 + a2*u1*y1*y2*y3*y4*y5 + (a2^2*u1^2*u2*y3*y4*y5)/(2*a1) + (a3^2*u1^2*u2*y3*y4*y5)/(2*a1) + (a4^2*u1^2*u2*y3*y4*y5)/(2*a1) + (a5^2*u1^2*u2*y3*y4*y5)/(2*a1) + (d2^2*u1^2*u2*y3*y4*y5)/(2*a1) - (d3^2*u1^2*u2*y3*y4*y5)/(2*a1) - (d4^2*u1^2*u2*y3*y4*y5)/(2*a1) - (d5^2*u1^2*u2*y3*y4*y5)/(2*a1) - (a3*a4*u1^2*u2*u3*u4*y5)/a1 - (a3*a5*u1^2*u2*u3*u5*y4)/a1 - (a4*a5*u1^2*u2*u4*u5*y3)/a1 + (a2*a3*u1^2*u3*y2*y4*y5)/a1 + (a2*a4*u1^2*u4*y2*y3*y5)/a1 + (a2*a5*u1^2*u5*y2*y3*y4)/a1;
f11_1=(k14_3+k19_3-k16_3-k18_3)- (k14_2+k19_2-k16_2-k18_2);
f11_2=2*(k14_1+k19_1-k16_1-k18_1);
f11_3=(k14_2+k19_2-k16_2-k18_2)+ (k14_3+k19_3-k16_3-k18_3);
f12_1=(-2*k13_3+2*k17_3)- (-2*k13_2+2*k17_2);
f12_2=2*(-2*k13_1+2*k17_1);
f12_3=(-2*k13_2+2*k17_2)+(-2*k13_3+2*k17_3);
f13_1=(-k14_3-k16_3+k18_3+k19_3)- (-k14_2-k16_2+k18_2+k19_2);
f13_2=2*(-k14_1-k16_1+k18_1+k19_1);
f13_3=(-k14_2-k16_2+k18_2+k19_2)+ (-k14_3-k16_3+k18_3+k19_3);
f14_1=(2*k15_3-2*k12_3)- (2*k15_2-2*k12_2);
f14_2=2*(2*k15_1-2*k12_1);
f14_3=(2*k15_2-2*k12_2)+ (2*k15_3-2*k12_3);
f15_1=(4*k11_3)- (4*k11_2);
f15_2=2*(4*k11_1);
f15_3=(4*k11_2)+ (4*k11_3);
f16_1=(2*k12_3+2*k15_3)- (2*k12_2+2*k15_2);
f16_2=2*(2*k12_1+2*k15_1);
f16_3=(2*k12_2+2*k15_2)+ (2*k12_3+2*k15_3);
f17_1=(-k14_3+k16_3+k19_3-k18_3)- (-k14_2+k16_2+k19_2-k18_2);
f17_2=2*(-k14_1+k16_1+k19_1-k18_1);
f17_3=(-k14_2+k16_2+k19_2-k18_2)+ (-k14_3+k16_3+k19_3-k18_3);
f18_1=(2*k13_3+2*k17_3)- (2*k13_2+2*k17_2);
f18_2=2*(2*k13_1+2*k17_1);
f18_3=(2*k13_2+2*k17_2)+(2*k13_3+2*k17_3);
f19_1=(k14_3+k16_3+k18_3+k19_3)- (k14_2+k16_2+k18_2+k19_2);
f19_2=2*(k14_1+k16_1+k18_1+k19_1);
f19_3=(k14_2+k16_2+k18_2+k19_2)+ (k14_3+k16_3+k18_3+k19_3);
f21_1=(k24_3+k29_3-k26_3-k28_3)- (k24_2+k29_2-k26_2-k28_2);
f21_2=2*(k24_1+k29_1-k26_1-k28_1);
f21_3=(k24_2+k29_2-k26_2-k28_2)+ (k24_3+k29_3-k26_3-k28_3);
f22_1=(-2*k23_3+2*k27_3)- (-2*k23_2+2*k27_2);
f22_2=2*(-2*k23_1+2*k27_1);
f22_3=(-2*k23_2+2*k27_2)+ (-2*k23_3+2*k27_3);
f23_1=(-k24_3-k26_3+k28_3+k29_3)- (-k24_2-k26_2+k28_2+k29_2);
f23_2=2*(-k24_1-k26_1+k28_1+k29_1);
f23_3=(-k24_2-k26_2+k28_2+k29_2)+(-k24_3-k26_3+k28_3+k29_3);
f24_1=(2*k25_3-2*k22_3)- (2*k25_2-2*k22_2);
f24_2=2*(2*k25_1-2*k22_1);
f24_3=(2*k25_2-2*k22_2)+ (2*k25_3-2*k22_3);
f25_1=(4*k21_3)- (4*k21_2);
f25_2=2*(4*k21_1);
f25_3=(4*k21_2)+ (4*k21_3);
f26_1=(2*k22_3+2*k25_3)- (2*k22_2+2*k25_2);
f26_2=2*(2*k22_1+2*k25_1);
f26_3=(2*k22_2+2*k25_2)+ (2*k22_3+2*k25_3);
f27_1=(-k24_3+k26_3+k29_3-k28_3)- (-k24_2+k26_2+k29_2-k28_2);
f27_2=2*(-k24_1+k26_1+k29_1-k28_1);
f27_3=(-k24_2+k26_2+k29_2-k28_2)+(-k24_3+k26_3+k29_3-k28_3);
f28_1=(2*k23_3+2*k27_3)- (2*k23_2+2*k27_2);
f28_2=2*(2*k23_1+2*k27_1);
f28_3=(2*k23_2+2*k27_2)+ (2*k23_3+2*k27_3);
f29_1=(k24_3+k26_3+k28_3+k29_3)- (k24_2+k26_2+k28_2+k29_2);
f29_2=2*(k24_1+k26_1+k28_1+k29_1);
f29_3=(k24_2+k26_2+k28_2+k29_2)+ (k24_3+k26_3+k28_3+k29_3);
f31_1=(k34_3+k39_3-k36_3-k38_3)- (k34_2+k39_2-k36_2-k38_2);
f31_2=2*(k34_1+k39_1-k36_1-k38_1);
f31_3=(k34_2+k39_2-k36_2-k38_2)+(k34_3+k39_3-k36_3-k38_3);
f32_1=(-2*k33_3+2*k37_3)- (-2*k33_2+2*k37_2);
f32_2=2*(-2*k33_1+2*k37_1);
f32_3=(-2*k33_2+2*k37_2)+ (-2*k33_3+2*k37_3);
f33_1=(-k34_3-k36_3+k38_3+k39_3)- (-k34_2-k36_2+k38_2+k39_2);
f33_2=2*(-k34_1-k36_1+k38_1+k39_1);
f33_3=(-k34_2-k36_2+k38_2+k39_2)+ (-k34_3-k36_3+k38_3+k39_3);
f34_1=(2*k35_3-2*k32_3)- (2*k35_2-2*k32_2);
f34_2=2*(2*k35_1-2*k32_1);
f34_3=(2*k35_2-2*k32_2)+ (2*k35_3-2*k32_3);
f35_1=(4*k31_3)- (4*k31_2);
f35_2=2*(4*k31_1);
f35_3=(4*k31_2)+ (4*k31_3);
f36_1=(2*k32_3+2*k35_3)- (2*k32_2+2*k35_2);
f36_2=2*(2*k32_1+2*k35_1);
f36_3=(2*k32_2+2*k35_2)+ (2*k32_3+2*k35_3);
f37_1=(-k34_3+k36_3+k39_3-k38_3)- (-k34_2+k36_2+k39_2-k38_2);
f37_2=2*(-k34_1+k36_1+k39_1-k38_1);
f37_3=(-k34_2+k36_2+k39_2-k38_2)+ (-k34_3+k36_3+k39_3-k38_3);
f38_1=(2*k33_3+2*k37_3)- (2*k33_2+2*k37_2);
f38_2=2*(2*k33_1+2*k37_1);
f38_3=(2*k33_2+2*k37_2)+ (2*k33_3+2*k37_3);
f39_1=(k34_3+k36_3+k38_3+k39_3)- (k34_2+k36_2+k38_2+k39_2);
f39_2=2*(k34_1+k36_1+k38_1+k39_1);
f39_3=(k34_2+k36_2+k38_2+k39_2)+ (k34_3+k36_3+k38_3+k39_3);
f41_1=(k44_3+k49_3-k46_3-k48_3)- (k44_2+k49_2-k46_2-k48_2);
f41_2=2*(k44_1+k49_1-k46_1-k48_1);
f41_3=(k44_2+k49_2-k46_2-k48_2)+ (k44_3+k49_3-k46_3-k48_3);
f42_1=(-2*k43_3+2*k47_3)- (-2*k43_2+2*k47_2);
f42_2=2*(-2*k43_1+2*k47_1);
f42_3=(-2*k43_2+2*k47_2)+ (-2*k43_3+2*k47_3);
f43_1=(-k44_3-k46_3+k48_3+k49_3)- (-k44_2-k46_2+k48_2+k49_2);
f43_2=2*(-k44_1-k46_1+k48_1+k49_1);
f43_3=(-k44_2-k46_2+k48_2+k49_2)+ (-k44_3-k46_3+k48_3+k49_3);
f44_1=(2*k45_3-2*k42_3)- (2*k45_2-2*k42_2);
f44_2=2*(2*k45_1-2*k42_1);
f44_3=(2*k45_2-2*k42_2)+ (2*k45_3-2*k42_3);
f45_1=(4*k41_3)- (4*k41_2);
f45_2=2*(4*k41_1);
f45_3=(4*k41_2)+ (4*k41_3);
f46_1=(2*k42_3+2*k45_3)- (2*k42_2+2*k45_2);
f46_2=2*(2*k42_1+2*k45_1);
f46_3=(2*k42_2+2*k45_2)+ (2*k42_3+2*k45_3);
f47_1=(-k44_3+k46_3+k49_3-k48_3)- (-k44_2+k46_2+k49_2-k48_2);
f47_2=2*(-k44_1+k46_1+k49_1-k48_1);
f47_3=(-k44_2+k46_2+k49_2-k48_2)+ (-k44_3+k46_3+k49_3-k48_3);
f48_1=(2*k43_3+2*k47_3)- (2*k43_2+2*k47_2);
f48_2=2*(2*k43_1+2*k47_1);
f48_3=(2*k43_2+2*k47_2)+ (2*k43_3+2*k47_3);
f49_1=(k44_3+k46_3+k48_3+k49_3)- (k44_2+k46_2+k48_2+k49_2);
f49_2=2*(k44_1+k46_1+k48_1+k49_1);
f49_3=(k44_2+k46_2+k48_2+k49_2)+ (k44_3+k46_3+k48_3+k49_3);
f51_1=(k54_3+k59_3-k56_3-k58_3)- (k54_2+k59_2-k56_2-k58_2);
f51_2=2*(k54_1+k59_1-k56_1-k58_1);
f51_3=(k54_2+k59_2-k56_2-k58_2)+ (k54_3+k59_3-k56_3-k58_3);
f52_1=(-2*k53_3+2*k57_3)- (-2*k53_2+2*k57_2);
f52_2=2*(-2*k53_1+2*k57_1);
f52_3=(-2*k53_2+2*k57_2)+ (-2*k53_3+2*k57_3);
f53_1=(-k54_3-k56_3+k58_3+k59_3)- (-k54_2-k56_2+k58_2+k59_2);
f53_2=2*(-k54_1-k56_1+k58_1+k59_1);
f53_3=(-k54_2-k56_2+k58_2+k59_2)+ (-k54_3-k56_3+k58_3+k59_3);
f54_1=(2*k55_3-2*k52_3)- (2*k55_2-2*k52_2);
f54_2=2*(2*k55_1-2*k52_1);
f54_3=(2*k55_2-2*k52_2)+ (2*k55_3-2*k52_3);
f55_1=(4*k51_3)- (4*k51_2);
f55_2=2*(4*k51_1);
f55_3=(4*k51_2)+ (4*k51_3);
f56_1=(2*k52_3+2*k55_3)- (2*k52_2+2*k55_2);
f56_2=2*(2*k52_1+2*k55_1);
f56_3=(2*k52_2+2*k55_2)+ (2*k52_3+2*k55_3);
f57_1=(-k54_3+k56_3+k59_3-k58_3)- (-k54_2+k56_2+k59_2-k58_2);
f57_2=2*(-k54_1+k56_1+k59_1-k58_1);
f57_3=(-k54_2+k56_2+k59_2-k58_2)+ (-k54_3+k56_3+k59_3-k58_3);
f58_1=(2*k53_3+2*k57_3)- (2*k53_2+2*k57_2);
f58_2=2*(2*k53_1+2*k57_1);
f58_3=(2*k53_2+2*k57_2)+ (2*k53_3+2*k57_3);
f59_1=(k54_3+k56_3+k58_3+k59_3)- (k54_2+k56_2+k58_2+k59_2);
f59_2=2*(k54_1+k56_1+k58_1+k59_1);
f59_3=(k54_2+k56_2+k58_2+k59_2)+ (k54_3+k56_3+k58_3+k59_3);
f61_1=(k64_3+k69_3-k66_3-k68_3)- (k64_2+k69_2-k66_2-k68_2);
f61_2=2*(k64_1+k69_1-k66_1-k68_1);
f61_3=(k64_2+k69_2-k66_2-k68_2)+ (k64_3+k69_3-k66_3-k68_3);
f62_1=(-2*k63_3+2*k67_3)- (-2*k63_2+2*k67_2);
f62_2=2*(-2*k63_1+2*k67_1);
f62_3=(-2*k63_2+2*k67_2)+ (-2*k63_3+2*k67_3);
f63_1=(-k64_3-k66_3+k68_3+k69_3)- (-k64_2-k66_2+k68_2+k69_2);
f63_2=2*(-k64_1-k66_1+k68_1+k69_1);
f63_3=(-k64_2-k66_2+k68_2+k69_2)+ (-k64_3-k66_3+k68_3+k69_3);
f64_1=(2*k65_3-2*k62_3)- (2*k65_2-2*k62_2);
f64_2=2*(2*k65_1-2*k62_1);
f64_3=(2*k65_2-2*k62_2)+ (2*k65_3-2*k62_3);
f65_1=(4*k61_3)- (4*k61_2);
f65_2=2*(4*k61_1);
f65_3=(4*k61_2)+ (4*k61_3);
f66_1=(2*k62_3+2*k65_3)- (2*k62_2+2*k65_2);
f66_2=2*(2*k62_1+2*k65_1);
f66_3=(2*k62_2+2*k65_2)+ (2*k62_3+2*k65_3);
f67_1=(-k64_3+k66_3+k69_3-k68_3)- (-k64_2+k66_2+k69_2-k68_2);
f67_2=2*(-k64_1+k66_1+k69_1-k68_1);
f67_3=(-k64_2+k66_2+k69_2-k68_2)+ (-k64_3+k66_3+k69_3-k68_3);
f68_1=(2*k63_3+2*k67_3)- (2*k63_2+2*k67_2);
f68_2=2*(2*k63_1+2*k67_1);
f68_3=(2*k63_2+2*k67_2)+ (2*k63_3+2*k67_3);
f69_1=(k64_3+k66_3+k68_3+k69_3)- (k64_2+k66_2+k68_2+k69_2);
f69_2=2*(k64_1+k66_1+k68_1+k69_1);
f69_3=(k64_2+k66_2+k68_2+k69_2)+ (k64_3+k66_3+k68_3+k69_3);
%Constructing M matrix
A=[f11_1,f12_1,f13_1,f14_1,f15_1,f16_1,f17_1,f18_1,f19_1,0,0,0;f21_1,f22_1,f23_1,f24_1,f25_1,f26_1,f27_1,f28_1,f29_1,0,0,0;f31_1,f32_1,f33_1,f34_1,f35_1,f36_1,f37_1,f38_1,f39_1,0,0,0;f41_1,f42_1,f43_1,f44_1,f45_1,f46_1,f47_1,f48_1,f49_1,0,0,0;f51_1,f52_1,f53_1,f54_1,f55_1,f56_1,f57_1,f58_1,f59_1,0,0,0;f61_1,f62_1,f63_1,f64_1,f65_1,f66_1,f67_1,f68_1,f69_1,0,0,0;0,0,0,f11_1,f12_1,f13_1,f14_1,f15_1,f16_1,f17_1,f18_1,f19_1;0,0,0,f21_1,f22_1,f23_1,f24_1,f25_1,f26_1,f27_1,f28_1,f29_1;0,0,0,f31_1,f32_1,f33_1,f34_1,f35_1,f36_1,f37_1,f38_1,f39_1;0,0,0,f41_1,f42_1,f43_1,f44_1,f45_1,f46_1,f47_1,f48_1,f49_1;0,0,0,f51_1,f52_1,f53_1,f54_1,f55_1,f56_1,f57_1,f58_1,f59_1;0,0,0,f61_1,f62_1,f63_1,f64_1,f65_1,f66_1,f67_1,f68_1,f69_1];
B=[f11_2,f12_2,f13_2,f14_2,f15_2,f16_2,f17_2,f18_2,f19_2,0,0,0;f21_2,f22_2,f23_2,f24_2,f25_2,f26_2,f27_2,f28_2,f29_2,0,0,0;f31_2,f32_2,f33_2,f34_2,f35_2,f36_2,f37_2,f38_2,f39_2,0,0,0;f41_2,f42_2,f43_2,f44_2,f45_2,f46_2,f47_2,f48_2,f49_2,0,0,0;f51_2,f52_2,f53_2,f54_2,f55_2,f56_2,f57_2,f58_2,f59_2,0,0,0;f61_2,f62_2,f63_2,f64_2,f65_2,f66_2,f67_2,f68_2,f69_2,0,0,0;0,0,0,f11_2,f12_2,f13_2,f14_2,f15_2,f16_2,f17_2,f18_2,f19_2;0,0,0,f21_2,f22_2,f23_2,f24_2,f25_2,f26_2,f27_2,f28_2,f29_2;0,0,0,f31_2,f32_2,f33_2,f34_2,f35_2,f36_2,f37_2,f38_2,f39_2;0,0,0,f41_2,f42_2,f43_2,f44_2,f45_2,f46_2,f47_2,f48_2,f49_2;0,0,0,f51_2,f52_2,f53_2,f54_2,f55_2,f56_2,f57_2,f58_2,f59_2;0,0,0,f61_2,f62_2,f63_2,f64_2,f65_2,f66_2,f67_2,f68_2,f69_2];
C=[f11_3,f12_3,f13_3,f14_3,f15_3,f16_3,f17_3,f18_3,f19_3,0,0,0;f21_3,f22_3,f23_3,f24_3,f25_3,f26_3,f27_3,f28_3,f29_3,0,0,0;f31_3,f32_3,f33_3,f34_3,f35_3,f36_3,f37_3,f38_3,f39_3,0,0,0;f41_3,f42_3,f43_3,f44_3,f45_3,f46_3,f47_3,f48_3,f49_3,0,0,0;f51_3,f52_3,f53_3,f54_3,f55_3,f56_3,f57_3,f58_3,f59_3,0,0,0;f61_3,f62_3,f63_3,f64_3,f65_3,f66_3,f67_3,f68_3,f69_3,0,0,0;0,0,0,f11_3,f12_3,f13_3,f14_3,f15_3,f16_3,f17_3,f18_3,f19_3;0,0,0,f21_3,f22_3,f23_3,f24_3,f25_3,f26_3,f27_3,f28_3,f29_3;0,0,0,f31_3,f32_3,f33_3,f34_3,f35_3,f36_3,f37_3,f38_3,f39_3;0,0,0,f41_3,f42_3,f43_3,f44_3,f45_3,f46_3,f47_3,f48_3,f49_3;0,0,0,f51_3,f52_3,f53_3,f54_3,f55_3,f56_3,f57_3,f58_3,f59_3;0,0,0,f61_3,f62_3,f63_3,f64_3,f65_3,f66_3,f67_3,f68_3,f69_3];
O=zeros(12);
I=eye(12);
M=[O,I;-inv(A)*C,-inv(A)*B];
%Eigenvalues and eigenvectors of M matrices
[V,D]=eig(M);
D=diag(D);
for i=1:length(D)
     if isreal(D(i))==1
         %q3
         x3=D(i);
         q3= 2*atand(D(i));
         v=V(:,i);
         %q4,q5
         if abs(x3)<=1
            v=v(1:12);
            if max(abs(v))==abs(v(1))
                q4=2*atand(v(1)/v(4));
                q5=2*atand(v(1)/v(2));
            elseif max(abs(v))==abs(v(3))
                q4=2*atand(v(3)/v(6));
                q5=2*atand(v(2)/v(3));
            elseif max(abs(v))==abs(v(10))
                q4=2*atand(v(7)/v(10));
                q5=2*atand(v(10)/v(11));
            elseif max(abs(v))==abs(v(12))
                q4=2*atand(v(9)/v(12));
                q5=2*atand(v(11)/v(12));
            end
         else
            v=v(13:24);
            if max(abs(v))==abs(v(1))
                q4=2*atand(v(1)/v(4));
                q5=2*atand(v(1)/v(2));
            elseif max(abs(v))==abs(v(3))
                q4=2*atand(v(3)/v(6));
                q5=2*atand(v(2)/v(3));
            elseif max(abs(v))==abs(v(10))
                q4=2*atand(v(7)/v(10));
                q5=2*atand(v(10)/v(11));
            elseif max(abs(v))==abs(v(12))
                q4=2*atand(v(9)/v(12));
                q5=2*atand(v(11)/v(12));
            end
         end
         %q1
         c3=cosd(q3);
         s3=sind(q3);
         c4=cosd(q4);
         s4=sind(q4);
         c5=cosd(q5);
         s5=sind(q5);
         m1= s5*u5;
         m2= c5*y4*u5+u4*y5;
         m3= -c5*u4*u5+y4*y5;
         g1= c5*a5+a4;
         g2= -s5*y4*a5+u4*d5;
         g3= s5*u4*a5+y4*d5+d4;
         r1= c4*m1+s4*m2;
         r2= -y3*(s4*m1-c4*m2)+u3*m3;
         r3= u3*(s4*m1-c4*m2)+y3*m3;
         f1= c4*g1+s4*g2+a3;
         f2= -y3*(s4*g1-c4*g2)+u3*g3;
         f3= u3*(s4*g1-c4*g2)+y3*g3+d3;
         v= my*u6+ny*y6;
         s1 =(u1*v*(d2 + f3*y2 + y1*(d1 - y) - u2*(c3*f2 - f1*s3)) + q*u1*(w*y1 - r3*y2 + u2*(c3*r2 - r1*s3)))/z;
         q1=asind(s1);
         re1=[];
         re1(length(re1)+1)=q1;
         re1(length(re1)+1)=180-q1;
         re6=[];
         for i=1:length(re1)
            q1=re1(i);
            %q6
            c1=cosd(q1);
            s1=sind(q1);
            n= u6*(nz*y1-c1*ny*u1+nx*s1*u1)- y6*(mz*y1-c1*my*u1+mx*s1*u1);
            m= lz*y1-c1*ly*u1+lx*s1*u1;
            l= s5*(u4*(y2*y3 - c3*u2*u3) + c4*y4*(u3*y2 + c3*u2*y3) - s3*s4*u2*y4) + c5*(s4*(u3*y2 + c3*u2*y3) + c4*s3*u2);
            x61=(n-sqrt(m*m+n*n-l*l))/(l+m);
            x62=(n+sqrt(m*m+n*n-l*l))/(l+m);
            q61=2*atand(x61);
            q62=2*atand(x62);
            re6(length(re6)+1)=q61;
            re6(length(re6)+1)=q62;
            for i=1:length(re6)
                q6=re6(i);
                if isreal(q6)==1
                     %q2
                     c6=cosd(q6);
                     s6=sind(q6);
                     n= ly*s1*u2+lx*c1*u2;
                     m= lx*s1*u2*y1-lz*u1*u2-ly*c1*u2*y1;
                     l= s6*(u5*(y3*y4-c4*u3*u4)+c5*y5*(u4*y3+c4*u3*y4)-s4*s5*u3*y5)+c6*(s5*(u4*y3+c4*u3*y4)+c5*s4*u3)-lz*y1*y2+ ly*c1*u1*y2-lx*s1*u1*y2;
                     x21=(n+sqrt(m*m+n*n-l*l))/(l+m);
                     q21=2*atand(x21);
                     if isreal(q21)==1
                        ret=vpa([q1,q21,q3,q4,q5,q6],7);
                     end
                     x22=(n-sqrt(m*m+n*n-l*l))/(l+m);                
                     q22=2*atand(x22);
                     if isreal(q21)==1
                        ret=vpa([q1,q22,q3,q4,q5,q6],7)
                     end
                end
           end
         end
     end
end




