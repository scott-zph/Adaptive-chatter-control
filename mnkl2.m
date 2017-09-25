function [m,n,k,l2] = mnkl2(M,C,K,F,u,omega,t,l,x1,x2,x3,x4,y1,y2)
m = y1*g2(t,omega,l,x1,x2,x3,x4);
n = y2*g2(t,omega,l,x1,x2,x3,x4);
k = [y1;y2];
l2 = -M^(-1)*C*[y1;y2] -M^(-1)*K*[x2;x4] + M^(-1)*F + M^(-1)*u;