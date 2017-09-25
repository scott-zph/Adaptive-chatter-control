function [m,n,k,l1] = mnkl1(M,C,K,F,u,omega,t,l,x1,x2,y1,y2)
m = y1*g1(t,omega,l,x1,x2);
n = y2*g1(t,omega,l,x1,x2);
k = [y1;y2];
l1 = -M^(-1)*C*[y1;y2] -M^(-1)*K*[x1;x2] + M^(-1)*F + M^(-1)*u;