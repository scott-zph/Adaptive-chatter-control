function g2 = g2(t,omega,l,x1,x2,x3,x4)
h = 1;
    for i = 1:l
        h = [h;cos(i*omega*t);sin(i*omega*t)];
    end
g2 = [(x1-x2)*h;(x3-x4)*h];
    