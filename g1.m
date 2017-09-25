function g1 = g1(t,omega,l,x1,x2)
h = 1;
    for i = 1:l
        h = [h;cos(i*omega*t);sin(i*omega*t)];
    end
g1 = [x1*h;x2*h];    
    