function X = Hl(t,mega_omega,N,phi_e,phi_a,Kt,Kn)
    Hxx = 0;
    Hxy = 0;
    Hyx = 0;
    Hyy = 0;

    for j = 1:N
        phi_j = mod(2 * pi * mega_omega * t / 60 + 2 * pi * (j - 1) / N,2*pi);   
        s_phi_j = 1.*((phi_j > phi_e)&(phi_j < phi_a));
        Hxx = Hxx + s_phi_j * (Kt * cos(phi_j) + Kn * sin(phi_j)) * sin(phi_j);
        Hxy = Hxy + s_phi_j * (Kt * cos(phi_j) + Kn * sin(phi_j)) * cos(phi_j);
        Hyx = Hyx + s_phi_j * (-Kt * sin(phi_j) + Kn * cos(phi_j)) * sin(phi_j);
        Hyy = Hyy + s_phi_j * (-Kt * sin(phi_j) + Kn * cos(phi_j)) * cos(phi_j);
    end  
    X =[Hxx,Hxy;Hyx,Hyy];