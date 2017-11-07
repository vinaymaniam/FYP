function [c0]=get_c0(phi,T_s,T,alpha_vec)
    T_s = T_s / T;
    t_kernel=(0:1:length(phi)-1)';
    Phi0=[];
    phi=phi;
    
    for counter=1:length(alpha_vec)
        Phi0=[Phi0;phi*exp(-T_s*alpha_vec(counter)*t_kernel)]; %laplace transform of phi
    end;
    c0=1./Phi0;