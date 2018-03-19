function [t,lambda,c_m_n_mdf]=setup_filter(filter,N,level,a)
    if level>1
        [phi,~,~,~,~] = wavefun(filter,level);
    else
        [phi,~,~,~]=wfilters(filter);
    end
    tmp1=find(phi~=0,1);
    tmp2=find(phi~=0,1,'last');
    phi=phi(tmp1:tmp2);
    shift=ceil((length(phi)-1)/2/(2^level));
    shift=0;
    phi_len=shift*2*2^level+1;
    m=round((phi_len-length(phi))/2);
    phi=[zeros(1,m) phi zeros(1,m)];
    shift_test=((length(phi)-1)/2/(2^level));
    dwtmode('per','nodisplay');
    P=N-1;    
    T=1/N;
    T_s = T/(2^level);

    if mod(N, 2) == 0
        n1 = -N/2;
        n2 = N/2 - 1;
    else
        n1 = -(N-1)/2;
        n2 = (N-1)/2;
    end
    n_vec = (n1:n2)';
    t1    =  n1 * T;
    t2    = (n2+1) * T - T_s;
    t     = (t1:T_s:t2)';

    m = 0:P;

    L_den = a*(P+1);
    alpha_0 = -1j * pi * P / L_den ;
    lambda  = 2* 1j * pi / L_den ;
    alpha_vec = alpha_0 + lambda * m;
    [nphi, t_nphi] = generate_b_spline(0, T_s, T);
    newphi = (T_s) * conv(phi, nphi);%

    [c_m_0]=get_c0(newphi(:)',T_s,T,alpha_vec);
    [X,Y]=meshgrid(n_vec,(alpha_vec));
    c_m_n_mdf = diag(c_m_0)*exp(X.*Y);
    c_m_n_mdf = circshift(c_m_n_mdf.',shift).';
end