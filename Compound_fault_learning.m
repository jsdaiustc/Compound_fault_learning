function Y_res=Compound_fault_learning(Y, Fs, Fb)

%%  normalize Y
N_orginal=size(Y,1);
Norm_y=norm(Y,'fro')/sqrt(N_orginal);
Y=Y/Norm_y;

K=length(Fb);
R_all=Fs./Fb;
F_num=K;
L_all= floor(N_orginal./R_all);
R_all=floor(R_all);

index_record_all=zeros(N_orginal,K);
for kk=1:K
    Ini_pints=[0:L_all(kk)-1].*R_all(kk)+1;
    Ini_pints=round(Ini_pints);
    R_all(kk)=round(R_all(kk));
    index_record=[];
    for ii=1:length(Ini_pints)
        index_record=[index_record; [Ini_pints(ii):Ini_pints(ii)+  R_all(kk)-1]'  ];
    end
    record_num(kk)=length(index_record);
    index_record_all(1:record_num(kk),kk)=index_record;
end



%%  Initialization
a=1e-20;
b=a;
N=N_orginal;
gamma_with_tau=ones(round(N*1.1),K);
G=eye(K);

tau=ones(max(L_all),K);
gamma_core=zeros(max(R_all),K);
Z=ones(round(N*1.1),3,K)/3;
alpha=1;
iter = 0;
maxiter=100;
rho=0.5;
converged=false;

mu=zeros(round(N),F_num);
Sigma=zeros(round(N),F_num);
Y_res=zeros(N_orginal,F_num);

while ~converged
    %%  Update x
    for ff=1:F_num
        y_new=  Y- sum(Y_res,2) +  Y_res(:,ff);
        gamma_with_tau_temp=[];
        for kk=1:K
            index_active=1 :record_num(kk);
            gamma_with_tau_kk(index_record_all(index_active,kk),kk)=gamma_with_tau( index_active ,kk);
            gamma_with_tau_kk=gamma_with_tau_kk(:,kk);
            deltal=[gamma_with_tau_kk(2:length(gamma_with_tau_kk) );gamma_with_tau_kk(1)];
            deltar=[gamma_with_tau_kk(end);gamma_with_tau_kk(1:length(gamma_with_tau_kk)-1)];
            gamma_with_tau_temp(1:length(gamma_with_tau_kk),kk)=G(ff,kk)*( Z(1:length(gamma_with_tau_kk),1,kk).*deltar ...
                + Z(1:length(gamma_with_tau_kk),2,kk).*gamma_with_tau_kk +  Z(1:length(gamma_with_tau_kk),3,kk).*deltal );
        end
        sum_gamma_with_tau=sum(gamma_with_tau_temp,2);
        Sigma(1:length(sum_gamma_with_tau),ff)=1./(alpha + sum_gamma_with_tau );
        mu(1:length(sum_gamma_with_tau),ff)= alpha*( Sigma(1:length(sum_gamma_with_tau),ff).* y_new(1:length(sum_gamma_with_tau) ) );
        Y_res(1:length(sum_gamma_with_tau),ff)=mu(1:length(sum_gamma_with_tau),ff);
    end
    
    %% Update alpha
    alpha_old=alpha;
    resid=Y-sum(Y_res,2);
    alpha=( 1*N + 2*a )/( 2*b +  norm(resid, 'fro')^2  + 1*  sum(Sigma(:))  );
    alpha_rho=0.5;
    alpha= alpha*(1-alpha_rho)  + alpha_old*alpha_rho; %% damping
    
    
    %% Update gamma
    for kk=1:K
        L_kk= L_all(kk);
        R_kk= R_all(kk);
        index_active=1 :record_num(kk);
        Z_combine=[];mu2_combine=[];
        for ff=1:F_num
            sum_mu=sum( mu(:,ff).*conj(mu(:,ff)), 2);
            mu2=sum_mu + real(Sigma(:,ff));
            mu2l=[mu2(2:length(mu2));mu2(1)];
            mu2r=[mu2(end);mu2(1:length(mu2)-1)];
            Z_2=Z(1:length(mu2),2,kk);
            Z_1l=[Z(2:length(mu2),1,kk);Z(1,1,kk)];
            Z_3r=[Z(end,3,kk);Z(1:length(mu2)-1,3,kk)];
            Z_combine(1:length(mu2),ff)=G(ff,kk)*(Z_1l + Z_2 + Z_3r);
            mu2_combine(1:length(mu2),ff)=G(ff,kk)*(Z_1l.*mu2l + Z_2.*mu2 + Z_3r.*mu2r);
        end
        sum_Z=sum(Z_combine,2);
        sum_Z=sum_Z( index_record_all(index_active,kk ) );
        sum_mu2=sum(mu2_combine,2);
        sum_mu2=sum_mu2( index_record_all(index_active,kk) );
        
        sum_Z_reshape=reshape(sum_Z,R_kk,L_kk);
        sum_mu2_reshape=reshape(sum_mu2,R_kk,L_kk);
        c_k=sum(sum_Z_reshape,2)+2*a;
        d_k=sum(sum_mu2_reshape*diag(tau(1:L_kk,kk)) , 2)+2*b;
        gamma_core(1:R_kk,kk)=c_k./d_k;
        gamma_with_tau(1:L_kk*R_kk,kk)=kron( tau(1:L_kk,kk) , gamma_core(1:R_kk,kk) );
        ln_gamma_core= psi( c_k ) -  log(  d_k  );
        ln_delta_with_ones(1:L_kk*R_kk,kk)=kron(ones( L_kk,1),ln_gamma_core ) ;
        
    end
    
    
    if iter>10
        %% update F
        var_y=var(Y_res);
        position=find(var_y/max(var_y)<1/50);
        if ~isempty(position)
            Y_res(:,position)=[];
            Sigma(:,position)=[];
            mu(:,position)=[];
            G(position,:)=[];
            G(:,position)=1e-10;
            F_num=F_num-length(position);
        end
        %%  update G
        eg=zeros(F_num,K);
        for ff=1:F_num
            for kk=1:K
                L_kk= L_all(kk);
                R_kk= R_all(kk);
                index_active=1 :record_num(kk);
                
                ln_delta_with_ones_kg(index_record_all(index_active,kk),kk)=ln_delta_with_ones( index_active ,kk);
                ln_delta_with_ones_kg=ln_delta_with_ones_kg(:,kk);
                ln_deltal=[ln_delta_with_ones_kg(2:end);ln_delta_with_ones_kg(1)];
                ln_deltar=[ln_delta_with_ones_kg(end);ln_delta_with_ones_kg(1:end-1)];
                
                gamma_with_tau_k3(index_record_all(index_active,kk),kk)=gamma_with_tau( index_active ,kk);
                gamma_with_tau_k3=gamma_with_tau_k3(:,kk);
                deltal=[gamma_with_tau_k3(2:length(gamma_with_tau_k3));gamma_with_tau_k3(1)];
                deltar=[gamma_with_tau_k3(end);gamma_with_tau_k3(1:length(gamma_with_tau_k3)-1)];
                
                sum_mu=sum( mu(:,ff).*conj(mu(:,ff)), 2);
                mu2=sum_mu + real(Sigma(:,ff));
                tt1=sum(Z(1:length(gamma_with_tau_k3),1,kk).*ln_deltar+ ...
                    Z(1:length(gamma_with_tau_k3),2,kk).*ln_delta_with_ones_kg+ Z(1:length(gamma_with_tau_k3),3,kk).*ln_deltal)  ;
                tt1=tt1   + sum(  kron( log(tau(1:L_kk,kk) ), ones(R_kk,1) ) );
                tt2=sum( Z(1:length(gamma_with_tau_k3),1,kk).*deltar.*mu2(1:length(gamma_with_tau_k3)) + ...
                    Z(1:length(gamma_with_tau_k3),2,kk).*gamma_with_tau_k3.*mu2(1:length(gamma_with_tau_k3)) ...
                    +  Z(1:length(gamma_with_tau_k3),3,kk).*deltal.*mu2(1:length(gamma_with_tau_k3)) );
                eg(ff,kk)=(tt1-tt2)/2 ;
            end
        end
        G_old=G;
        max_temp=max(eg')';
        eg=eg-max_temp*ones(1,K);
        temp1=exp(eg);
        G=(diag(  1./sum(temp1,2) ) *   temp1);
        G= (1-0.9)*G + 0.9*G_old; %% damping
        
        
    end
    
    %% Update Z
    for kk=1:K
        index_active=1 :record_num(kk);
        ln_delta_with_ones_kk(index_record_all(index_active,kk),kk)=ln_delta_with_ones( index_active ,kk);
        ln_delta_with_ones_kk=ln_delta_with_ones_kk(:,kk);
        ln_deltal=[ln_delta_with_ones_kk(2:end);ln_delta_with_ones_kk(1)];
        ln_deltar=[ln_delta_with_ones_kk(end);ln_delta_with_ones_kk(1:end-1)];
        
        gamma_with_tau_kk1(index_record_all(index_active,kk),kk)=gamma_with_tau( index_active ,kk);
        gamma_with_tau_kk1=gamma_with_tau_kk1(:,kk);
        deltal=[gamma_with_tau_kk1(2:length(gamma_with_tau_kk1));gamma_with_tau_kk1(1)];
        deltar=[gamma_with_tau_kk1(end);gamma_with_tau_kk1(1:length(gamma_with_tau_kk1)-1)];
        temp_t1=[];temp_t2=[];temp_t3=[];
        for ff=1:F_num
            sum_mu=sum( mu(:,ff).*conj(mu(:,ff)), 2);
            mu2=sum_mu + real(Sigma(:,ff));
            temp_t1(:,ff) =G(ff,kk)*(ln_deltar -mu2(1:length(ln_delta_with_ones_kk)).*deltar) ;
            temp_t2(:,ff) =G(ff,kk)*(ln_delta_with_ones_kk - mu2(1:length(ln_delta_with_ones_kk)).*gamma_with_tau_kk1);
            temp_t3 (:,ff)=G(ff,kk)*(ln_deltal - mu2(1:length(ln_delta_with_ones_kk)).*deltal) ;
        end
        t1=sum(temp_t1,2);t2=sum(temp_t2,2);t3=sum(temp_t3,2);
        et=[t1,t2,t3]/2;
        max_et=max(et')';
        et=et-max_et*ones(1,3);
        temp_p= exp(et);
        Z_kk_old=Z(1:length(temp_p),:,kk);
        Z_kk= ((  1./sum(temp_p,2) )*ones(1,3)).*  temp_p;
        Z(1:length(temp_p),:,kk)=(1-rho)*Z_kk + rho*Z_kk_old;  % damping
    end
    
    %% Update tau
    for kk=1:K
        L_kk= L_all(kk);
        R_kk= R_all(kk);
        index_active=1 :record_num(kk);
        
        if iter>=20
            [~,ind_sort]=sort(gamma_core(1:R_kk,kk),'descend');
            gamma_core(ind_sort(1:round(end*0.8)),kk)=1e3;
            if iter>=maxiter-10
                gamma_core(ind_sort(1:round(end*0.8)),kk)=1e5;
            end
        end
        
        delta_11=kron( ones(L_kk,1),gamma_core(1:R_kk,kk) );
        delta_11_kk(index_record_all(index_active,kk),kk )=delta_11( index_active);
        delta_11_kk=delta_11_kk(:,kk);
        deltal=[delta_11_kk(2:length(delta_11_kk));delta_11_kk(1)];
        deltar=[delta_11_kk(end);delta_11_kk(1:length(delta_11_kk)-1)];
        delta_temp=[];
        for ff=1:F_num
            sum_mu=sum( mu(:,ff).*conj(mu(:,ff)), 2);
            mu2=sum_mu + real(Sigma(:,ff));
            delta_temp(:,ff)=G(ff,kk)*( Z(1:length(delta_11_kk),1,kk).*deltar ...
                + Z(1:length(delta_11_kk),2,kk).*delta_11_kk + Z(1:length(delta_11_kk),3,kk).*deltal ).*mu2(1:length(delta_11_kk));
        end
        sum_delta_temp=sum(delta_temp,2);
        sum_delta_temp=sum_delta_temp( index_record_all(index_active,kk ) );
        temp2=reshape(sum_delta_temp, R_kk, L_kk );
        cl= sum(temp2,1)';
        tau_old= tau(1:L_kk,kk);
        tau(1:L_kk,kk)= ( R_kk*sum(G(:,kk)) )./( cl );
        tau(1:L_kk,kk)=rho* tau(1:L_kk,kk)+ (1-rho)*tau_old; % damping
        gamma_with_tau(1:L_kk*R_kk,kk)=kron( tau(1:L_kk,kk),gamma_core(1:R_kk,kk) );
        
    end
    
    if iter >= maxiter
        converged = true;
    end
    iter = iter + 1;
    
end
Y_res=Y_res*Norm_y;



