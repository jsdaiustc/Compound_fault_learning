function [x_RTEA] = RTEA( y,B,lam1)



K0=4;
L=2;
for i = 1 : L
    x(:, i) = y;
end


beta1=0.325;
beta=0.325;
eta=0.1;
% lam_0=eta*beta;
    lam_0=0.1;
for i=1:2
    lam(i)=0.5*(1-eta)*beta1;
   lam(i)=lam1;
end

converged=false;
iter = 0;
maxiter=50;

while ~converged   
         
    temp=ones(1,K0);
    r_0=sqrt(conv(abs( sum(x,2)).^2, temp ));
    rr_0 = conv(1./(r_0 + 1e-10), temp , 'valid');
%       rr_0=0;
    
                         
    for i = 1 : L
        r = sqrt(conv(abs(x(:, i)).^2, B{1, i}));
        rr(:,i) = conv(1./(r + 1e-10), B{1, i} , 'valid');
        p(:,i)=2+2*lam_0*rr_0 + lam(i)*rr(:,i);
    end
    
    q1=y+ ( 1+lam_0*rr_0 ).*( x(:,1)-x(:,2) );
    q2=y+ ( 1+lam_0*rr_0 ).*( x(:,2)-x(:,1) );
    x(:,1)=q1./p(:,1);
    x(:,2)=q2./p(:,2);

           
    if iter >= maxiter
        converged = true;
    end
    iter=iter+1;

end

x_RTEA=x;



