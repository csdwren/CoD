function [ x,iter] = ALMCoD(y,H,lambda,L,eps)
%  y observed image
%  iter iterative times
%  H blur module
[m,n]=size(y);
% tic;
x=y;
p=zeros(m,n);

%%
H_FFT=psf2otf(H,[m,n]);
HC_FFT = conj(H_FFT);

AsDy=HC_FFT.*fft2(y);
ATA_FFT=H_FFT.*HC_FFT;
%%
delta=1e-4;

%  lef=abs(norm(x(:)-x_pre(:))/norm(x(:)));
%      func(1)=lef;
% delta_max=1000;

maxIter=100;
perturb=20e-3;
%%
for i=1:maxIter
    if L==0
        %         u=tvdenoise(x-p,delta/lambda,5);
        pars.MAXITER=5;
        pars.print=0;
        pars.tv='iso';
        u=denoise_bound(x-p,lambda/delta,-inf,inf,pars);
    elseif L==10
        u=CoorDenoiseC_v2(x-p,lambda/delta,5);
    elseif L==1
        %         if i==3
        u=CoorDenoiseLC_cyc(x-p,lambda/delta,3,L,10e-3);
        %         perturb=perturb/1.2;
        %         else
        %              u=CoorDenoiseLC_cyc(x-p,lambda/delta,10,L,0e-3,x);
        %         end
    elseif L==3
        u=CoorDenoiseLC_cyc(x-p,lambda/delta,3,L,1e-3);
    end
    %   [u,pux,puy] = chambolle_prox_TV_stop(real(x-bu), 'lambda', threshold, 'maxiter', 10, 'dualvars',[pux puy]);
    x_pre=x;
    x=real(ifft2((AsDy+delta*fft2(u+p))./(ATA_FFT+delta)));
    
    normdx=norm(x_pre(:)-x(:));
    if  normdx/norm(x(:))<=eps;
        %          fprintf('iter: %d\n',i);
        iter=i;
        break;
    end
    %%
    p=p+u-x;
    
    if (normdx)*delta/norm(u(:))<1e-3
        rao=1.9;
    else
        rao=1;
    end
    delta=rao*delta;
    
end

if i==maxIter
    iter=i;
end


