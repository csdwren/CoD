

addpath(genpath('.'));

clear;

%%

    x = imread(['imgs/1.tif']);
    if(length(size(x))==3)
        x=im2double(rgb2gray(x));
    else
        x=im2double(x);
    end
    
 
    
    sigma=0.001;
    miu=1e-3;
    lambda=5e-5;
    
    
    
    [m, n] = size(x);
    
    %%
    %%get the oberverd image
    d=15;
    
    H = fspecial('gaussian',[d,d],7);%the gaussian kernel


    
    %get the convolution matrix in frequency domain
    H_FFT=psf2otf(H,[m,n]);
    HC_FFT = conj(H_FFT);
    
    %curl operator
    A = @(x) real(ifft2(H_FFT.*fft2(x)));
    AT = @(x) real(ifft2(HC_FFT.*fft2(x)));
    
    % y = A(x) + sigma*randn(m,n);% obeserved image
    y=imfilter(x,H,'circular','conv')+ sigma*randn(m,n);
    % y=imnoise(y,'salt & pepper',0.1);
    %%

    [result,iter]=ALMCoD(y,H,lambda,3,1e-4);
    
    imshow(result)

    
    
    
