

x = imread('imgs/1.tif');
if(length(size(x))==3)
    x=im2double(rgb2gray(x));
else
    x=im2double(x);
end


sigma=0.1;
lambda=0.3;

[m, n] = size(x);

%%


randn('seed',314);
y = x + sigma*randn(m,n);


%%

    result=CoorDenoiseLC_cyc(y,lambda,10,3,1e-3);
%         result=CoorDenoiseC_med(y,lambda,45,1e-1);
   
%%
figure,
imshow(result);

