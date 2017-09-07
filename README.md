This code is for the paper: Deng, Hong, et al. "A Coordinate Descent Method for Total Variation Minimization."

If any problem, please contact with Hong Deng (denghong_hit@163.com ) or Wangmeng Zuo (cswmzuo@gmail.com)

1. Denoising

Please run "test_CoD.m". 

Our core code is implemented in C/C++, and we provide the executable file compiled under Windows 64bit. For ohter platforms, please run "compile_CoD.m".

For isotropic TV, please use "CoorDenoiseLC_cyc" with L = 3. 
For anisotropic TV, we recommend to use "CoorDenoiseC_med" which is faster than "CoorDenoiseLC_cyc" with L = 1. 

2. Deblurring

Please run "test_deblur.m".
