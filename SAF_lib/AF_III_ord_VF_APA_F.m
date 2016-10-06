function  [VF, y, e] = AF_III_ord_VF_APA_F(VF, x, d)

% This function implements the adaptation of a III-order Volterra AF  
% structure by using the APA adaptive algorithm.
%
% Details can be found in:
% 
%
%
% USAGE:
%  [VF, y, e] = AF_III_ord_VF_APA_F(VF, x, d)
%
% Input Arguments:
%   F       : Volterra AF structurs
%   x       : array of input signal
%   d       : array of desired signal
%  
% Output Arguments:
%   F       : Volterra AF structure
%   y       : array of output signal   
%   e       : array of error signal    
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2012/01/22$
% $Revision: 1.1$  $Date: 2016/07/31$
% -------------------------------------------------------------------------
% License to use and modify this code is granted freely without warranty to
% all, as long as the original authors are referenced and attributed as such.
% The original authors maintain the right to be solely associated with this
% work.
% -------------------------------------------------------------------------
% © 2012
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------


K = VF.K;    % APA order
M = VF.M;    % Number of coefficients

% Regression delay-line update (I ord Volterra Kernel) --------------------
VF.xh1(2 : M) = VF.xh1(1 : M-1);    % Shift the input delay-line
VF.xh1(1) = x;                      % Load a new input into the delay-line  
    
% Kronecker product for II order Volterra kernel generation ---------------
kk = 1;
MM = VF.bl(M,1); 
for i = 1:M
    VF.xh2(kk : kk + MM - VF.bl(i,1) ) = VF.xh1(i)*VF.xh1( VF.bl(i,1):MM );
    kk = kk + MM - VF.bl(i,1) +  1  ;
end
 
% Kronecker product for III order Volterra kernel generation --------------
kk = 1;
MM = VF.bl(M+1,2); 
for i = 1:M
    VF.xh3(kk : kk + MM - VF.bl(i,2) - 1 ) = VF.xh1(i)*VF.xh2( VF.bl(i,2)+1:MM );
    kk = kk + MM - VF.bl(i,2);
end

% Full buffer composition -------------------------------------------------
VF.xh = [VF.xh1; VF.xh2; VF.xh3 ];

% Update X covariance matrix for APA agorithm -----------------------------
M = VF.MP;       % Now M is the Volterra entire buffer length 
for i = 1:K-1    % Shift the covarince matrix X
    VF.X(i,1:M) = VF.X(i+1,1:M); 
end
VF.X(K,1:M) = VF.xh(1:M);  % Fill the covariance matrix X

% Desired-output buffer for APA -------------------------------------------
VF.xd(1:K-1) = VF.xd(2:K);  % Shift the input delay-line  
VF.xd(K) = d ;              % Load a new desired output into the delay-line      
yy = VF.X*VF.h;             % (K,1) = (K,M) x (M,1) output array    
ee = VF.xd - yy;            % (K,1) = (K,1) - (K,1) error array     
VF.h = VF.h + VF.mu*(VF.X'/(VF.dI + VF.X*VF.X'))*ee;  % (M,1) = (M,1) + (M,K) x (K,1)
y = yy(VF.K);               % Output sample
e = ee(VF.K);               % Error sample
% -------------------------------------------------------------------------  

% -------------------------------------------------------------------------
% © 2012
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------