function  [VF, y] =  FW_III_ord_VF_F(VF, x)

% This function evaluate the output of a III-order Volterra nonlinear 
% structure. 
%
% Details can be found in:
% 
%
%
% USAGE:
%   [VF, y(n)] = FW_III_ord_VF_F(VF, x(n))
%
% Input Arguments:
%   VF      : ADAPTIVE VOLTERRA FILTER STRUCT
%   x       : x[n] input signal sample
%
% Output Arguments:
%   VF      : ADAPTIVE VOLTERRA FILTER STRUCT
%   y       : y[n] output signal sample
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2012/01/27$
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


M = VF.M;    % Number of coefficients

% Regression delay-line update (I ord Volterra Kernel) --------------------
VF.xh1(2 : M) = VF.xh1(1 : M-1);    % Shift of the input delay-line
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
VF.xh = [VF.xh1; VF.xh2; VF.xh3 ];     % Full buffer composition
    
y = VF.h.'*VF.xh;   % (K,1) = (K,M) x (M,1) output array    
% -------------------------------------------------------------------------  

% -------------------------------------------------------------------------
% © 2012
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------