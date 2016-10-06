function Volterra_3rd_f = create_III_ord_Volterra_filter_1(M, K, mu, delta)

% This function creates and initializes the structure implementing a III-order 
% Volterra AF architecture and adapted by the LMS adaptive algorithm.
%
% Details can be found in:
% 
%
%
% USAGE:
% Volterra_3rd_f = create_III_ord_Volterra_filter_1(M, K, mu, delta)
%   - M: number of filter coefficients
%   - K: APA order
%   - mu: step-size of the lineaar filter
%   - delta: regularization parameter
%   - Volterra_3rd_f: structure of the created III-order Volterra filter
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


M2 = M*(M+1)/2;          % Copmuted II  degree Volterra filter length            
M3 = M*(M+2)*(M+1)/6;    % Computed III degree Volterra filter length            
MP = M + M2 + M3;        % Full Volterra filter length   
dI  = delta*eye(K);      % Regularizing diagonal matrix
bl  = zeros(M+1,2);      % Upper block pointers address  
xd  = zeros(K,1);        % Buffer for desired samples
xh1 = zeros(M,1);        % Buffer for I   degree Volterra filter status
xh2 = zeros(M2,1);       % Buffer for II  degree Volterra filter status
xh3 = zeros(M3,1);       % Buffer for III degree Volterra filter status
xh  = [xh1;  xh2; xh3];  % Volterra filter status  
h   = zeros(MP,1);       % Volterra AF taps  
X  = zeros(K,MP);        % Total covariance data matrix for Volterra APA
X1 = zeros(K,M);         % I   ord cov. data matrix for Volterra APA
X2 = zeros(K,M2);        % II  ord cov. data matrix for Volterra APA
X3 = zeros(K,M3);        % III ord cov. data matrix for Volterra APA
h1 = zeros(M,1);         % Buffer for I   degree Volterra filter
h2 = zeros(M2,1);        % Buffer for II  degree Volterra filter
h3 = zeros(M3,1);        % Buffer for III degree Volterra filter


% II order vector-block of pointers to the Kronecker product --------------
for k=1:M 
    bl(k,1)  = k ;   % Upper block addres (1 : bl(k)) 
end;

% III order vector-block of pointers to the Kronecker product -------------
kk = 1;
for i = 1:M
    kk = kk + bl(M,1) - bl(i,1) +  1 ;
    bl(i+1,2) = kk - 1;  % Dummy pointer array for next order block addressing
end
bl(1,2) = 0;

Circuit = 'III order full Volterra FIlter';   % Architecture name

% -------------------------------------------------------------------------
Volterra_3rd_f = struct('Circuit', Circuit, ... 
    'M',M,'M2',M2,'M3',M3,'MP',MP, ...                    % Memory lengths
    'xh1',xh1,'xh2',xh2,'xh3',xh3,'xh',xh, ...            % filters status
    'bl',bl, ...                                          % pointers matrix for the Kronecker product    
    'K',K,'mu',mu,'dI',dI, ...                            % APA order, learning rate regularizing matrix  
    'h1',h1,'h2',h2,'h3',h3,'h', h, ...                   % array of Volterra coeffs.
    'X', X, 'X1', X1, 'X2', X2, 'X3', X3, 'xd',xd);       % APA matricies and statusd status
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2012
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------