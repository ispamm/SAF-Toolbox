% WSAF_IIR_Volterra_compare_demo.m
% 
% DEMO on comparisons of different nonlinear system identification architectures
% based on an IIR Wiener Spline Adaptive Filter (IIR_WSAF), WSAF, a IIR 
% polynomial filter and a III-order Volterra filter.
%
% This demo file implements the second experimental test shown in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear System
% Identification using IIR Spline Adaptive Filters", Signal Processing, 
% Vol. 108, pp. 30-35, March 2015.
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2014/04/14$
% $Revision: 1.1$  $Date: 2016/07/31$
% -------------------------------------------------------------------------
% License to use and modify this code is granted freely without warranty to
% all, as long as the original authors are referenced and attributed as such.
% The original authors maintain the right to be solely associated with this
% work.
% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------

clear all; 
close all;
disp('WSAF_IIR_Volterra_compare_demo');
% -------------------------------------------------------------------------

%%  Parameters setting

% Input parameters --------------------------------------------------------
Lx = 50000;                  % Length of input signal
nRun = 10;                   % Number of runs
out_noise_level_dB = Inf;    % SNR
out_noise_level = 10^(-out_noise_level_dB/20);  % Noise level

x = zeros(Lx,1);   % x (Nx1) input signal array definition

% Colored signal generation -----------------------------------------------
a = 0.95;
b = sqrt(1-a^2); 

% Adaptive filter definition ----------------------------------------------
% IIR_WSAF
M =  4;          % MA filter length
N =  3;          % AR filter length
K =  1;          % Affine Projection order; K=1 => NLMS
mu = 0.02;       % Learning rate for the linear filter
mQ = 0.02;       % Learning rate for control points
% WSAF
Mw = 5;          % Linear filter length
mu0 = 0.1;       % Learning rate APA
mQ0 = 0.1;       % Learning rate for control points
delta = 1e-2;    % APA regularization parameters
% Polynomial
muP = 0.005;     % Learning rate for the linear filter
mQP = 0.005;     % Learning rate for polynomial coeficients
% Volterra AF
MV1 = 15;        % Number of Volterra coefficients
Kv  = 1;         % Affine projection order
muV = 0.01;      % Learning rate for Volterra AF


% Spline activation function definition and initialization ----------------
afinit = 0;     % Init act. func. -1 0 ... (ONLY -1, bip.sig. or  0,linear ) 
aftype = 5;     % Kind act. f. -1 0 1 2 3 4 5  
aftype2 = 3;    % Kind act. f. -1 0 1 2 3 4 5  
Slope = 1;      % Slope
DeltaX = 0.2;   % Delta X
x_range = 2.0;  % Range limit
Pord = 5;       % Olynomial order

% Creating the nonlinearity -----------------------------------------------
af1 = create_activation_function( afinit, aftype, DeltaX, x_range, Slope, M);
af2 = create_activation_function( afinit, aftype, DeltaX, x_range, Slope, Mw);
af3 = create_activation_function( afinit, aftype2, DeltaX, x_range, Slope, M, Pord);


%% Initialization

% --- Model definition ----------------------------------------------------
apa_af1 = create_SPL_lms_iir_adaptive_filter_1(M,N,mu0,mQ0,delta,af1);  % IIR_WSAF LMS 
apa_af2 = create_SPL_apa_adaptive_filter_1(M,K,mu,mQ,delta,af2);        % WSAF APA  
%apa_af2 = create_SPL_lms_adaptive_filter_1(M,mu,mQ,delta,af);          % WSAF LMS
apa_af3 = create_SPL_lms_iir_adaptive_filter_1(M,N,muP,mQP,delta,af3);  % IIR Polynomial LMS
VF      = create_III_ord_Volterra_filter_1(MV1, Kv, muV, delta);        % III Order Volterra AF

% Initialize --------------------------------------------------------------
Nn = Lx + M + K;   % Total number of samples
for i = Lx+1:Nn
    x(i)=0;
end

e1 = zeros(Nn,1);    % Error array IIR_WSAF
e2 = zeros(Nn,1);    % Error array III-order Volterra AF
e3 = zeros(Nn,1);    % Error array WSAF
e4 = zeros(Nn,1);    % Error array IIR_Polynomial AF
y  = zeros(Nn,1);    % Output array 
dn = zeros(Nn,1);    % Noise output array 
d  = zeros(Nn,1);    % Desired signal array


%% Main loop --------------------------------------------------------------
disp('WSAF_IIR_Volterra_compare algorithm start ... ');
t = clock;

for n = 0 : nRun-1
    fprintf('Test nr. %d/%d\n', n+1, nRun);
    
    % Back & Tsoi NARMA system signal generation --------------------------
    % Back, A. D., Tsoi, A. C. (1993). A simplified gradient algorithm
    % for IIR synapse Multilayer Perceptron. Neural Networks, 5, 456-462.
    x = filter( sqrt(1-b^2), [1 -b], randn( size(x)));  % H(z) = a*z^-1/1-b*z^-1
    x = x/max(abs(x));     % Normalize x
    z = filter([0.0154 0.0462 0.0465 0.0154], [1 -1.99 1.572 -0.4583], x);
    d = sin(z);  %d = sin((pi/2)*z);

    % Noise generation ----------------------------------------------------
    dn = out_noise_level * randn( size(x) ) ;
    
    % SAF I.C. ------------------------------------------------------------
    apa_af1.b    = zeros(M,1);
    apa_af1.b(1) = 1;
    apa_af1.a    = zeros(N,1);
    apa_af1.w    = [apa_af1.b.' apa_af1.a.']';
    apa_af1.af   = create_activation_function(afinit, aftype, DeltaX, x_range, Slope, M);
    VF           = create_III_ord_Volterra_filter_1(MV1, Kv, muV, delta);  % III Order Volterra
    apa_af2.w(:) = 0;
    apa_af2.w(1) = 1;
    apa_af2.af   = create_activation_function(afinit, aftype, DeltaX, x_range, Slope, Mw);
    apa_af3.b    = zeros(M,1);
    apa_af3.b(1) = 1;
    apa_af3.a    = zeros(N,1);
    apa_af3.w    = [apa_af1.b.' apa_af1.a.']';
    
    % Set Nonlinear Func I.C. ---------------------------------------------
    apa_af3.af   = create_activation_function(afinit, aftype2, DeltaX, x_range, Slope, M, Pord);
   
    % Algorithms Comparison -----------------------------------------------   
    for k = 1 : Lx
        % Updating WSAF_IIR -----------------------------------------------
        [apa_af1, y(k), ee] = AF_LMS_WSPL_IIR_F(apa_af1, x(k), d(k) + dn(k));  % SAF 
        e1(k) = e1(k) +  ee^2;
        
        % Updating III-order Volterra AF ----------------------------------
        [VF, yDummy, ee] = AF_III_ord_VF_APA_F( VF, x(k), d(k) + dn(k) );      % III Order Volterra AF
        e2(k) = e2(k) +  ee^2;
        
        % Updating WSAF ---------------------------------------------------
        [apa_af2, y(k), ee] = AF_APA_WSPL_F(apa_af2, x(k), d(k) + dn(k));      % WSAF APA
        %[apa_af2, y(k), ee] = AF_LMS_WSPL_F(apa_af2, x(k), d(k) + dn(k));     % WSAF LMS
        e3(k) = e3(k) +  ee^2;
        
        % Updating IIR Polynomilal AF -------------------------------------
        [apa_af3, y(k), ee] = AF_LMS_WPOLY_IIR_F(apa_af3, x(k), d(k) + dn(k)); % IIR Polynomial AF
        e4(k) = e4(k) +  ee^2;
    end

end

% Mean square errors ------------------------------------------------------
e1 = e1/nRun;   % IIR_WSAF
e2 = e2/nRun;   % III-order Volterra AF
e3 = e3/nRun;   % WSAF
e4 = e4/nRun;   % IIR Polynomial AF

fprintf('\n');


%% Resuts

% Plot Spline functions ---------------------------------------------------
yLIM = 1.5;
xLIM = 3.0;
fig1 = figure('PaperSize',[20.98 29.68]);
box('on');
hold on;
hold('all');
ylim([-yLIM yLIM]);   
xlim([-xLIM xLIM]);   
grid on;
KK = 500;
yy1 = zeros(1,KK);
yy2 = zeros(1,KK);
yy3 = zeros(1,KK);
yy4 = zeros(1,KK);
xa1 = zeros(1,KK);
dx = 2*xLIM/KK;
xx = -xLIM;   
for k = 1:KK
    yy1(k) = ActFunc(xx, apa_af1.af);      % IIR WSAF   
    [VF,yy2(k)] = FW_III_ord_VF_F(VF, xx); % III-order Volterra AF
    yy3(k) = ActFunc(xx, apa_af2.af);      % WSAF
    yy4(k) = ActFunc(xx, apa_af3.af);      % IIR Polynomial AF
    xa1(k) = xx;
    xx = xx + dx;
end  
title('Profile of spline nonlinearity after learning'); 
xlabel('Linear combiner output {\its}[{\itn}] ');
ylabel('SAF output {\ity}[{\itn}]');
plot(xa1,yy1,'-.k','LineWidth',2);  
plot(xa1,yy2,'r','LineWidth',2); 
plot(xa1,yy3,'b','LineWidth',2); 
plot(xa1,yy4,'g','LineWidth',2); 
legend('IIR WSAF', 'Volterra', 'WSAF', 'IIR Poly');

% MSE dB ------------------------------------------------------------------
fig2 = figure('PaperSize',[10 15]);
box('on');hold on;hold('all');
ylim([-45 0]);   
grid on;
[bb,aa] = butter(3, 0.02 );
edb1 = 10*log10( e1 ); 
edb2 = 10*log10( e2 ); 
edb3 = 10*log10( e3 );
edb4 = 10*log10( e4 );
plot( filter(bb,aa,edb1),'-k','LineWidth',2);   % IIR_WSAF
plot( filter(bb,aa,edb2),'-.b','LineWidth',2);  % III-order Volterra AF
plot( filter(bb,aa,edb3),'--r','LineWidth',2);  % WSAF
plot( filter(bb,aa,edb4),':m','LineWidth',2);   % IIR Polynomial AF
title('Back and Tsoi NARMA model','FontSize', 12, 'FontWeight', 'demi');
xlabel('Samples','FontSize', 12, 'FontWeight', 'demi');
ylabel('MSE [dB]   10log({\itJ}({\bfw},{\bfq}))','FontSize', 12, 'FontWeight', 'demi');
legend('IIR WSAF', 'Volterra', 'WSAF', 'IIR Poly');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);

% Filter coefficients -----------------------------------------------------
fig3 = figure('PaperSize',[20.98 29.68]);
hold on;
plot(apa_af1.w,'k','LineWidth',2);
plot(apa_af2.w,'r','LineWidth',2);
plot(apa_af3.w,'m','LineWidth',2);
title('Filter coefficients');
xlabel('time {\itn}');
ylabel('Linear combiner coefficients');
legend('IIR WSAF', 'WSAF', 'IIR Poly');
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------