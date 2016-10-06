% WSAF_IIR_demo.m
% 
% DEMO on nonlinear system identification (w0,q0) based on an IIR Wiener
% Spline Adaptive Filter (IIR_WSAF).
%
% This demo file implements the experimental test shown in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear System
% Identification using IIR Spline Adaptive Filters", Signal Processing, 
% Vol. 108, pp. 30-35, March 2015.
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2014/04/13$
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
disp('WSAF_IIR_demo');
% -------------------------------------------------------------------------

%%  Parameters setting

% Input parameters --------------------------------------------------------
Lx = 30000;               % Length of input signal
nRun = 10;                % Number of runs
out_noise_level_dB = 30;  % SNR
out_noise_level = 10^(-out_noise_level_dB/20);  % Noise level

x = zeros(Lx,1);          % x (Nx1) input signal array definition

% Colored signal generation -----------------------------------------------
a = 0.9;
b = sqrt(1-a^2); 
% x = filter( b, [1 -a], randn( size(x))) ; % H(z) = b/(1+a*z^-1)

% Adaptive filter definition ----------------------------------------------
M =  2;            % MA filter length
N =  3;            % AR filter length
K =  1;            % Affine Projection order; K=1 => NLMS
mu0  = 0.01;       % Learning rate for the linear IIR filter
mQ0  = 0.01;       % Learning rate for control points
if Lx < 30000,     % Batch for evaluating MSE
        B = 100;
else
        B = 4000;
end

% Spline activation function definition and initialization ----------------
afinit = 0;    % Init act. func. -1 0 ... (ONLY -1, bip.sig. or  0 =linear ) 
aftype = 5;    % Kind act. f. -1 0 1 2 4 5; (4 = CR-spline, 5 = B-spline )
Slope  = 1;    % Slope
DeltaX = 0.2;  % Delta X
x_range = 2;   % Range limit

% Creating the nonlinearity -----------------------------------------------
af0 = create_activation_function( afinit, aftype, DeltaX, x_range, Slope, 2);
af1 = create_activation_function( afinit, aftype, DeltaX, x_range, Slope, M);
 
% Printing the control points ---------------------------------------------
fprintf('af.lut_len = %d \n', af0.lut_len ); 
for j=1 : af0.lut_len       
    fprintf('%16.2f \n', af0.Q(j) ); 
end


%% Initialization

% --- Model definition ----------------------------------------------------
% apa_af0 = create_SPL_apa_iir_adaptive_filter_1(M,N,K,mu0,mQ0,1e-2,af0); % Model APA
apa_af0 = create_SPL_lms_iir_adaptive_filter_1(2,3,mu0,mQ0,1e-2,af0); % Model LMS

% TARGET: Nonlinear memoryless function implemented by Spline interpolated LUT
apa_af0.af.Q = [
           -2.20 
           -2.00 
           -1.80 
           -1.60 
           -1.40 
           -1.20 
           -1.00 
           -0.80 
           -0.91    
           -0.40   
           -0.20   
            0.05  
            0.0    
           -0.15
            0.58  
            1.00 
            1.00 
            1.20 
            1.40 
            1.60 
            1.80 
            2.00
            2.20 
];
QL = length(apa_af0.af.Q);   % Number of control points
 
% Linear filter -----------------------------------------------------------
apa_af0.b =  [0.6  -0.4].';                % MA model
apa_af0.a =  [0.2  -0.1  0.1].';           % AR model
apa_af0.w =  [apa_af0.b.'  apa_af0.a.'].'; % ARMA model

% --- SAF definition ------------------------------------------------------
apa_af1 = create_SPL_lms_iir_adaptive_filter_1(M,N,mu0,mQ0,1e-2,af1);  % SAF 

% Initialize --------------------------------------------------------------
Nn = Lx + M + K;    % Total samples
for i = Lx+1:Nn 
    x(i)=0;
end

dn = zeros(Nn,1);       % Noise output array 
d  = zeros(Nn,1);       % Desired signal array
y  = zeros(Nn,1);       % Output array 
e  = zeros(Nn,1);       % Error array
em = zeros(Nn,1);       % Mean square error array 
wm   = zeros(M+N,1);    % Mean value of w 
varW = zeros(M+N,1);    % Variance value of w 
qm   = zeros( QL, 1 );  % Mean value Spline coeff
varQ = zeros( QL, 1 );  % Variance value Spline coeff 


%% Main loop --------------------------------------------------------------
disp('WSAF_IIR algorithm start ... ');
t = clock;

for n = 0 : nRun-1
    fprintf('Test nr. %d/%d\n', n+1, nRun);
    x = filter( b, [1 -a], randn( size(x))) ;   % H(z) = b/(1+a*z^-1)
    dn = out_noise_level * randn( size(x) ) ;   % Noise
    
    % SAF I.C. ------------------------------------------------------------
    apa_af1.b    = zeros(M,1);
    apa_af1.b(1) = 1;
    apa_af1.a    = zeros(N,1);
    apa_af1.w    = [apa_af1.b.' apa_af1.a.']';
    
    % Set Nonlinear Func I.C. ---------------------------------------------
    apa_af1.af  = create_activation_function(afinit, aftype, DeltaX, x_range, Slope, M);
 
    % WSAF_IIR Evaluation -------------------------------------------------
    for k = 1 : Lx
        % Computing the desired output ------------------------------------
        [apa_af0, d(k)] = FW_WSPL_IIR_F(apa_af0, x(k));    % Model
        
        % Updating WSAF_IIR -----------------------------------------------
        [apa_af1, y(k), e(k)] = AF_LMS_WSPL_IIR_F(apa_af1, x(k), d(k) + dn(k));  % SAF LMS (Eqs. (17) and (19) )
    end
    
    em   = em + e.^2;  % Squared error
    
    % SAF run-time mean and variance estimation ---------------------------
    wm   = (1/(n+1))*apa_af1.w + (n/(n+1))*wm;
    varW = varW + (n/(n+1))*((apa_af1.w - wm).^2);
    
    qm   = (1/(n+1))*apa_af1.af.Q + (n/(n+1))*qm;
    varQ = varQ + (n/(n+1))*((apa_af1.af.Q - qm).^2);   
end
em  = em/nRun;  % MSE
apa_af1.af.Q = qm;

%--------------------------------------------------------------------------
% Average MSE evaluation
mse = mean(em(end-B-M-1:end-M-1));  % Average MSE
%--------------------------------------------------------------------------
fprintf('\n');


%% Results

%--------------------------------------------------------------------------
% Print tables
% -------------------------------------------------------------------------
fprintf('\n');
fprintf('Number of iterations = %d\n',nRun); 
fprintf('Learning rates:  muW = %5.3f   muQ = %5.3f\n', mu0, mQ0);
fprintf('a = %4.2f  b = %4.2f\n',a, b ); 
fprintf('Number of filter weights: M = %d and N = %d\n', M, N); 
fprintf('Number of control points = %d\n', QL); 
fprintf('AF type = %d\n',aftype); 
fprintf('DeltaX  = %4.2f\n',DeltaX); 
fprintf('SNR_dB  = %4.2f dB\n',out_noise_level_dB);
fprintf('Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse,10*log10(mse));
fprintf('\n');
fprintf('Mean and Variance Tables -----------------------------------------\n');
for i=1:QL
    fprintf('i=%2d  q0 =%5.2f    qm =%9.6f   varQ = %10.3e \n', i,  apa_af1.af.Q(i), qm(i), varQ(i) );
end
fprintf('\n');
fprintf('------------------------------------------------------------------\n');
for i=1:M
    fprintf('i=%d  w0 =%5.2f    wm =%9.6f   varW = %10.3e \n', i, apa_af1.w(i), wm(i), varW(i) );
end
fprintf('------------------------------------------------------------------\n');


% -------------------------------------------------------------------------
% Plotting figures
% -------------------------------------------------------------------------

% Plot Spline functions ---------------------------------------------------
yLIM = 1.5;
xLIM = 3.0;
figure1 = figure('PaperSize',[20.98 29.68]);
box('on');
hold on;
hold('all');
ylim([-yLIM yLIM]);   
xlim([-xLIM xLIM]);   
grid on;
KK = 500;
yy1 = zeros(1,KK);
yy2 = zeros(1,KK);
xa1 = zeros(1,KK);
dx = 2*xLIM/KK;
xx = -xLIM;   
for k = 1:KK
    yy1(k) = ActFunc(xx, apa_af0.af);   % Model 
    yy2(k) = ActFunc(xx, apa_af1.af);   % Adapted
    xa1(k)=xx;
    xx = xx + dx;
end  
xlabel('Linear combiner output {\its}[{\itn}] ','FontSize', 12, 'FontWeight', 'demi');
ylabel('SAF output {\ity}[{\itn}]','FontSize', 12, 'FontWeight', 'demi');
title('Profile of model and adapted nonlinearity \varphi(s[n])','FontSize', 12, 'FontWeight', 'demi');
plot(xa1,yy1,'-.r','LineWidth',2);  
plot(xa1,yy2,'k','LineWidth',2);
legend('Target','Adapted','Location','SouthEast');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');

% Filter coefficients -----------------------------------------------------
figure2 = figure('PaperSize',[20.98 29.68]);
hold on;
plot(apa_af0.w,'LineWidth',2);
plot(apa_af1.w,'m','LineWidth',2);
xlabel('time {\itn}','FontSize', 12, 'FontWeight', 'demi');
ylabel('Linear combiner coefficients','FontSize', 12, 'FontWeight', 'demi');
legend('Model','Adapted');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');

% MSE dB ------------------------------------------------------------------
figure3 = figure('PaperSize',[10 15]);
box('on');hold on;hold('all');
ylim([-out_noise_level_dB-5 10]);
grid on;
edb = 10*log10( em );
[bb,aa] = butter(2, 0.02 );
plot( filter(bb,aa,edb ),'Color',[1 0 0],'LineWidth',2);
noiseLevel(1: length(edb)-1 ) = -out_noise_level_dB;
plot( noiseLevel,'--','Color',[0 0 1],'LineWidth',2 );
title('Wiener IIR SAF convergence test','FontSize', 12, 'FontWeight', 'demi');
xlabel('Samples','FontSize', 12, 'FontWeight', 'demi');
ylabel('MSE [dB]   10log({\itJ}({\bfw},{\bfQ}))','FontSize', 12, 'FontWeight', 'demi');
legend('MSE', 'NoiseLevel' );
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);
% -------------------------------------------------------------------------

fprintf('END WSAF_demo ----------------------------------------------------\n');

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------