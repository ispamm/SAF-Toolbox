% S1SAF_demo.m
% 
% DEMO on nonlinear system identification based on a Sandwich 1 Spline
% Adaptive Filter (S1SAF).
%
% This demo file implements the experimental test shown in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Novel Cascade 
% Spline Architectures for the Identification of Nonlinear Systems", 
% IEEE Transactions on Circuits and Systems---I: Regular Papers, Vol. 62,
% No. 7, pp. 1825-1835, July 2015.
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2014/10/15$
% $Revision: 1.1$  $Date: 2016/08/18$
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
disp('S1SAF_demo');
% -------------------------------------------------------------------------

%%  Parameters setting

% Input parameters --------------------------------------------------------
Lx = 30000;                % Length of input signal
nRun = 10;                 % Number of runs
out_noise_level_dB = 30;   % SNR
out_noise_level = 10^(-out_noise_level_dB/20);  % Noise level

x = zeros(Lx,1);           % x (Nx1) input signal array definition

% Colored signal generation -----------------------------------------------
a = 0.0;
b = sqrt(1-a^2);
% x = filter( b, [1 -a], randn( size(x))) ; % H(z) = b/(1+a*z^-1)
disp('...... ');

% Adaptive filter definition ----------------------------------------------
M =  7;            % Length of the linear filter
K =  1;            % Affine Projection order; K=1 => NLMS
mu0 = 0.05;        % Learning rate for linear filter
mQ1 = 0.05;        % Learning rate for ctrl points
mQ2 = 0.05;        % Learning rate for ctrl points
if Lx < 30000,     % Batch for evaluating MSE
        B = 100;
else
        B = 4000;
end


% Spline activation function definition and initialization ----------------
% NL 1
afinit1 = 0;     % Init act. func. -1 0 ... (ONLY -1, bip.sig. or 0 = linear) 
aftype1 = 4;     % Kind act. f. -1 0 1 2 4 5  
Slope1  = 1;     % Slope
DeltaX1 = 0.2;   % Delta X
x_range1 = 2;    % Range limit

% NL 2
afinit2 = 0;     % Init act. func. -1 0 ... (ONLY -1, bip.sig. or 0 = linear) 
aftype2 = 4;     % Kind act. f. -1 0 1 2 4 5  
Slope2 = 1;      % Slope
DeltaX2 = 0.2;   % Delta X
x_range2 = 2;    % Range limit


% Creating the nonlinearity -----------------------------------------------
af01 = create_activation_function(afinit1, aftype1, DeltaX1, x_range1, Slope1, M); % NL 1
af02 = create_activation_function(afinit2, aftype2, DeltaX2, x_range2, Slope2, M); % NL 2


%% Initialization

% --- Model definition ----------------------------------------------------
%H1 = create_Sandwich1SPL_apa_adaptive_filter_1(M,K,mu0,mQ1,mQ2,1e-2,af01,af02); % Linear Filter Model APA
H1 = create_Sandwich1SPL_lms_adaptive_filter_1(M,mu0,mQ1,mQ2,1e-2,af01,af02); % Linear Filter Model LMS

% --- Target Definition ---------------------------------------------------
%TH1 = create_Sandwich1SPL_apa_adaptive_filter_1(M,K,mu0,mQ1,mQ2,1e-2,af01,af02); % Target h Model APA
TH1 = create_Sandwich1SPL_lms_adaptive_filter_1(M,mu0,mQ1,mQ2,1e-2,af01,af02); % Target h Model LMS 

% TARGET: Nonlinear memoryless function implemented by HW-Spline interpolated LUT
% NL 1
Q01 = [    -2.20 
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
            0.00 
           -0.40
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

% NL 2
Q02 = [    -2.20 
           -2.00 
           -1.80 
           -1.60 
           -1.40 
           -1.20 
           -1.00 
           -0.80 
           -0.60 
           -0.10 
           -0.20 
            0.00
            0.02
            0.40
            0.60
            0.80
            1.35 
            1.20 
            1.40 
            1.60 
            1.80 
            2.00
            2.20 
];

TH1.af1.Q = Q01;
TH1.af2.Q = Q02;
QL1 = length(Q01);  % Number of control points
QL2 = length(Q02);  % Number of control points
 
% Linear filter -----------------------------------------------------------
TH1.w  = [0.6  -0.4  0.25  -0.15  0.1  -0.05  0.001]';  % MA system to be identified

% --- SAF definition ------------------------------------------------------
 
% Initialize --------------------------------------------------------------  
N = Lx + M + K;     % Total samples
for i = Lx+1:N 
    x(i) = 0;
end

dn    = zeros(N,1);        % Noise output array 
d     = zeros(N,1);        % Desired signal array
y     = zeros(N,1);        % Output array
e     = zeros(Lx,1);       % Error array 
em    = zeros(Lx,1);       % Mean square error 
wm    = zeros( M,1 );      % Mean value of w 
varW  = zeros( M,1 );      % Variance value of w 
qm1   = zeros( QL1, 1 );   % Mean value Spline 1 coeff
varQ1 = zeros( QL1, 1 );   % Variance value Spline 1 coeff 
qm2   = zeros( QL2, 1 );   % Mean value Spline 2 coeff
varQ2 = zeros( QL2, 1 );   % Variance value Spline 2 coeff 


%% Main loop --------------------------------------------------------------
disp('S1SAF algorithm start ... ');
t = clock;

for n = 0 : nRun-1
    fprintf('Test nr. %d/%d\n', n+1, nRun);
    x = filter( b, [1 -a], randn( size(x))) ; % H(z) = b/(1+a*z^-1)
    dn = out_noise_level * randn( size(x) ) ;
       
    % SAF I.C. ------------------------------------------------------------
    H1.w (:) = 0;
    H1.w (1) = 0.1;    % Set filter I.C. 
       
    % Set Activation Func I.C. --------------------------------------------
    H1.af1.Q = af01.Q;
    H1.af2.Q = af02.Q;
    
    % S1SAF Evaluation ----------------------------------------------------
    for k = 1 : Lx
        % Computing the desired output ------------------------------------
        [TH1, d(k), s, r] = FW_Sandwich1SPL_F(TH1, x(k));   % Sandwich 1 model        
        
        % Updating S1SAF --------------------------------------------------
        %[H1, y(k), e(k)] = AF_APA_Sandwich1SPL_F(H1, x1(k), d(k) + dn(k) );  % SandwichSAF APA 
        [H1, y(k), e(k)] = AF_LMS_Sandwich1SPL_F(H1, x(k), d(k) + dn(k) );   % SandwichSAF LMS      
    end
    
    em  = em  + (e.^2);     % Squared Error


    % SAF run-time mean and variance estimation ---------------------------
    wm   = (1/(n+1))*H1.w + (n/(n+1))*wm;
    varW = varW + (n/(n+1))*((TH1.w - wm).^2);
    
    qm1   = (1/(n+1))*H1.af1.Q + (n/(n+1))*qm1;
    varQ1 = varQ1 + (n/(n+1))*((TH1.af1.Q - qm1).^2);  
    
    qm2   = (1/(n+1))*H1.af2.Q + (n/(n+1))*qm2;
    varQ2 = varQ2 + (n/(n+1))*((TH1.af2.Q - qm2).^2);  
    
end

em = em/nRun;  % MSE
H1.af1.Q = qm1;
H1.af2.Q = qm2;

%--------------------------------------------------------------------------
% Average MSE evaluation
mse = mean(em(end-B-M-1:end-M-1));  % Average MSE
%--------------------------------------------------------------------------
fprintf('\n');


%% Results

% -------------------------------------------------------------------------
% Print table
% -------------------------------------------------------------------------
fprintf('\n');
fprintf('Number of iterations = %d\n',nRun); 
fprintf('Learning rates:  muW = %5.3f   muQ1 = %5.3f   muQ2\n', mu0, mQ1, mQ2);
fprintf('a = %4.2f  b = %4.2f\n',a, b ); 
fprintf('Number of filter weights = %d\n', M); 
fprintf('Number of control points = %d and %d\n', QL1, QL2); 
fprintf('AF type 1 = %d   AF type 2 = %d\n',aftype1, aftype2); 
fprintf('DeltaX1 = %4.2f    DeltaX2 = %4.2f\n',DeltaX1, DeltaX2); 
fprintf('SNR_dB  = %4.2f dB\n',out_noise_level_dB);
fprintf('Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse,10*log10(mse));
fprintf('\n');
fprintf('Mean and Variance Tables -----------------------------------------\n');
for i=1:QL1
    fprintf('i=%2d  q0 =%5.2f    qm =%9.6f   varQ = %10.3e \n', i,  Q01(i), qm1(i),  varQ1(i) );
end
fprintf('\n');
for i=1:QL2
    fprintf('i=%2d  q0 =%5.2f    qm =%9.6f   varQ = %10.3e \n', i,  Q02(i), qm2(i),  varQ2(i) );
end
fprintf('\n');
for i=1:M
    fprintf('i=%d  w0 =%5.2f    wm =%9.6f   varW = %10.3e \n', i, TH1.w(i), wm(i),  varW(i)  );
end
fprintf('------------------------------------------------------------------\n');


% -------------------------------------------------------------------------
% Plotting figures
% -------------------------------------------------------------------------

% Plot Spline functions ---------------------------------------------------
yLIM = 1.5;
xLIM = 3.0;
figure1 = figure('PaperSize',[10 15]);
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
    yy1(k) = ActFunc(xx, TH1.af1);     
    yy2(k) = ActFunc(xx, TH1.af2);
    yy3(k) = ActFunc(xx, H1.af1);
    yy4(k) = ActFunc(xx, H1.af2);
    xa1(k) = xx;
    xx = xx + dx;
end  
title('Profile of spline nonlinearity after learning','FontSize', 12, 'FontWeight', 'demi');
xlabel('Linear combiner output {\its}[{\itn}] ','FontSize', 12, 'FontWeight', 'demi');
ylabel('SAF output {\ity}[{\itn}]','FontSize', 12, 'FontWeight', 'demi');
plot(xa1,yy1,'-.k','LineWidth',2);  % reference 1
plot(xa1,yy2,'-.b','LineWidth',2);  % reference 2
plot(xa1,yy3,'r','LineWidth',2);    % adapted 1
plot(xa1,yy4,'m','LineWidth',2);    % adapted 2
legend('Target 1','Target 2','SPL 1','SPL 2');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);

% Filter coefficients -----------------------------------------------------
figure2 = figure('PaperSize',[20.98 29.68]);
hold on;
grid on;
plot(TH1.w,'-k','LineWidth',2.5);
plot(H1.w,'--r','LineWidth',2.5);
title('Comparison of model and adapted linear filter','FontSize', 12, 'FontWeight', 'demi');
xlabel('samples {\itn}','FontSize', 12, 'FontWeight', 'demi');
ylabel('Linear combiner coefficients','FontSize', 12, 'FontWeight', 'demi');
legend('Model','Adapted');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);

% MSE dB ------------------------------------------------------------------
figure3 = figure('PaperSize',[10 15]);
box('on');
hold on;
grid on;
ylim([-out_noise_level_dB-5 0]);   
[bb,aa] = butter(3, 0.01 );
edb1 = 10*log10( em ); 
plot( filter(bb,aa,edb1 ),'Color',[1 0 0],'LineWidth',2);
noiseLevel(1: length(edb1)-1 ) = -out_noise_level_dB;
plot( noiseLevel,'--','Color',[0 0 1],'LineWidth',2 );
title('Sandwich SAF convergence test','FontSize', 12, 'FontWeight', 'demi');
xlabel('Samples','FontSize', 12, 'FontWeight', 'demi');
ylabel('MSE [dB]   10log({\itJ}({\bfw},{\bfQ}))','FontSize', 12, 'FontWeight', 'demi');
legend('MSE', 'NoiseLevel' );
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------