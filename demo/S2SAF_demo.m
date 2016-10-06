% S2SAF_demo.m
% 
% DEMO on nonlinear system identification based on a Sandwich 2 Spline
% Adaptive Filter (S2SAF).
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
disp('S2SAF_demo');
% -------------------------------------------------------------------------

%%  Parameters setting

% Input parameters --------------------------------------------------------
Lx = 30000;               % Length of the input signal
nRun = 10;                % Number of runs
out_noise_level_dB = 30;  % SNR
out_noise_level = 10^(-out_noise_level_dB/20);  % Noise level

x = zeros(Lx,1);   % x (Nx1) input signal array definition

% Colored signal generation -----------------------------------------------
a = 0.0;
b = sqrt(1-a^2);
% x = filter( b, [1 -a], randn( size(x))) ; % H(z) = b/(1+a*z^-1)
disp('...... ');


% Adaptive filter definition ----------------------------------------------
M1 = 7;            % Length of the first linear filter
M2 = 7;            % Length of the second linear filter
mu1 = 0.005;       % Learning rate of the first linear filter
mu2 = 0.009;       % Learning rate of the second linear filter
mQ0 = 0.004;       % Learning rate for ctrl points
if Lx < 30000,     % Batch for evaluating MSE
        B = 100;
else
        B = 4000;
end


% Spline activation function definition and initialization ----------------
afinit = 0;     % Init act. func. -1 0 ... (ONLY -1, bip.sig. or 0 = linear)
aftype = 4;     % Kind act. f. -1 0 1 2 4 5  
Slope  = 1;     % Slope
DeltaX = 0.2;   % Delta X
x_range = 2;    % Limit range


% Creating the nonlinearity -----------------------------------------------
af01 = create_activation_function(afinit, aftype, DeltaX, x_range, Slope, M2); % i.c.

%% Initialization

% --- Model definition ----------------------------------------------------
H1 = create_Sandwich2SPL_lms_adaptive_filter_1(M1,M2,mu1,mu2,mQ0,1e-2,af01); % Linear Filter Model LMS

% --- Target Definition ---------------------------------------------------
TH1 = create_Sandwich2SPL_lms_adaptive_filter_1(M1,M2,mu1,mu2,mQ0,1e-2,af01); % Target h Model LMS 

% TARGET: Nonlinear memoryless function implemented by HW-Spline interpolated LUT
Q0  = [    -2.20 
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

TH1.af.Q = Q0;
QL = length(Q0);  % Number of control points
 
% Linear filters ----------------------------------------------------------
TH1.w1  = [0.6  -0.4  0.25  -0.15  0.1  -0.05  0.001]';  % MA system to be identified
%TH1.w1  = [1    0    0      0     0     0     0    ]';  % Ideal impulse 
TH1.w2  = [1    0.5 -0.25   0.15  0.25  -0.10  0.050]';  % MA system to be identified
%TH1.w2  = [1    0    0      0     0     0     0    ]';  % Ideal impulse 

% --- SAF definition ------------------------------------------------------

% Initialize --------------------------------------------------------------  
N = Lx + M1 + 1;    % Total number of samples
for i = Lx+1:N
    x(i)=0;
end

d     = zeros(N,1);      % Desired signal array
dn    = zeros(N,1);      % Noise output array 
y     = zeros(N,1);      % Output array
e     = zeros(Lx,1);     % Error array 
em    = zeros(Lx,1);     % Mean square error array 
wm1   = zeros(M1,1);     % Mean value of w 
varW1 = zeros(M1,1);     % Variance value of w 
wm2   = zeros(M2,1);     % Mean value of w 
varW2 = zeros(M2,1);     % Variance value of w 
qm    = zeros(QL,1);     % Mean value Spline coeff
varQ  = zeros(QL,1);     % Variance value Spline coeff 


%% Main loop --------------------------------------------------------------
disp('S2SAF algorithm start ... ');
t = clock;

for n = 0 : nRun-1
    fprintf('Test nr. %d/%d\n', n+1, nRun);
    x = filter( b, [1 -a], randn( size(x))) ; % H(z) = b/(1+a*z^-1)
    dn = out_noise_level * randn( size(x) ) ;
       
    % SAF I.C. ------------------------------------------------------------
    H1.w1 (:) = 0;  
    H1.w1 (1) = 1;   % Set filter i.c.
    H1.w2 (:) = 0; 
    H1.w2 (1) = 1;   % Set filter i.c.
    
    % Set Activation Func I.C. --------------------------------------------
    H1.af.Q = af01.Q;
        
    % S2SAF Evaluation ----------------------------------------------------
    for k = 1 : Lx
        % Computing the desired output ------------------------------------
        [TH1, d(k)] = FW_Sandwich2SPL_F(TH1, x(k));      % Sandwich 2 model        

        % Updating S2SAF --------------------------------------------------
        [H1, y(k), e(k)] = AF_LMS_Sandwich2SPL_F(H1, x(k), d(k) + dn(k));   % Sandwich2SAF LMS
    end
    
    em  = em  + (e.^2);  % Squared error


    % SAF run-time mean and variance estimation ---------------------------
    wm1   = (1/(n+1))*H1.w1 + (n/(n+1))*wm1;
    varW1 = varW1 + (n/(n+1))*((TH1.w1 - wm1).^2);
    
    wm2   = (1/(n+1))*H1.w2 + (n/(n+1))*wm2;
    varW2 = varW2 + (n/(n+1))*((TH1.w2 - wm2).^2);
    
    qm   = (1/(n+1))*H1.af.Q + (n/(n+1))*qm;
    varQ = varQ + (n/(n+1))*((TH1.af.Q - qm).^2);  
     
end

em  = em/nRun;  % MSE
H1.af.Q = qm;

%--------------------------------------------------------------------------
% Average MSE evaluation
mse = mean(em(end-B-M1-1:end-M1-1));  % Average MSE
%--------------------------------------------------------------------------
fprintf('\n');


%% Results

% -------------------------------------------------------------------------
% Print table
% -------------------------------------------------------------------------
fprintf('\n');
fprintf('Number of iterations = %d\n',nRun); 
fprintf('Learning rates:  muW1 = %5.3f   muw2 = %5.3f   muQ\n', mu1, mu2, mQ0);
fprintf('a = %4.2f  b = %4.2f\n',a, b ); 
fprintf('Number of filter weights M1 = %d   M2 = %d\n', M1, M2); 
fprintf('Number of control points = %d\n', QL); 
fprintf('AF type = %d\n',aftype); 
fprintf('DeltaX = %4.2f\n',DeltaX); 
fprintf('SNR_dB  = %4.2f dB\n',out_noise_level_dB);
fprintf('Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse,10*log10(mse));
fprintf('\n');
fprintf('Mean and Variance Tables -----------------------------------------\n');
for i=1:QL
    fprintf('i=%2d  q0 =%5.2f    qm =%9.6f   varQ = %10.3e \n', i,  Q0(i), qm(i),  varQ(i) );
end
fprintf('\n');
for i=1:M1
    fprintf('i=%d  w01 =%5.2f    wm1 =%9.6f   varW1 = %10.3e \n', i, TH1.w1(i), wm1(i),  varW1(i)  );
end
fprintf('\n');
for i=1:M2
    fprintf('i=%d  w02 =%5.2f    wm2 =%9.6f   varW2 = %10.3e \n', i, TH1.w2(i), wm2(i),  varW2(i)  );
end
fprintf('END Sandwich 2 SAF test ------------------------------------------\n');
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Plotting figures
% -------------------------------------------------------------------------
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
xa1 = zeros(1,KK);
dx = 2*xLIM/KK;
xx = -xLIM;    
for k = 1:KK
    yy1(k) = ActFunc(xx, TH1.af);     
    yy2(k) = ActFunc(xx, H1.af);
    xa1(k) = xx;
    xx = xx + dx;
end  
title('Profile of spline nonlinearity after learning','FontSize', 12, 'FontWeight', 'demi');
xlabel('Linear combiner output {\its}[{\itn}] ','FontSize', 12, 'FontWeight', 'demi');
ylabel('SAF output {\ity}[{\itn}]','FontSize', 12, 'FontWeight', 'demi');
plot(xa1,yy1,'-.k','LineWidth',2);  % reference
plot(xa1,yy2,'r','LineWidth',2);    % adapted
legend('Target','Adapted');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);

% Filter coefficients -----------------------------------------------------
figure2 = figure('PaperSize',[20.98 29.68]);
title('Comparison of model and adapted linear filter','FontSize', 12, 'FontWeight', 'demi');
hold on;
grid on;
plot(TH1.w1,'-k','LineWidth',2.5);
plot(TH1.w2,'-r','LineWidth',2.5);
plot(H1.w1,'--b','LineWidth',2.5);
plot(H1.w2,'--m','LineWidth',2.5);
xlabel('samples {\itn}','FontSize', 12, 'FontWeight', 'demi');
ylabel('Linear combiner coefficients','FontSize', 12, 'FontWeight', 'demi');
legend('Model 1','Model 2','Adapted 1','Adapted 2');
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