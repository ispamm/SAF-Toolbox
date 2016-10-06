% S2SAF_compare_demo.m
% 
% DEMO of comparison results on nonlinear system identification based on a 
% Sandwich 2 Spline Adaptive Filter (S2SAF).
%
% This demo file implements the comparison results shown in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Novel Cascade 
% Spline Architectures for the Identification of Nonlinear Systems", 
% IEEE Transactions on Circuits and Systems---I: Regular Papers, Vol. 62,
% No. 7, pp. 1825-1835, July 2015.
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2014/10/15$
% $Revision: 1.1$  $Date: 2016/09/01$
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
disp('S2SAF_compare_demo');
% -------------------------------------------------------------------------

%%  Parameters setting

% Input parameters --------------------------------------------------------
Lx = 30000;                % Length of the input signal
nRun = 10;                 % Number of runs
out_noise_level_dB = Inf;  % SNR
out_noise_level = 10^(-out_noise_level_dB/20);   % Noise level

x = randn(Lx,1);           % x (Nx1) input signal array definition

% Colored signal generation -----------------------------------------------
a = 0.1;
b = sqrt(1-a^2);
% x = filter( b, [1 -a], randn( size(x)));  % H(z) = b/(1+a*z^-1)
disp('...... ');


% Adaptive filter definition ----------------------------------------------
M1 = 5;            % Length of the first linear filter
M2 = 5;            % Length of the second linear filter
MV = 5;            % Length of the Volterra filter
M3 = 5;            % Length of the S1SAF linear filter
K =  1;            % Affine Projection order; K=1 => NLMS
Kv = 1;            % Volterra Affine Projection order; K=1 => NLMS
Pord = 3;          % Polynomial order
delta = 1e-2;      % APA regularization parameters
mu1 = 0.005;       % Learning rate for the first linear filter
mu2 = 0.009;       % Learning rate for the second linear filter
mQ0 = 0.04;        % Learning rate for the spline ctrl points
muV = 0.1;         % Learning rate for the Volterra filter
muM = 0.00001;     % Learning rate for the linear part of Mathews architecture
mQM = 0.00001;     % Learning rate for the nonlinear part of Mathews architecture
muJ1 = 0.0001;     % Learning rate for the first part of LNL Jenkins architecture
muJ2 = 0.0001;     % Learning rate for the second part of LNL Jenkins architecture
muJ3 = 0.001;      % Learning rate for the third part of LNL Jenkins architecture
muS  = 0.05;       % Learning rate for the linear part of S1SAF
mQ1S = 0.05;       % Learning rate for the first nonlinearity of S1SAF
mQ2S = 0.005;      % Learning rate for the second nonlinearity of S1SAF 
if Lx < 30000,     % Batch for evaluating MSE
        B = 100;
else
        B = 4000;
end


% Spline activation function definition and initialization ----------------
afinit = 0;     % Init act. func. -1 0 ... (ONLY -1, bip.sig. or 0 = linear)
aftype = 4;     % Kind act. f. -1 0 1 2 3 4 5  
aftype2 = 3;    % Kind act. f. -1 0 1 2 3 4 5  
Slope  = 1;     % Slope
DeltaX = 0.2;   % Delta X
x_range = 2;    % Range limit


% --- Nonlinearity definition ---------------------------------------------
af01 = create_activation_function(afinit, aftype, DeltaX, x_range, Slope, M2);         % S2SAF
af02 = create_activation_function(afinit, aftype2, DeltaX, x_range, Slope, M2, Pord);  % Polynomial Mathews
af03 = create_activation_function(afinit, aftype, DeltaX, x_range, Slope, M3);         % First S1SAF nonlinearity
af04 = create_activation_function(afinit, aftype, DeltaX, x_range, Slope, M3);         % Second S1SAF nonlinearity


% --- Models definition ---------------------------------------------------
H1 = create_Sandwich2SPL_lms_adaptive_filter_1(M1,M2,mu1,mu2,mQ0,delta,af01);          % S2SAF 
H2 = create_SPL_lms_adaptive_filter_1(M1,mu1,mQ0,delta,af01);                          % WSAF
H3 = create_SPL_lms_adaptive_filter_1(M2,mu2,mQ0,delta,af01);                          % HSAF
H4 = create_SPL_lms_adaptive_filter_1(M2,muM,mQM,delta,af02);                          % Polynomial Mathews
H5 = create_LNLJenkins_lms_adaptive_filter_1(M1,Pord,M2,muJ1,muJ2,muJ3,delta,af02);    % Hedge et al.
VF = create_III_ord_Volterra_filter_1(MV, Kv, muV, delta);                             % III Order Volterra
A1 = create_Sandwich1SPL_lms_adaptive_filter_1(M3,muS,mQ1S,mQ2S,delta,af03,af04);      % S1SAF 


% --- SAF definition ------------------------------------------------------
 
% Initialize --------------------------------------------------------------  
N = Lx + M1 + K;
for i = Lx+1:N
    x(i)=0;
end

dn  = zeros(N,1);        % Noise output array 
d   = zeros(N,1);        % Desired signal array
y1  = zeros(N,1);        % Output array
y2  = zeros(N,1);        % Output array
y3  = zeros(N,1);        % Output array
y4  = zeros(N,1);        % Output array
y5  = zeros(N,1);        % Output array
y6  = zeros(N,1);        % Output array
y7  = zeros(N,1);        % Output array
e1  = zeros(Lx,1);       % Error array 
e2  = zeros(Lx,1);       % Error array 
e3  = zeros(Lx,1);       % Error array 
e4  = zeros(Lx,1);       % Error array 
e5  = zeros(Lx,1);       % Error array 
e6  = zeros(Lx,1);       % Error array 
e7  = zeros(Lx,1);       % Error array
em1 = zeros(Lx,1);       % Mean square error 
em2 = zeros(Lx,1);       % Mean square error
em3 = zeros(Lx,1);       % Mean square error
em4 = zeros(Lx,1);       % Mean square error
em5 = zeros(Lx,1);       % Mean square error
em6 = zeros(Lx,1);       % Mean square error
em7 = zeros(Lx,1);       % Mean square error


%% Main loop --------------------------------------------------------------
disp('S2SAF_compare_demo start ... ');
t = clock;

for n = 0 : nRun-1
    fprintf('Test nr. %d/%d\n', n+1, nRun);
    
    % Nonlinear system signal generation ----------------------------------
    % 1998_Sicuranza-Mathews HPA Model ------------------------------------
    x = randn(Lx,1);           % x (Nx1) input signal array definition
    if a>0,
        x = filter(b, [1 -a], x) ;    % H(z) = a/1-b*z^-1
    end
    
    y = filter(conv([0.2851 0.5704 0.2851], [0.2851 0.5701 0.2851]), ...
               conv([1     -0.1024 0.4475],  [1    -0.0736 0.0408]), x);  % First linear filter 

    z = y./(1+abs(y.^2));      % Nonlinearity
        
    d = filter(conv([0.2025 0.288 0.2025],[0.2025 0.0034 0.2025]), ...
               conv([1     -1.01  0.5861],[1     -0.6591 0.1498]), z);  % Second linear filter
    
    dn = out_noise_level*randn( size(x) );    % Noisy desired signal
    
    
    % SAF I.C. ------------------------------------------------------------
    H1.w1 (:) = 0;
    H1.w1 (1) = 1;  % Set filter i.c.
    H1.w2 (:) = 0;
    H1.w2 (1) = 1;  % Set filter i.c.
    H2.w (:)  = 0;
    H2.w  (1) = 1;  % Set filter i.c.
    H3.w (:)  = 0;
    H3.w  (1) = 1;  % Set filter i.c.
    H4.w (:)  = 0;
    H4.w  (1) = 1;  % Set filter i.c.
    H5.w1(:)  = 0;
    H5.w1 (1) = 1;  % Set filter i.c.
    H5.w2(:)  = 0;
    H5.w2 (1) = 1;  % Set filter i.c.
    H5.w3(:)  = 0;
    H5.w3 (1) = 1;  % Set filter i.c.
    
    VF = create_III_ord_Volterra_filter_1(MV, Kv, muV, delta);  % III Order Volterra
    A1.w (:) = 0;
    A1.w (1) = 1;
    
    % Set Activation Func I.C. --------------------------------------------
    H1.af.Q = af01.Q;
    H2.af.Q = af01.Q;
    H3.af.Q = af01.Q;
    H4.af.Q = af02.Q;
    H5.af.Q = af02.Q;
    A1.af1.Q = af03.Q;
    A1.af2.Q = af04.Q;
        
    % Set Activation Func I.C. --------------------------------------------
    for k = 1 : Lx
        % Updating S2SAF --------------------------------------------------
        [H1, y1(k), e1(k)] = AF_LMS_Sandwich2SPL_F(H1, x(k), d(k) + dn(k) );   % S2SAF LMS
        
        % Updating VAF ----------------------------------------------------
        [VF, y2(k), e2(k)] = AF_III_ord_VF_APA_F( VF, x(k), d(k) + dn(k) );    % III Order Volterra
        
        % Updating WSAF ---------------------------------------------------
        [H2, y3(k), e3(k)] = AF_LMS_WSPL_F(H2, x(k), d(k) + dn(k) );           % WSAF LMS 
        
        % Updating HSAF ---------------------------------------------------
        [H3, y4(k), e4(k)] = AF_LMS_HSPL_F(H3, x(k), d(k) + dn(k) );           % HSAF LMS
        
        % Updating Polynomial Mathews -------------------------------------
        [H4, y5(k), e5(k)] = AF_LMS_MATHEWS_POLY_F(H4, x(k), d(k) + dn(k) );   % Polyomial Mathews 

        % Updating LNL Jenkins --------------------------------------------
        [H5, y6(k), e6(k)] = AF_LMS_LNLJenkins_F(H5, x(k), d(k) + dn(k) );     % LNL Jenkins 

        % Updating S1SAF --------------------------------------------------
        [A1, y7(k), e7(k)] = AF_LMS_Sandwich1SPL_F(A1, x(k), d(k) + dn(k) );   % S1SAF LMS
    end
    
    % Squared Error -------------------------------------------------------
    em1  = em1 + (e1.^2);
    em2  = em2 + (e2.^2);
    em3  = em3 + (e3.^2);
    em4  = em4 + (e4.^2);
    em5  = em5 + (e5.^2);
    em6  = em6 + (e6.^2);
    em7  = em7 + (e7.^2);
    
end

% MSE ---------------------------------------------------------------------
em1  = em1/nRun;
em2  = em2/nRun;
em3  = em3/nRun;
em4  = em4/nRun;
em5  = em5/nRun;
em6  = em6/nRun;
em7  = em7/nRun;

%--------------------------------------------------------------------------
% Average MSE evaluations
mse1 = mean(em1(end-B-M1-1:end-M1-1));  % Average MSE
mse2 = mean(em2(end-B-M1-1:end-M1-1));  % Average MSE
mse3 = mean(em3(end-B-M1-1:end-M1-1));  % Average MSE
mse4 = mean(em4(end-B-M1-1:end-M1-1));  % Average MSE
mse5 = mean(em5(end-B-M1-1:end-M1-1));  % Average MSE
mse6 = mean(em6(end-B-M1-1:end-M1-1));  % Average MSE
mse7 = mean(em7(end-B-M1-1:end-M1-1));  % Average MSE
%--------------------------------------------------------------------------
fprintf('\n');


%% Results

% -------------------------------------------------------------------------
% Print table
% -------------------------------------------------------------------------

fprintf('\n');
fprintf('Mean Square Errors -----------------------------------------------\n');
fprintf('  S2SAF Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse1,10*log10(mse1));
fprintf('    VAF Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse2,10*log10(mse2));
fprintf('   WSAF Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse3,10*log10(mse3));
fprintf('   HSAF Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse4,10*log10(mse4));
fprintf('Mathews Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse5,10*log10(mse5));
fprintf('Jenkins Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse6,10*log10(mse6));
fprintf('  S1SAF Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse7,10*log10(mse7));
fprintf('------------------------------------------------------------------\n');


% -------------------------------------------------------------------------
% Plotting figures
% -------------------------------------------------------------------------

yLIM = 1.5;
xLIM = 3.0;

% Without target
figure1_1 = figure('PaperSize',[10 15]);
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
yy5 = zeros(1,KK);
yy6 = zeros(1,KK);
yy7 = zeros(1,KK);
yy8 = zeros(1,KK);
xa1 = zeros(1,KK);
dx = 2*xLIM/KK;
xx = -xLIM;   
for k = 1:KK
    yy1(k) = ActFunc(xx, H1.af);            % S2SAF
    yy2(k) = xx - 0.5*xx.^3 + 0.02*xx.^5;   % Target
    yy3(k) = ActFunc(xx, H2.af);            % WSAF
    yy4(k) = ActFunc(xx, H3.af);            % HSAF
    [VF,yy5(k)] = FW_III_ord_VF_F(VF, xx);  % Volterra
    yy6(k) = ActFunc(xx, H4.af);            % Mathews
    yy7(k) = 2*xx./(1 + abs(xx.^2));        % Target
    xa1(k) = xx;
    xx = xx + dx;
end  
title('Profile of spline nonlinearity after learning','FontSize', 12, 'FontWeight', 'demi');
xlabel('Linear combiner output {\its}[{\itn}] ','FontSize', 12, 'FontWeight', 'demi');
ylabel('SAF output {\ity}[{\itn}]','FontSize', 12, 'FontWeight', 'demi');
plot(xa1,yy1,'-k','LineWidth',2);  % S2SAF
plot(xa1,yy2,'-.r','LineWidth',2); % Target
plot(xa1,yy3,'-g','LineWidth',2);  % WSAF
plot(xa1,yy4,'-c','LineWidth',2);  % HSAF
plot(xa1,yy5,'-m','LineWidth',2);  % Volterra
plot(xa1,yy6,'-y','LineWidth',2);  % Mathews
legend('S2SAF', 'Target', 'WSAF', 'HSAF', 'Volterra', 'Mathews');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);

% With target
figure1_2 = figure('PaperSize',[10 15]);
box('on');
hold on;
grid on;
ylim([-yLIM yLIM]);   
xlim([-xLIM xLIM]);   
plot(xa1,yy1,'-k','LineWidth',2);  % S2SAF
plot(xa1,yy7,'-.r','LineWidth',2); % Target
legend('S2SAF', 'Target');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);


% MSE dB ------------------------------------------------------------------
figure2 = figure('PaperSize',[10 15]);
box('on');
hold on;
grid on;
ylim([-out_noise_level_dB-5 0]);   
[bb,aa] = butter(3, 0.01 );
edb1 = 10*log10( em1 ); 
edb2 = 10*log10( em2 ); 
edb3 = 10*log10( em3 ); 
edb4 = 10*log10( em4 );
edb5 = 10*log10( em5 );
edb6 = 10*log10( em6 );
edb7 = 10*log10( em7 );
plot(filter(bb,aa,edb1),'k-','LineWidth',2);   % S2SAF
plot(filter(bb,aa,edb2),'b-.','LineWidth',2);  % Volterra 
plot(filter(bb,aa,edb3),'r--','LineWidth',2);  % WSAF
plot(filter(bb,aa,edb4),'m:','LineWidth',2);   % HSAF
plot(filter(bb,aa,edb5),'c-.','LineWidth',2);  % Mathews 
plot(filter(bb,aa,edb6),'g--','LineWidth',2);  % Jenkins
plot(filter(bb,aa,edb7),'y-','LineWidth',2);   % S1SAF
if out_noise_level_dB ~= Inf
    noiseLevel(1: length(edb1)-1 ) = -out_noise_level_dB;  % Noise level
    plot( noiseLevel,'--','Color',[0 0 1],'LineWidth',2 );
end
title('Comparisons of Sandwich SAF architectures','FontSize', 12, 'FontWeight', 'demi');
xlabel('Samples','FontSize', 12, 'FontWeight', 'demi');
ylabel('MSE [dB]   10log({\itJ}({\bfw},{\bfQ}))','FontSize', 12, 'FontWeight', 'demi');
if out_noise_level_dB ~= Inf
    legend('S2SAF', 'Volterra', 'WSAF', 'HSAF', 'Jeraj & Mathews', 'Hegde et al.', 'S1SAF', 'NoiseLevel' );
else
    legend('S2SAF', 'Volterra', 'WSAF', 'HSAF', 'Jeraj & Mathews', 'Hegde et al.', 'S1SAF' );
end
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);


% Filter coefficients -----------------------------------------------------
figure3 = figure('PaperSize',[20.98 29.68]);
hold on;
grid on;
title('Adapted linear filters','FontSize', 12, 'FontWeight', 'demi');
plot(H1.w1,'-b','LineWidth',2.5);     % First linear filter of S2SAF
plot(H1.w2,'-k','LineWidth',2.5);     % Second linear filter of S2SAF
plot(H2.w ,'-r','LineWidth',2.5);     % WSAF
plot(H3.w ,'-g','LineWidth',2.5);     % HSAF
plot(H4.w ,'-c','LineWidth',2.5);     % Mathews
xlabel('samples {\itn}','FontSize', 12, 'FontWeight', 'demi');
ylabel('Linear combiner coefficients','FontSize', 12, 'FontWeight', 'demi');
legend('S2SAF w1', 'S2SAF w2', 'WSAF', 'HSAF', 'Mathews');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------