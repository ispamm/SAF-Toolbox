% WSAF_Volterra_compare_demo.m
% 
% DEMO of comparison results on nonlinear system identification based on a 
% Wiener Spline Filter (WSAF) with a III-order Volterra Adaptive Filter.
%
% This demo file implements the experimental test shown in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear Spline 
% Adaptive Filtering", Signal Processing, Vol. 93, No. 4, pp. 772-783, 
% April 2013.
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2012/05/15$
% $Revision: 1.1$  $Date: 2016/09/06$
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

clear all; 
close all;
disp('WSAF_Volterra_compare_demo');
% -------------------------------------------------------------------------

%%  Parameters setting

% Input parameters --------------------------------------------------------
Lx = 100000;                   % Length of the input signal
nRun = 10;                     % Number of runs
out_noise_level_dB = Inf;      % SNR
out_noise_level = 10^(-out_noise_level_dB/20);   % Noise level

x = zeros(Lx,1);               % x (Nx1) input signal array definition

% Colored signal generation -----------------------------------------------
a = 0.9;
b = sqrt(1-a^2);
% x = filter(b, [1 -a], randn( size(x)));   % H(z) = b/(1+a*z^-1)
disp('...... ');

% Adaptive filter definition ----------------------------------------------
M   = 15;           % Length of the linear filter
MV  = 15;           % Length of Volterra filter
K   = 1;            % Affine Projection order; K=1 => NLMS
KV  = 1;            % Volterra Affine Projection order; K=1 => NLMS
mu  = 0.001;        % Learning rate for the linear filter
mQ  = 0.001;        % Learning rate for ctrl points
muV = 0.001;        % Learning rate for Volterra
delta = 1e-2;       % APA regularization parameters
if Lx < 30000,      % Batch for evaluating MSE
        B = 100;
else
        B = 4000;
end

% Spline activation function definition and initialization ----------------
afinit  = 0;       % Kind of init act. f. -1 0 1 2 */
aftype  = 5;       % Kind act. f. -1 0 1 2 4 5 */
Slope   = 1;       % Slope
DeltaX  = 0.2;     % DeltaX
x_range = 8;       % Range limit


%% Initialization

% Creating the nonlinearity -----------------------------------------------
af = create_activation_function(afinit, aftype, DeltaX, x_range, Slope, M); 

% --- Models definition ---------------------------------------------------
apa_af = create_SPL_apa_adaptive_filter_1(M, K, mu, mQ, delta, af);   % WSAF APA
VF     = create_III_ord_Volterra_filter_1(MV, KV, muV, delta);        % III Order Volterra


% Initialize --------------------------------------------------------------
N = Lx + M + K;
for i = Lx+1:N
    x(i)=0;
end

dn = zeros(N,1);        % Noise output array 
d  = zeros(N,1);        % Desired signal array
y1  = zeros(N,1);       % Output array 
y2  = zeros(N,1);       % Output array 
e1  = zeros(Lx,1);      % Error array
e2  = zeros(Lx,1);      % Error array 
em1 = zeros(Lx,1);      % Mean square error 
em2 = zeros(Lx,1);      % Mean square error


%% Main loop --------------------------------------------------------------
disp('WSAF_Volterra_compare algorithm start ... ');
t = clock;

for n = 0 : nRun-1
    fprintf('Test nr. %d/%d\n', n+1, nRun);
    
    % Nonlinear system signal generation ----------------------------------
    x = randn(Lx,1);           % x (Nx1) input signal array definition
    if a>0,
        x = filter( a, [1 -a], x);               % Colored input
    end
    y = filter(conv([0.2851 0.5704 0.2851], [0.2851 0.5701 0.2851]), ...
               conv([1     -0.1024 0.4475],  [1     -0.0736 0.0408]), x) ; 
    z = y./(1+abs(y.^2));   
    d = filter(conv([0.2025 0.288 0.2025],[0.2025 0.0034 0.2025]), ...
               conv([1     -1.01  0.5861],[1     -0.6591 0.1498]), z);
    dn = out_noise_level*randn( size(x) );       % Noisy desired signal
      
    % SAF I.C. ------------------------------------------------------------
    apa_af.w (:) = 0;
    apa_af.w (1) = 0.1;
    
    % Set Activation Func I.C. --------------------------------------------
    apa_af.af.Q = af.Q;


    % Models Evaluation ----------------------------------------------------
    for k = 1 : Lx
    
        % Updating WSAF ---------------------------------------------------
        [apa_af, y1(k), e1(k)] = AF_APA_WSPL_F(apa_af, x(k), d(k) + dn(k) );    % WSAF
        
        % Updating VAF ----------------------------------------------------
        [VF, y2(k), e2(k)] = AF_III_ord_VF_APA_F( VF, x(k), d(k) + dn(k) );     % III-order Volterra
    end
    
    em1  = em1 + (e1.^2);   % Squared error
    em2  = em2 + (e2.^2);   % Squared error
    
end

em1 = em1/nRun;      % MSE
em2 = em2/nRun;      % MSE

%--------------------------------------------------------------------------
% Average MSE evaluation
mse1 = mean(em1(end-B-M-1:end-M-1));  % Average MSE
mse2 = mean(em2(end-B-M-1:end-M-1));  % Average MSE
%--------------------------------------------------------------------------
fprintf('\n');


%% Results

% -------------------------------------------------------------------------
% Print table
% -------------------------------------------------------------------------

fprintf('\n');
fprintf('Mean Square Errors -----------------------------------------------\n');
fprintf('WSAF Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse1,10*log10(mse1));
fprintf(' VAF Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse2,10*log10(mse2));
fprintf('------------------------------------------------------------------\n');


% -------------------------------------------------------------------------
% Plotting figures
% -------------------------------------------------------------------------

% Plot Spline functions ---------------------------------------------------
 

% MSE dB ------------------------------------------------------------------
figure2 = figure('PaperSize',[20.98 29.68]);
box('on');
hold on;
grid on;
ylim([-out_noise_level_dB-5 0]);   
[bb,aa] = butter(2, 0.005);
edb1 = 10*log10( em1 );
edb2 = 10*log10( em2 );
plot(filter(bb,aa,edb1),'Color',[0 0 0],'LineWidth',2);       % WSAF
plot(filter(bb,aa,edb2),'Color',[1 0 0],'LineWidth',2);       % Volterra 
title('Comparisons of MSE','FontSize', 12, 'FontWeight', 'demi');
xlabel('Samples','FontSize', 12, 'FontWeight', 'demi');
ylabel('MSE [dB]   10log({\itJ}({\bfw},{\bfQ}))','FontSize', 12, 'FontWeight', 'demi');
legend('WSAF','3-rd Order Volterra');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2012
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------