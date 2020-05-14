%% Heat2_Apr19_figs - Creation of baseline comparison figures
% Includes straight baseline, spline baseline connecting linear sections,
% Peak area analysis figures using both methods

%% Load data
% Control samples
% Editable filename...
path_main = 'C:\Users\Rachel\Documents\Rachel\MIT\MIT Research\Short Lab\';
path_specific = 'PTFE\Irradiated April 2019\CSV\';
s1 = 'ptfe_postirr_ctrl_';
s2 = '1';
s3 = '_heat2.csv';
sample = [s1 s2 s3];
file = [path_main path_specific sample];

SampleNum = 20;
[T, ~] = importDataCSV(file);
x = length(T);
T = zeros(length(T),SampleNum);
Q = T;

for i = 1:SampleNum
    c = num2str(i);
    sample = [s1 c s3];
    file = [path_main path_specific sample];
    [T_temp, Q_temp] = importDataCSV(file);
    T(:,i) = T_temp(1:x);
    Q(:,i) = Q_temp(1:x); 
end

T_min = 260;  
T_max = 360; 
N = length(find(T(:,1) > T_min & T(:,1) < T_max));
T_int = zeros(N,SampleNum);
Q_int = T_int;
for i = 1:SampleNum
    run = [T(:,i), Q(:,i)];
    [T_temp, Q_temp] = Data_interp(run, T_min, T_max,N);
    T_int(:,i) = T_temp;
    Q_int(:,i) = Q_temp;
end


% 1e08 samples
% Editable filename...
s1 = 'ptfe_postirr_1e08_';
SampleNum = 12;

for i = 1:SampleNum
    c = num2str(i);
    sample = [s1 c s3];
    file = [path_main path_specific sample];
    [T_temp, Q_temp] = importDataCSV(file);
    T(:,i) = T_temp(1:x);
    Q(:,i) = Q_temp(1:x); 
end

for i = 1:SampleNum
    run = [T(:,i), Q(:,i)];
    [T_temp, Q_temp] = Data_interp(run, T_min, T_max,N);
    T_int(:,i+20) = T_temp;
    Q_int(:,i+20) = Q_temp;
end


% 1e09 samples
% Editable filename...
s1 = 'ptfe_postirr_1e09_';
SampleNum = 12;

for i = 1:SampleNum
    c = num2str(i);
    sample = [s1 c s3];
    file = [path_main path_specific sample];
    [T_temp, Q_temp] = importDataCSV(file);
    T(:,i) = T_temp(1:x);
    Q(:,i) = Q_temp(1:x); 
end

for i = 1:SampleNum
    run = [T(:,i), Q(:,i)];
    [T_temp, Q_temp] = Data_interp(run, T_min, T_max,N);
    T_int(:,i+32) = T_temp;
    Q_int(:,i+32) = Q_temp;
end

% 1e10 samples
% Editable filename...
s1 = 'ptfe_postirr_1e10_';
SampleNum = 12;

for i = 1:SampleNum
    c = num2str(i);
    sample = [s1 c s3];
    file = [path_main path_specific sample];
    [T_temp, Q_temp] = importDataCSV(file);
    T(:,i) = T_temp(1:x);
    Q(:,i) = Q_temp(1:x); 
end

for i = 1:SampleNum
    run = [T(:,i), Q(:,i)];
    [T_temp, Q_temp] = Data_interp(run, T_min, T_max,N);
    T_int(:,i+44) = T_temp;
    Q_int(:,i+44) = Q_temp;
end

% 1e11 samples
% Editable filename...
s1 = 'ptfe_postirr_1e11_';
SampleNum = 12;

for i = 1:SampleNum
    c = num2str(i);
    sample = [s1 c s3];
    file = [path_main path_specific sample];
    [T_temp, Q_temp] = importDataCSV(file);
    T(:,i) = T_temp(1:x);
    Q(:,i) = Q_temp(1:x); 
end

for i = 1:SampleNum
    run = [T(:,i), Q(:,i)];
    [T_temp, Q_temp] = Data_interp(run, T_min, T_max,N);
    T_int(:,i+56) = T_temp;
    Q_int(:,i+56) = Q_temp;
end

%% Spline Baseline
% Centered around peak

% Need to edit all following code and then repeat for cooling curves.

% Modify Heat analysis: 
% Add cooling data to Load Data
% Review differences between heat and cool peaks
% Repeat analysis for cooling data

% 1. find peak
% 2. Isolate N=-30/+15*60 point around peak
% 3. Run spline fit on what's left
% 4. Check how it did visually
% 5. Check relative error in final result
% 6. Compare to previous
% 7. Repeat for cooling
% 8. Make figures

% What I did before was fit a straight a line to the segments first

% 60 data points per degree
left_end = 30; % degrees from peak
right_end = 15; % degrees from peak
f = 2;

% Setup correct ranges

[~,index] = min(Q_int(:,1));
left = index - left_end*60;
right = index + right_end*60;
num = right-left+1;
mid_left = index - left_end*60/f;
mid_right = index + right_end*60/f;  
num_fit = (right-mid_right)+(mid_left-left)+2;

% Initialize holders for all correctly ranged data
T_peak = zeros(1,68);
T_corrected = zeros(num,68);
Q_temp = zeros(num,68);
T_fit = zeros(num_fit,68);
Q_fit = zeros(num_fit,68);
Q_baseline = zeros(num,68);
Q_corrected = zeros(num,68);

for j = 1:68
    % Find range
    [~,index] = min(Q_int(:,j));
    left = index - left_end*60;
    right = index + right_end*60;
    mid_left = index - left_end*60/f;
    mid_right = index + right_end*60/f;  
  
    % Save range and peak temperature
    Q_temp(:,j) = Q_int(left:right,j);
    T_corrected(:,j) = T_int(left:right,j);
    T_peak(j) = T_int(index,j);
    
    % Fit left and right straights to line
    T_right = T_int(mid_right:right,j);
    T_left = T_int(left:mid_left,j);
    T_fit(:,j) = [T_right; T_left];
    Q_right = Q_int(mid_right:right,j);
    Q_left = Q_int(left:mid_left,j);
    fit_right = [ones(length(T_right),1) T_right]\Q_right;
    fit_left = [ones(length(T_left),1) T_left]\Q_left;
    Q_fit(:,j) = [fit_right(1)+fit_right(2)*T_right; fit_left(1)+fit_left(2)*T_left];
    
    % Create spline
    Q_baseline(:,j) = spline(T_fit(:,j),Q_fit(:,j),T_corrected(:,j));
    
    % Subtract baseline from data
    Q_corrected(:,j) = Q_temp(:,j) - Q_baseline(:,j);    
end

%% Plot figures
f_data = figure;
f_corrected = figure;

figure(f_data)
hold on
% plot(T_corrected(:,1),Q_temp(:,1),'b','LineWidth',0.25)
% plot(T_fit(:,1),Q_fit(:,1),'b.','LineWidth',0.75)
% plot(T_corrected(:,1),Q_baseline(:,1),'b--','LineWidth',0.5)
plot(T_corrected(:,1:20),Q_temp(:,1:20),'b','LineWidth',0.25)
plot(T_corrected(:,21:32),Q_temp(:,21:32),'r','LineWidth',0.25)
plot(T_corrected(:,33:44),Q_temp(:,33:44),'g','LineWidth',0.25)
plot(T_corrected(:,45:56),Q_temp(:,45:56),'k','LineWidth',0.25)
plot(T_corrected(:,57:68),Q_temp(:,57:68),'m','LineWidth',0.25)
plot(T_fit(:,1:20),Q_fit(:,1:20),'b.','LineWidth',0.75)
plot(T_fit(:,21:32),Q_fit(:,21:32),'r.','LineWidth',0.75)
plot(T_fit(:,33:44),Q_fit(:,33:44),'g.','LineWidth',0.75)
plot(T_fit(:,45:56),Q_fit(:,45:56),'k.','LineWidth',0.75)
plot(T_fit(:,57:68),Q_fit(:,57:68),'m.','LineWidth',0.75)
plot(T_corrected(:,1:20),Q_baseline(:,1:20),'b','LineWidth',0.5)
plot(T_corrected(:,21:32),Q_baseline(:,21:32),'r','LineWidth',0.5)
plot(T_corrected(:,33:44),Q_baseline(:,33:44),'g','LineWidth',0.5)
plot(T_corrected(:,45:56),Q_baseline(:,45:56),'k','LineWidth',0.5)
plot(T_corrected(:,57:68),Q_baseline(:,57:68),'m','LineWidth',0.5)

title('Creating Baseline')
xlabel('Temperature [^{\circ}C]')
ylabel('Heat Flow [W/g]')
h = zeros(5, 1);
h(1) = plot(NaN,NaN,'b', 'LineWidth', 1);
h(2) = plot(NaN,NaN,'r', 'LineWidth', 1);
h(3) = plot(NaN,NaN,'g', 'LineWidth', 1);
h(4) = plot(NaN,NaN,'k', 'LineWidth', 1);
h(5) = plot(NaN,NaN,'m', 'LineWidth', 1);
legend(h, 'Unirradiated','10^{8} ions/cm^2','10^{9} ions/cm^2','10^{10} ions/cm^2','10^{11} ions/cm^2');


figure(f_corrected)
hold on
plot(T_corrected(:,1:20),Q_corrected(:,1:20),'b','LineWidth',0.5)
plot(T_corrected(:,21:32),Q_corrected(:,21:32),'r','LineWidth',0.5)
plot(T_corrected(:,33:44),Q_corrected(:,33:44),'g','LineWidth',0.5)
plot(T_corrected(:,45:56),Q_corrected(:,45:56),'k','LineWidth',0.5)
plot(T_corrected(:,57:68),Q_corrected(:,57:68),'m','LineWidth',0.5)

title('Baseline Corrected Data')
xlabel('Temperature [^{\circ}C]')
ylabel('Heat Flow [W/g]')
h = zeros(5, 1);
h(1) = plot(NaN,NaN,'b', 'LineWidth', 1);
h(2) = plot(NaN,NaN,'r', 'LineWidth', 1);
h(3) = plot(NaN,NaN,'g', 'LineWidth', 1);
h(4) = plot(NaN,NaN,'k', 'LineWidth', 1);
h(5) = plot(NaN,NaN,'m', 'LineWidth', 1);
legend(h, 'Unirradiated','10^{8} ions/cm^2','10^{9} ions/cm^2','10^{10} ions/cm^2','10^{11} ions/cm^2');

%% Plot figures
f_data = figure;
f_corrected = figure;

figure(f_data)
hold on
plot(T_corrected(:,1:20),Q_temp(:,1:20),'b','LineWidth',1)
plot(T_corrected(:,21:32),Q_temp(:,21:32),'r','LineWidth',1)
plot(T_corrected(:,33:44),Q_temp(:,33:44),'g','LineWidth',1)
% plot(T_corrected(:,45:56),Q_temp(:,45:56),'k','LineWidth',1)
% plot(T_corrected(:,57:68),Q_temp(:,57:68),'m','LineWidth',1)

title('Raw Data')
xlabel('Temperature [^{\circ}C]')
ylabel('Heat Flow [W/g] (exo up)')
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'b', 'LineWidth', 1);
h(2) = plot(NaN,NaN,'r', 'LineWidth', 1);
h(3) = plot(NaN,NaN,'g', 'LineWidth', 1);
% h(4) = plot(NaN,NaN,'k', 'LineWidth', 1);
% h(5) = plot(NaN,NaN,'m', 'LineWidth', 1);
% legend(h, 'Unirradiated');
% legend(h, 'Unirradiated','10^{8} ions/cm^2');
legend(h, 'Unirradiated','10^{8} ions/cm^2','10^{9} ions/cm^2');
% legend(h, 'Unirradiated','10^{8} ions/cm^2','10^{9} ions/cm^2','10^{10} ions/cm^2');
% legend(h, 'Unirradiated','10^{8} ions/cm^2','10^{9} ions/cm^2','10^{10} ions/cm^2','10^{11} ions/cm^2');


% figure(f_corrected)
% hold on
% plot(T_corrected(:,1:20),Q_corrected(:,1:20),'b','LineWidth',1)
% plot(T_corrected(:,21:32),Q_corrected(:,21:32),'r','LineWidth',1)
% plot(T_corrected(:,33:44),Q_corrected(:,33:44),'g','LineWidth',1)
% plot(T_corrected(:,45:56),Q_corrected(:,45:56),'k','LineWidth',1)
% plot(T_corrected(:,57:68),Q_corrected(:,57:68),'m','LineWidth',1)
% 
% title('Baseline Subtracted from Data')
% xlabel('Temperature [^{\circ}C]')
% ylabel('Heat Flow [W/g]')
% h = zeros(5, 1);
% h(1) = plot(NaN,NaN,'b', 'LineWidth', 1);
% h(2) = plot(NaN,NaN,'r', 'LineWidth', 1);
% h(3) = plot(NaN,NaN,'g', 'LineWidth', 1);
% h(4) = plot(NaN,NaN,'k', 'LineWidth', 1);
% h(5) = plot(NaN,NaN,'m', 'LineWidth', 1);
% legend(h, 'Unirradiated','10^{8} ions/cm^2','10^{9} ions/cm^2','10^{10} ions/cm^2','10^{11} ions/cm^2');

%% Explanatory figure
figure;
hold on
plot(T_corrected(:,2),Q_baseline(:,2),'k','LineWidth',1.5)
plot(T_corrected(:,58),Q_baseline(:,58),'k','LineWidth',1.5)
plot(T_corrected(:,2),Q_temp(:,2),'b','LineWidth',1.5)
plot(T_corrected(:,58),Q_temp(:,58),'r','LineWidth',1.5)

title('Representative DSC signal')
xlabel('Temperature [^{\circ}C]')
ylabel('Heat Flow [W/g], exo up')

h = zeros(2, 1);
h(1) = plot(NaN,NaN,'b', 'LineWidth', 1);
h(2) = plot(NaN,NaN,'r', 'LineWidth', 1);
legend(h,'Control','Irradiated')

%% Explanatory figure
figure;
hold on
% plot(T_corrected(:,2),Q_baseline(:,2),'k','LineWidth',1.5)
plot(T_corrected(:,2),Q_temp(:,2),'b','LineWidth',1.5)

title('Representative DSC signal')
xlabel('Temperature [^{\circ}C]')
ylabel('Heat Flow [W/g], exo up')

%% Spline baseline integrals

data = 68; 
Int = zeros(1,data);
for j = 1:data
    Int(j) = trapz(T_corrected(:,j),Q_corrected(:,j));
end

% Mean:
Int_mean = zeros(1,5);
Int_mean(1) = mean(Int(1:20));
Int_mean(2) = mean(Int(21:32));
Int_mean(3) = mean(Int(33:44));
Int_mean(4) = mean(Int(45:56));
Int_mean(5) = mean(Int(57:68));

% Std:
Int_std = zeros(1,5);
Int_std(1) = std(Int(1:20));
Int_std(2) = std(Int(21:32));
Int_std(3) = std(Int(33:44));
Int_std(4) = std(Int(45:56));
Int_std(5) = std(Int(57:68));

coord = [1 2 3 4 5];
figure;
errorbar(coord, Int_mean, Int_std,'o-','LineWidth',1)
title('Mean Peak Areas (heating)')
xlabel('Fluence [ions/cm^2]')
xticks([1 2 3 4 5]);
xticklabels({'Ctrl', '10^{8}','10^{9}','10^{10}','10^{11}'});
xlim([0.75, 5.25])
ax = gca;
ax.FontSize =10;
ylabel('Peak Area [W/g ^{\circ}C]')

%%
% coord = [1 2 3 4 5];
% figure(1)
% hold on
% plot(coord, -Int_std./Int_mean,'bo-')
% legend off
% h = zeros(1,7);
% h(1) = plot(NaN,NaN,'b', 'LineWidth', 1);
% h(2) = plot(NaN,NaN,'r', 'LineWidth', 1);
% h(3) = plot(NaN,NaN,'k', 'LineWidth', 1);
% h(4) = plot(NaN,NaN,'g', 'LineWidth', 1);
% h(5) = plot(NaN,NaN,'m', 'LineWidth', 1);
% h(6) = plot(NaN,NaN,'y', 'LineWidth', 1);
% h(7) = plot(NaN,NaN,'c', 'LineWidth', 1);
% legend(h, 'f=1.7','Original Spline','Spline mod f=2','f=3','f=2.5','f=4','f=1.5');
% 

%% Peak location baseline integrals

Peak = T_peak;

% Mean:
Peak_mean = zeros(1,5);
Peak_mean(1) = mean(Peak(1:20));
Peak_mean(2) = mean(Peak(21:32));
Peak_mean(3) = mean(Peak(33:44));
Peak_mean(4) = mean(Peak(45:56));
Peak_mean(5) = mean(Peak(57:68));

% Std:
Peak_std = zeros(1,5);
Peak_std(1) = std(Peak(1:20));
Peak_std(2) = std(Peak(21:32));
Peak_std(3) = std(Peak(33:44));
Peak_std(4) = std(Peak(45:56));
Peak_std(5) = std(Peak(57:68));

coord = [1 2 3 4 5];
figure;
errorbar(coord, Peak_mean, Peak_std,'o-','LineWidth',1)
title('Mean Peak Temperatures (heating)')
xlabel('Fluence [ions/cm^2]')
xticks([1 2 3 4 5]);
xticklabels({'Ctrl', '10^{8}','10^{9}','10^{10}','10^{11}'});
xlim([0.75, 5.25])
ax = gca;
c = ax.FontSize;
ax.FontSize =10;
ylabel('Peak Temperature [^{\circ}C]')

%% Plot results with error bands: patch

Int_plus1S = Int_mean+Int_std;
Int_minus1S = Int_mean-Int_std;
Int_plus2S = Int_mean+2*Int_std;
Int_minus2S = Int_mean-2*Int_std;
v1 = [coord, coord; Int_mean, Int_plus1S].';
v2 = [coord, coord; Int_mean, Int_minus1S].';
v3 = [coord, coord; Int_plus1S, Int_plus2S].';
v4 = [coord, coord; Int_minus1S, Int_minus2S].';
f = [1 6 7 2; 2 7 8 3; 3 8 9 4; 4 9 10 5];

figure;
hold on
patch('Faces',f,'Vertices', v1,'FaceColor',[153/255 204/255 255/255],'EdgeColor','none','FaceAlpha', 1)
patch('Faces',f,'Vertices', v2,'FaceColor',[153/255 204/255 255/255],'EdgeColor','none','FaceAlpha', 1)
patch('Faces',f,'Vertices', v3,'FaceColor',[204/255 229/255 255/255],'EdgeColor','none','FaceAlpha', 1)
patch('Faces',f,'Vertices', v4,'FaceColor',[204/255 229/255 255/255],'EdgeColor','none','FaceAlpha', 1)
plot(coord, Int_mean, 'ko-','LineWidth',2)

title('Melting Peak Area after Recrystallization at 10^{\circ}C/min')
xlabel('4.5 MeV \alpha Fluence [ions/cm^2]')
xticks([1 2 3 4 5]);
xticklabels({'Ctrl', '10^{8}','10^{9}','10^{10}','10^{11}'});
xlim([0.75, 5.25])
ax = gca;
ax.FontSize =10;
set(gcf,'Position',[360 100 700 525])
ylabel('Peak Area [W^{\circ}C/g]')

h = zeros(3, 1);
h(1) = plot(NaN,NaN,'ko-', 'LineWidth', 1);
h(2) = plot(NaN,NaN,'Color',[153/255 204/255 255/255], 'LineWidth', 1);
h(3) = plot(NaN,NaN,'Color',[204/255 229/255 255/255], 'LineWidth', 1);
legend(h, {'Mean','1\sigma', '2\sigma'},'Location', 'southwest',...
    'Orientation','horizontal','FontSize',12)




Peak_plus1S = Peak_mean+Peak_std;
Peak_minus1S = Peak_mean-Peak_std;
Peak_plus2S = Peak_mean+2*Peak_std;
Peak_minus2S = Peak_mean-2*Peak_std;
v5 = [coord, coord; Peak_mean, Peak_plus1S].';
v6 = [coord, coord; Peak_mean, Peak_minus1S].';
v7 = [coord, coord; Peak_plus1S, Peak_plus2S].';
v8 = [coord, coord; Peak_minus1S, Peak_minus2S].';
f = [1 6 7 2; 2 7 8 3; 3 8 9 4; 4 9 10 5];

figure;
hold on
patch('Faces',f,'Vertices', v5,'FaceColor',[153/255 204/255 255/255],'EdgeColor','none','FaceAlpha', 1)
patch('Faces',f,'Vertices', v6,'FaceColor',[153/255 204/255 255/255],'EdgeColor','none','FaceAlpha', 1)
patch('Faces',f,'Vertices', v7,'FaceColor',[204/255 229/255 255/255],'EdgeColor','none','FaceAlpha', 1)
patch('Faces',f,'Vertices', v8,'FaceColor',[204/255 229/255 255/255],'EdgeColor','none','FaceAlpha', 1)
plot(coord, Peak_mean, 'ko-','LineWidth',2)

title('Melting Peak Temperature after Recrystallization at 10^{\circ}C/min')
xlabel('4.5 MeV \alpha Fluence [ions/cm^2]')
xticks([1 2 3 4 5]);
xticklabels({'Ctrl', '10^{8}','10^{9}','10^{10}','10^{11}'});
xlim([0.75, 5.25])
ax = gca;
ax.FontSize =10;
set(gcf,'Position',[360 100 700 525])
ylabel('Peak Temperature [^{\circ}C/g]')

h = zeros(3, 1);
h(1) = plot(NaN,NaN,'ko-', 'LineWidth', 1);
h(2) = plot(NaN,NaN,'Color',[153/255 204/255 255/255], 'LineWidth', 1);
h(3) = plot(NaN,NaN,'Color',[204/255 229/255 255/255], 'LineWidth', 1);
legend(h, {'Mean','1\sigma', '2\sigma'},'Location', 'northwest',...
    'Orientation','horizontal','FontSize',12)

%% Check against XRD:

% Box number coordinating with DSC samples 1:20;1:12;1:12;1:12;1:12
boxNum = [1; 5; 9; 13; 17; 2; 6; 10; 14; 3; 7; 15; 19; 4; 23; 24; 20;
16; 12; 8; 1; 5; 13; 2; 6; 10; 12; 8; 4; 19; 11; 7; 1; 5; 13; 2; 14;
18; 16; 12; 8; 11; 7; 3; 1; 5; 9; 13; 17; 2; 4; 15; 11; 3; 14; 6; 1; 5;
9; 13; 6; 10; 12; 4; 11; 7; 3; 14];

label1 = {'ctrl 1','ctrl 2','ctrl 3'};
label2 = {'1e11 1', '1e11 2', '1e11 3'};

figure;
% subplot(2,1,1)
hold on
plot(boxNum(1:20),Int(1:20),'bo')
plot(boxNum([3,8,11]),Int([3,8,11]),'bo','LineWidth',2)
text(boxNum([3,8,11]),Int([3,8,11]),label1)
% title('control')

% subplot(2,1,2)
% hold on
plot(boxNum(57:68),Int(57:68),'mo')
plot(boxNum([59,62,66]),Int([59,62,66]),'mo','LineWidth',2)
text(boxNum([59,62,66]),Int([59,62,66]),label2)
% title('1e11')

XRD = [58.70212619, 50.02102675, 48.42222019,...
    62.32111333, 63.69459048, 61.70825028];
labelXRD = {'ctrl box 9','ctrl box 10','ctrl box 7', '1e11 box 9', '1e11 box 10', '1e11 box 7'};

figure;
plot(Int([3,8,11,59,62,66]),XRD,'o', 'LineWidth',2)
text(Int([3,8,11,59,62,66]),XRD,labelXRD,'VerticalAlignment','top')
xlabel('Peak Area (DSC)')
ylabel('Rel. % crystallinity (XRD)')
title('Comparison of XRD and DSC results')

%% More XRD

ctrl_box1 = [45.38 44.43];
ctrl_box3 = [53.87 60.21];
ctrl_box5 = [58.37];
ctrl_box6 = [54.04];
ctrl_box7 = [45.92 49.31];
ctrl_box9 = [56.53];
ctrl_box10 = [50.93];
ctrl_box15 = [62.76 67.79 58.08];
ctrl_box19 = [50.19 63.04];
e11_box1 = [62.58 59.71];
e11_box5 = [62.74 68.04 65.16];
e11_box6 = [65.48 73.83];
e11_box7 = [62.71 60.48];
e11_box9 = [59.22 61.41];
e11_box10 = [64.39];
e11_box11 = [69.12 63.61];
e11_box12 = [58.75];
e11_box13 = [52.32 67.71];

XRD_all = [ctrl_box1 ctrl_box3 ctrl_box5 ctrl_box6 ctrl_box7 ctrl_box9 ...
    ctrl_box10 ctrl_box15 ctrl_box19 e11_box1 e11_box5 e11_box6 e11_box7 ...
    e11_box9 e11_box10 e11_box11 e11_box12 e11_box13];
irrad_all = [ones(1,15),2*ones(1,17)];

inds_1 = [1,10,2,7,11,3,8,12,13,57,58,61,66,59,62,65,63,61];
inds_2 = [1,10,11,12,13,57,58,61,66,59,65,61];
inds_3 = [12,58];

XRD_1 = [ctrl_box1(1) ctrl_box3(1) ctrl_box5(1) ctrl_box6(1) ctrl_box7(1) ctrl_box9(1) ...
    ctrl_box10(1) ctrl_box15(1) ctrl_box19(1) e11_box1(1) e11_box5(1) e11_box6(1) e11_box7(1) ...
    e11_box9(1) e11_box10(1) e11_box11(1) e11_box12(1) e11_box13(1)];
XRD_2 = [ctrl_box1(2) ctrl_box3(2) ctrl_box7(2) ...
    ctrl_box15(2) ctrl_box19(2) e11_box1(2) e11_box5(2) e11_box6(2) e11_box7(2) ...
    e11_box9(2) e11_box11(2) e11_box13(2)];
XRD_3 = [ctrl_box15(3) e11_box5(3)];

irrad_1 = [ones(1,9),2*ones(1,9)];
irrad_2 = [ones(1,5),2*ones(1,7)];
irrad_3 = [ones(1,1),2*ones(1,1)];

labelXRD_1 = {'ctrl box1','ctrl box3','ctrl box5','ctrl box6','ctrl box7','ctrl box9',...
    'ctrl box10','ctrl box15','ctrl box19','1e11 box1', '1e11 box5', '1e11 box6'...
    ,'1e11 box7', '1e11 box9', '1e11 box10','1e11 box11', '1e11 box12', '1e11 box13'};
labelXRD_2 = {'ctrl box1','ctrl box3','ctrl box7',...
    'ctrl box15','ctrl box19','1e11 box1', '1e11 box5', '1e11 box6'...
    ,'1e11 box7', '1e11 box9', '1e11 box11', '1e11 box13'};
labelXRD_3 = {'ctrl box15','1e11 box5'};


mean_XRD_all = [mean([ctrl_box1 ctrl_box3 ctrl_box5 ctrl_box6 ctrl_box7 ctrl_box9 ...
    ctrl_box10 ctrl_box15 ctrl_box19]) mean([e11_box1 e11_box5 e11_box6 e11_box7 ...
    e11_box9 e11_box10 e11_box11 e11_box12 e11_box13])];
std_XRD_all = [std([ctrl_box1 ctrl_box3 ctrl_box5 ctrl_box6 ctrl_box7 ctrl_box9 ...
    ctrl_box10 ctrl_box15 ctrl_box19]) std([e11_box1 e11_box5 e11_box6 e11_box7 ...
    e11_box9 e11_box10 e11_box11 e11_box12 e11_box13])];

figure;
hold on
plot(irrad_1,XRD_1,'bo')
plot(irrad_2,XRD_2,'go')
plot(irrad_3,XRD_3,'ro')
errorbar([1,2], mean_XRD_all, std_XRD_all,'ko')
xlim([0,3])
xlabel('4.5 MeV \alpha Fluence [ions/cm^2]')
xticks([1 2]);
xticklabels({'Ctrl', '10^{11}'});
ylabel('Relative Crystallinity [%]')
title('All XRD results post-DSC')
h = zeros(4, 1);
h(1) = plot(NaN,NaN,'bo', 'LineWidth', 1);
h(2) = plot(NaN,NaN,'go', 'LineWidth', 1);
h(3) = plot(NaN,NaN,'ro', 'LineWidth', 1);
h(4) = plot(NaN,NaN,'ko', 'LineWidth', 1);
legend(h, '1st Measurement','2nd Measurement','3rd Measurement', 'Average')

%%

figure;
hold on
plot(Int(inds_1),XRD_1,'bo', 'LineWidth',2)
% text(Int(inds_1),XRD_1,labelXRD_1,'VerticalAlignment','top')
plot(Int(inds_2),XRD_2,'bo', 'LineWidth',2)
% text(Int(inds_2),XRD_2,labelXRD_2,'VerticalAlignment','top')
plot(Int(inds_3),XRD_3,'bo', 'LineWidth',2)
% text(Int(inds_3),XRD_3,labelXRD_3,'VerticalAlignment','top')
xlabel('Peak Area (DSC)')
ylabel('Rel. % crystallinity (XRD)')
title('Comparison of XRD and DSC results')
% h = zeros(3, 1);
% h(1) = plot(NaN,NaN,'bo', 'LineWidth', 1);
% h(2) = plot(NaN,NaN,'go', 'LineWidth', 1);
% h(3) = plot(NaN,NaN,'ro', 'LineWidth', 1);
% legend(h, '1st Measurement','2nd Measurement','3rd Measurement')















