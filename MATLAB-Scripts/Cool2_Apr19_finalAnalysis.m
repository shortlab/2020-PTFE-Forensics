%% Cool2_Apr19
% Includes spline baseline connecting linear sections,
% Peak area analysis figures 

%% Load data
% Control samples
% Editable filename...
path_main = 'C:\Users\Rachel\Documents\Rachel\MIT\MIT Research\Short Lab\';
path_specific = 'PTFE\Irradiated April 2019\CSV\';
s1 = 'ptfe_postirr_ctrl_';
s2 = '1';
s3 = '_cool2.csv';
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
    T(:,i) = flip(T_temp(1:x));
    Q(:,i) = flip(Q_temp(1:x)); 
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
    T(:,i) = flip(T_temp(1:x));
    Q(:,i) = flip(Q_temp(1:x)); 
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
    T(:,i) = flip(T_temp(1:x));
    Q(:,i) = flip(Q_temp(1:x)); 
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
    T(:,i) = flip(T_temp(1:x));
    Q(:,i) = flip(Q_temp(1:x)); 
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
    T(:,i) = flip(T_temp(1:x));
    Q(:,i) = flip(Q_temp(1:x)); 
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
left_end = 20; % degrees from peak
right_end = 15; % degrees from peak
f = 2;

% Setup correct ranges

[~,index] = max(Q_int(:,1));
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
    [~,index] = max(Q_int(:,j));
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

%% Explanatory figure
figure;
hold on
plot(T_corrected(:,2),Q_temp(:,2),'b','LineWidth',0.25)
plot(T_fit(:,2),Q_fit(:,2),'b.','LineWidth',0.75)
plot(T_corrected(:,2),Q_baseline(:,2),'b--','LineWidth',0.5)

title('Creating Baseline')
xlabel('Temperature [^{\circ}C]')
ylabel('Heat Flow [W/g]')

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
title('Mean Peak Areas (cooling)')
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
title('Mean Peak Temperatures (cooling)')
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

title('Crystallization Peak Area at 10^{\circ}C/min')
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
legend(h, {'Mean','1\sigma', '2\sigma'},'Location', 'northwest',...
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

title('Crystallization Peak Temperature at 10^{\circ}C/min')
xlabel('4.5 MeV \alpha Fluence [ions/cm^2]')
xticks([1 2 3 4 5]);
xticklabels({'Ctrl', '10^{8}','10^{9}','10^{10}','10^{11}'});
xlim([0.75, 5.25])
set(gca,'FontSize',10)
set(gcf,'Position',[360 100 700 525])
ylabel('Peak Temperature [^{\circ}C/g]')

h = zeros(3, 1);
h(1) = plot(NaN,NaN,'ko-', 'LineWidth', 1);
h(2) = plot(NaN,NaN,'Color',[153/255 204/255 255/255], 'LineWidth', 1);
h(3) = plot(NaN,NaN,'Color',[204/255 229/255 255/255], 'LineWidth', 1);
legend(h, {'Mean','1\sigma', '2\sigma'},'Location', 'northwest',...
    'Orientation','horizontal','FontSize',12)







