%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Nanowell Project: Image Analysis  %%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies: (1) Medical Imaging Toolbox
%               (2) Cellpose Trained models (in the filepath) Go to--> Get Add-ons--> Medical Imaging Toolbox Interface for Cellpose Library
%               (3) Deep Learning Toolbox
%               (4) Statistics and Machine Learning Toolbox
% $Author: Siddharth Rawat $    $Date: 2024/08/06 10:40:52 pm $    $Revision: 0.3 $
%                      $ Email: siddharth.rawat@unsw.edu.au
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathname1 = ('D:\Confluence_data\20240927\assay_488\New folder\'); %Folder where your images are saved
time_interval =1; % between image frames in mins.
total_frames =51; % total number of images in your folder 
coord_x =2828;  % initial x pixel location 1730
coord_y =2361;  % initial y pixel location 1810
ROI_width =190; % number of pixels to crop 
ROI_height =190; % number of pixels to crop 
format long;
format compact;
%% Single frame analysis
[fn, pn]=uigetfile('*.tiff','Load first frame'); % first frame 
vv1=[pn,fn];
x1=(imread(vv1)); % read first frame 
x1=mat2gray(im2double(x1(:,:,1))); % normalize
x1=medfilt2(x1);  % median filter
I2=imcrop(x1,[coord_x coord_y ROI_width ROI_height]); % crop the above frame
I3=padarray(I2,[20 20],0,'both');% zero padding

f = waitbar(0, 'Starting');

figure(1);imagesc(I3); title('Zero Padded ROI'); 
cp = cellpose(Model="nuclei"); % define the cellpose model of your choice e.g., cyto2 etc.
averageCellDiameter = 20; % in pixels
labelsDefault = segmentCells2D(cp,I3,ImageCellDiameter=averageCellDiameter); 
loverlayDefault = labeloverlay(I3,labelsDefault);
figure(2);imagesc(loverlayDefault); title('Cellpose Mask')

tic
L=logical(labelsDefault);
seg=L.*I3;
%seg=(seg-min(seg(:)))./(max(seg(:))-min(seg(:)));
figure(156);surf(seg); colormap jet; shading interp; axis tight;
figure(3);imshow(seg); colormap jet; shading interp; axis tight;  title('Cellpose Segmentation')
toc

s = regionprops(L, 'PixelIdxList'); % list containing raw pixel values for all segmented wells
s1 = regionprops(L,'Centroid'); % list containing centroid positions (x,y) for all segmented wells
k = 1:numel(s);
len_k=length(k);
figure(4); imshow(loverlayDefault./2);title(['Total number of wells in ROI =', ' ',num2str(len_k)])
hold on

for k = 1:numel(s)
    c = s1(k).Centroid;
    text(c(1), c(2), sprintf('%d', k), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','Color','w','FontSize',8);        
waitbar(k/numel(s), f, sprintf('Cellpose segmentation: %d %%', floor(k/numel(s)*100)));
end
hold off
close(f)
%% Create circular masks from centroids
f4 = waitbar(0,'Please wait...');

radius =8; % define the radius of well that you would like to crop
struct_centroid=struct2cell(s1'); % struct to cell
cell_centroid=cell2mat(struct_centroid); % centroids for all segmented wells
cell_centres(:,1)=cell_centroid(:,1,:); cell_centres(:,2)=cell_centroid(:,2,:);
radiis=(ones(1,length(cell_centres))).*radius; % vector
%%%% Double check if circles line up with image otherwise update cell centres
figure(31);
imagesc(I3)
h = viscircles(cell_centres,radiis);

%%%%% create a Mask of the region inside the circles 
wd= size(I3);
[columnsInImage, rowsInImage] = meshgrid(1:wd(1,2), 1:wd(1,1));

centerY=cell_centres(1:length(cell_centres),2);
centerX=cell_centres(1:length(cell_centres),1);
circlePixels = any((rowsInImage(:) - centerY').^2 ...
    + (columnsInImage(:) - centerX').^2 <= radius.^2, 2);
circlePixels = reshape(circlePixels, wd(1,1), wd(1,2));

% circlePixels is a 2D "logical" array.
% Now, display it.
figure(40);
image(circlePixels) ;
colormap([0 0 0; 1 1 1]);
hold on
scatter(centerX, centerY, 'r+');
hold off
title('Binary image of segmented wells');

figure(41);
cir_1=circlePixels.*I3.*2;
imagesc(cir_1) %% mask.* image

waitbar(.33,f4,'Generating Pixel List for well centres');

s2 = regionprops(circlePixels, 'PixelIdxList');
s22 =regionprops(circlePixels, 'PixelList');
title('Image of segmented wells');
%%%%% circles periphery: pixels values extraction
radius2=9; % radius of the bigger area= well area + extra area around well 
circlePixels1 = any((rowsInImage(:) - centerY').^2 ...
    + (columnsInImage(:) - centerX').^2 <= radius.^2, 2);
%circlePixels1 = reshape(circlePixels1, wd(1,1), wd(1,2));

circlePixels2 = any((rowsInImage(:) - centerY').^2 ...
    + (columnsInImage(:) - centerX').^2 <= radius2.^2, 2);
%circlePixels2 = reshape(circlePixels2, wd(1,1), wd(1,2));

circlePixels_peri = reshape(circlePixels2-circlePixels1, wd(1,1), wd(1,2));
circlePixels_peri=logical(circlePixels_peri);

waitbar(.67,f4,'Generating Pixel List for well periphery');

s3 = regionprops(circlePixels_peri, 'PixelIdxList'); % list containing raw pixel values for regions around all segmented wells
S4=regionprops(circlePixels_peri,'Centroid');
figure(43);
image(circlePixels_peri.*3) ;
colormap([0 0 0; 1 1 1]);
hold on
scatter(centerX, centerY, 'r+');
hold off
title('Binary image of the annular region around wells');

cir_2=circlePixels_peri.*I3.*2;
figure(45);
imagesc(cir_2)
title('Image of the annular region around wells');
close (f4)
%% Multi frame analysis (using seperate segemented mask for wells and local background subtraction for multiple images)
f5 = waitbar(0,'Please wait...');

v = VideoWriter('int.avi', 'MPEG-4');
v.FrameRate=1;
v.Quality=100;
open(v);
%int_vec=zeros(total_frames,length(s3)); %s2
for m=1:total_frames

    I=imread([pathname1,num2str(m),'.tiff']);  
    I2=(im2double(I(:,:,1))); 
    I2_1=imcrop(I2,[coord_x coord_y ROI_width ROI_height]);
    background = imopen(I2_1, strel('disk', 10));
    J1 = imsubtract(I2_1, background);
    
    I3=padarray(J1,[20 20],0,'both');% zero padding
    
   
    seg_mult=circlePixels.*I3;
    figure(45);imagesc(seg_mult); colormap jet; colorbar
    %pause(1)
 
    frame = getframe(gcf);
    writeVideo(v,frame);
    
    for ii=1:length(s3)  %% s2
        idx = s2(ii).PixelIdxList;
        zz=seg_mult(idx);
        idx_peri = s3(ii).PixelIdxList; 
        zz_peri=seg_mult(idx_peri);
        %int_vec(m,ii)=mean(zz(:))-mean(zz_peri(:));
        int_vec(m,ii)=mean(zz(:));
    end
    hold on;
    waitbar(m/total_frames, f5, sprintf('Multi Frame Analysis: %d %%', floor(m/total_frames*100)));

end
close (v)
close(f5)

%% Normalized intensity (plot starts with 1 then decays)
f1 = waitbar(0, 'Starting');
az=size(int_vec);
final_mat=zeros(total_frames,length(s2));%s2
for sdl=1:az(1,2)
    ax=max(int_vec(:,sdl));
    ax1=ax(:);
    final_mat(:,sdl)=int_vec(:,sdl)./ax;
end
%%%%% Plot
total_time =0:time_interval:(time_interval*m)-1;
total_time=total_time.*60;
jj1=size(final_mat); 
cmap = jet(jj1(1,2));
for kk1=1:jj1(1,2)
    figure(6);
    p=plot(total_time',final_mat(:,kk1),'-'); %--o
    p.Color=[cmap(kk1,1) cmap(kk1,2) cmap(kk1,3)];
    hold on;
    xlabel('Time (sec)')
    ylabel('Intensity (AU)')
    title(['Normalized Intensity Vs. Time plot for ',num2str(len_k),' ','wells.'])
    fontsize(16,"points");
    axis tight;
    waitbar(kk1/jj1(1,2), f1, sprintf('Normalized Intensity vs. time plots: %d %%', floor(kk1/jj1(1,2)*100)));
end
hold off;
close(f1)
%% Seperate normalized sub-plots


% [num_rows, num_cols] = size(final_mat);
% 
% % Generate x values (time steps)
% x = 1:num_rows;
% figure;
% 
% % Loop over each sample to create the stacked plots
% hold on;  % Hold on to plot multiple lines
% for col = 1:num_cols
%     y = final_mat(:, col);            % Get data for the current sample (column)
%     z = col * ones(size(x));     % Create a constant z value for each time step, stacking along the z-axis
%     plot3(x, y, z, 'LineWidth', 1.5);  % Plot the 2D data as a 3D line with thicker lines for visibility
% end

f23 = waitbar(0, 'Starting');
az=size(int_vec);
final_mat=zeros(total_frames,length(s2));%s2
for sdl=1:az(1,2)
    ax=max(int_vec(:,sdl));
    ax1=ax(:);
    final_mat(:,sdl)=int_vec(:,sdl)./ax;
end

[num_rows, num_cols] = size(final_mat);

% Number of subplots (one per column)
num_subplots = num_cols;
figure;
% Loop over each column to create a subplot
for col = 1:num_cols
    % Create subplot
    subplot(ceil(sqrt(num_subplots)), ceil(sqrt(num_subplots)), col);
    
    % Extract data for the current column
    y = final_mat(:, col);  % Data for the current column
    x = 1:num_rows;    % Time steps (or row indices)

    % Plot the data
    plot(x, y, 'LineWidth', 1.5);
    
    % Set labels and title for each subplot
    % xlabel('Time Steps');
    % ylabel('Data Value');
    % title(sprintf('Column %d', col));
    
    % Add grid for better readability
    grid on;
end
% Adjust figure layout
sgtitle('Normalized Intensity Vs. Time Plots of individual wells ');  % Super title for the whole figure
close(f23)

%% Normalized plot for a SIGNLE WELL
well_number=12;
az=size(int_vec);
xx=1:az(1,1);
plot(xx,final_mat(:,well_number),'--gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])

xlabel('Time (min.)')
ylabel('Normalized Intensity (au)')
title('Normalized Intensity vs. Time Plot')


%% Curve_fits1 (a*exp(-bx)) (OPTION 1)
f2 = waitbar(0, 'Starting');
total_time =0:time_interval:(time_interval*total_frames)-1;
total_time=total_time*60;
cmap = jet(length(final_mat));
jj=size(final_mat); 
for kk=1:jj(1,2)
    %f= fit(total_time',final_mat(:,kk),'exp1');
    f= fit(total_time',final_mat(:,kk),'exp1',Algorithm="Levenberg-Marquardt"); %% Algorithm="Levenberg-Marquardt"
    coef_c_1(:,kk)=f.a;
    coef_k1(:,kk) = f.b;
    %coef_c_2(:,kk)=f.c;
    %coef_k2(:,kk) = f.d;
    figure(7);
    p=plot(f); 
    p.Color=[cmap(kk,1) cmap(kk,2) cmap(kk,3)];  
    hold on;
    %k=plot(total_time',final_mat(:,kk),'o'); 
    %k.Color=[cmap(kk,1) cmap(kk,2) cmap(kk,3)];
    hold on;
    xlabel('Time (sec)')
    ylabel('Intensity (AU)')
    title(['Exponential Curve Fits for ',num2str(len_k),' ','wells.'])
   legend off;
   fontsize(16,"points");
   waitbar(kk/jj(1,2), f2, sprintf('exp1 curve fits (a*exp(-bx)): %d %%', floor(kk/jj(1,2)*100)));
end
close(f2)
%% Curve_fits2 (a*exp(-bx)+c)  (OPTION 2)
%%%% https://au.mathworks.com/matlabcentral/answers/592609-fit-curve-when-corresponding-to-one-x-value-we-have-multiple-y-value
%%%% https://au.mathworks.com/matlabcentral/answers/737197-exponential-decay-rate-constant
f6 = waitbar(0, 'Starting');
total_time =0:time_interval:(time_interval*total_frames)-1;
total_time=total_time*60;
tt=total_time';
cmap = jet(length(final_mat));
jj=size(final_mat); 
for kk=1:jj(1,2)
    yc=(final_mat(:,kk));
    tbl = table(tt(:), yc(:));
    modelfun = @(b,x) b(1) * exp(-b(2)*x(:, 1)) + b(3);  
    aGuessed = (mean(final_mat(:,kk))-min(final_mat(:,kk))); % Arbitrary sample values I picked.
    bGuessed = 0;
    cGuessed = min(final_mat(:,kk));
    beta0 = [aGuessed, bGuessed, cGuessed]; % Guess values to start with.  Just make your best guess.
    mdl = fitnlm(tbl, modelfun, beta0);
    coefficients = mdl.Coefficients{:, 'Estimate'};
    coef_a(:,kk)=coefficients(1);
    coef_b(:,kk)=coefficients(2);
    coef_c(:,kk)=coefficients(3);
    yFitted = coefficients(1) * exp(-coefficients(2)*tt) + coefficients(3);
    
    figure(7);
    p=plot(tt, yc, 'b.', 'MarkerSize', 5);
    p.Color=[cmap(kk,1) cmap(kk,2) cmap(kk,3)]; 
    grid on;
    k=plot(tt,yFitted, 'r-', 'LineWidth', 1); 
    k.Color=[cmap(kk,1) cmap(kk,2) cmap(kk,3)];  
    hold on;
    xlabel('Time (sec)')
    ylabel('Intensity (AU)')
    title(['Exponential Curve Fits for ',num2str(len_k),' ','wells.'])
    legend off;
    fontsize(16,"points");
    waitbar(kk/jj(1,2), f6, sprintf('exp1 curve fits (a*exp(-bx)+c) : %d %%', floor(kk/jj(1,2)*100)));
end
close(f6)
coef_k1=coef_b;
%% Curve_fits3 (a exp(-bx)+c)  (OPTION 3)
coef_b_star=zeros(1,length(final_mat));
f7 = waitbar(0, 'Starting');
total_time =0:time_interval:(time_interval*total_frames)-1;
total_time=total_time*60;
x=total_time;
cmap = jet(length(final_mat));
jj=size(final_mat); 
for kk=1:jj(1,2)
    y=final_mat(:,kk)';
    aGuessed = mean(y)-min(y); %sample values I picked.
    bGuessed = 0;
    cGuessed = min(y);
    f = @(b,x) b(1).*exp(-b(2).*x)+b(3);    % Objective Function
    B = fminsearch(@(b) norm(y - f(b,x)), [aGuessed; bGuessed; cGuessed]); % estimation
    coef_b_star(:,kk)=B(2,1);
    figure(7);
    p=plot(x, y, 'b.', 'MarkerSize', 5);
    p.Color=[cmap(kk,1) cmap(kk,2) cmap(kk,3)];  
    hold on
    k=plot(x, f(B,x), 'r-', 'LineWidth', 1); 
    k.Color=[cmap(kk,1) cmap(kk,2) cmap(kk,3)];  
    grid
   xlabel('Time (sec)')
   ylabel('Intensity (AU)')
   title(['Exponential Curve Fits for ',num2str(len_k),' ','wells.'])
   fontsize(16,"points");
   waitbar(kk/jj(1,2), f7, sprintf('exp1 curve fits (a*exp(-bx)+c): %d %%', floor(kk/jj(1,2)*100)));
end
close(f7)
coef_k1=coef_b_star;
%% S3 well locations

s3_1 = regionprops(circlePixels_peri,'Centroid');
kx = 1:numel(s3_1);
len_kx=length(kx);
figure(41); imshow(loverlayDefault);title(['Total number of wells in S3 ROI =', ' ',num2str(len_kx)])
hold on
for kx = 1:numel(s3)
    c1 = s3_1(kx).Centroid;
    text(c1(1), c1(2), sprintf('%d', kx), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','Color','w','FontSize',8);
end
hold off
%% Check which wells shows change (reaction rates vs. well number)
wd=1:length(coef_k1);
figure(42);
plot(wd,abs(coef_k1),'+r','LineWidth',2,'Color','r')
grid on;
hold on;
plot(wd,abs(coef_k1),'Color','r'); hold on;
legend off; axis tight;
xlabel('Well number')
ylabel('Reaction rate (K)')
title('Reaction rates Vs. Well number')


% for sd=1:length(coef_k1)
%     ddl(sd)=mean(final_mat(:,sd));
% end
% 
% plot(ddl,abs(coef_k1),'+r','LineWidth',2,'Color','r')
% grid on;
% hold on;
% plot(ddl,abs(coef_k1),'Color','r'); hold on;
% 
% 
% 
% Field_Size = wd;%rand(1,8); 
% tau        = (coef_k1);%rand(1,8);
% Field      = ddl; %1:8
% color      = jet(length(Field));
% markers    = 's'; %s^phvdxo
% sz         = 5:5:10; %6:2:20
% gsh = gscatter(Field_Size,tau,Field,color,markers,sz,'on','Field Size (TSCF)','\tau^{*} (Year)');
% % fill markers 
% for g = 1:length(gsh)
%     gsh(g).MarkerFaceColor = gsh(g).Color;
% end
% legend off;
% colormap jet

% plot(ddl,abs(coef_k1))
%% Simple Histogram of rate constant (k)
nbins=10;
figure(8); 
h=histogram(abs(coef_k1),nbins); 
h(1).FaceColor = [1 0.5 0];
h(1).LineWidth =0.3;
h(1).EdgeColor = 'g';
%h(2).Color = [0 0 0];
xlabel('Rate constant K (sec.^{-1})')
ylabel('Frequency')
title(['Rate constant K for ',num2str(len_k),' ','wells.'])
%set(gca,'YScale','log')

%% thresholded histogram
nbins= 30;
data=abs(coef_k1);
class1 = data(data <= 0.0025);
figure(8); 
h=histfit(class1,nbins,'kernel');
h(1).FaceColor = [1 0.5 0];
h(1).LineWidth =0.3;
h(1).EdgeColor = 'g';
h(2).Color = [0 0 0];
xlim([0 0.0025])
xlabel('Rate constant K (sec.^{-1})')
ylabel('Frequency')
title(['Rate constant K for ',num2str(length(class1)),' ','wells.'])

wd=1:length(class1);
figure(42);
plot(wd,class1,'+r','LineWidth',2,'Color','r')
grid on;
hold on;
plot(wd,class1,'Color','r'); hold on;
legend off; axis tight;
xlabel('Well number')
ylabel('Reaction rate (K)')
title('Reaction rates Vs. Well number')

%% Sum of Gaussian of above Histogram for rate constant (k) calculation
fac=length(coef_k1);
%nbins=round((nthroot(fac,3))*2); % Rice's rule https://www.statisticshowto.com/choose-bin-sizes-statistics/
nbins=30;
figure(8); 
h=histfit(abs(coef_k1),nbins,'kernel');
h(1).FaceColor = [1 0.5 0];
h(1).LineWidth =0.3;
h(1).EdgeColor = 'g';
h(2).Color = [0 0 0];
%ylim([0 100])
%xlim([0 0.01]) 
xlabel('Rate constant K (sec.^{-1})')
ylabel('Frequency')
title(['Rate constant K for ',num2str(len_k),' ','wells.'])
%% Choose 4 wells each with either 0, 1, 2, 3 alpha-Hla molecules and fit a function
%%%%%% Define parameters
Len_alpha=10*10^-9; % length of alpha-Hla in meters
N0_alpha=0.00001;% number of alpha-Hla pore: 0
N1_alpha=1;% number of alpha-Hla pore: 1
N2_alpha=2;% number of alpha-Hla pores: 2
N3_alpha=3;% number of alpha-Hla pores 3
N4_alpha=4;% number of alpha-Hla pores: 4
N5_alpha=5;% number of alpha-Hla pores: 5
N6_alpha=6;% number of alpha-Hla pores: 6
d=1*10^-9; % diameter of alpha-Hla pore in meters 
D=4*10^(-6); % diameter of well 
V=pi*((D/2)^2)*650*10^(-9); % volume of a single well in cubic meters
                            % volume in cubic meter (m^3) (1 m^3= 1000 litre)
xx= total_time'; % vector

%% Now choose the wells that you want to group together 
%%%% https://www.delftstack.com/howto/matlab/iterate-through-matrix-matlab/
%%%% wells with zero alpha-Hla pores
D0 = [167 275 799 1095]; % vector containing well numbers, update as per your need
[rows, columns] = size(D0);  %29 39 116 148 476 709
yy_0=zeros(az(1,1),length(D0)); %preallocation
for i = 1:rows
    for j = 1:columns
        yy_0(:,j) = final_mat(:,D0(i,j)); %array of wells with zero alpha-hla pores
        %disp(yy_test);
    end
end
% 125ug/mL 33 113 131 505 218 251 336 566 605 654

%%% wells with 1 alpha-Hla pores
D1 = [356 360 691 902]; % vector containing well numbers, update as per your need
[rows, columns] = size(D1); %198 224 229 235 268 304
yy_1=zeros(az(1,1),length(D1)); %preallocation
for ii = 1:rows
    for jj = 1:columns
        yy_1(:,jj) = final_mat(:,D1(ii,jj)); %array of wells with 1 alpha-hla pores
    end
end
% 125ug/mL 90 109 155 315 388 472 169 233 324 410

%%% wells with 2 alpha-Hla pores
D2 = [414 498 814 1237 ]; % vector containing well numbers, update as per your need
[rows, columns] = size(D0);
yy_2=zeros(az(1,1),length(D2)); %preallocation
for i = 1:rows
    for j = 1:columns
        yy_2(:,j) = final_mat(:,D2(i,j)); %array of wells with 1 alpha-hla pores
    end
end
% 18 269 398 446 507 579 686 292 476 477
%% Mean and STD for plotting errorbars 
%yy_mean_00=zeros(1,az(1,1)); % preallocation
%err_00=zeros(1,az(1,1)); % preallocation
for ze=1:az(1,1)
    yk_0=yy_0(ze,:);
    yy_mean_00(ze)=mean(yk_0(:));
    err_00(ze)=std(yy_0(ze,:));
end
yy_mean_0=yy_mean_00'; % vector
err_0=err_00';

%yy_mean_11=zeros(1,az(1,1)); % preallocation
%err_11=zeros(1,az(1,1)); % preallocation
for ze=1:az(1,1)
    yk_1=yy_1(ze,:);
    yy_mean_11(ze)=mean(yk_1(:));
    err_11(ze)=std(yy_1(ze,:));
end
yy_mean_1=yy_mean_11'; % vector
err_1=err_11';

yy_mean_22=zeros(1,az(1,1)); % preallocation
err_22=zeros(1,az(1,1)); % preallocation
for ze=1:jj(1,1)
    yk_2=yy_2(ze,:);
    yy_mean_22(ze)=mean(yk_2(:));
    err_22(ze)=std(yy_2(ze,:));
end
yy_mean_2=yy_mean_22'; % vector
err_2=err_22';

%% exponential fits 
%%%% Curve fit METHOD 1
%fk_0 = fit(xx,yy_mean_0,'exp1',Algorithm="Levenberg-Marquardt");
%fk_1 = fit(xx,yy_mean_1,'exp1',Algorithm="Levenberg-Marquardt");
%fk_2 = fit(xx,yy_mean_2,'exp1',Algorithm="Levenberg-Marquardt");

%%%% Curve fit METHOD 3
a0Guessed = mean(yy_mean_0)-min(yy_mean_0); %sample values I picked.
b0Guessed = 0.001;
c0Guessed = min(yy_mean_0);
f_0 = @(b0,xx) b0(1).*exp(-b0(2).*xx)+b0(3);    % Objective Function
B_0 = fminsearch(@(b0) norm(yy_mean_0 - f_0(b0,xx)), [a0Guessed; b0Guessed; c0Guessed]); % estimation
coef_f_0=B_0(2,1)

a1Guessed = mean(yy_mean_1)-min(yy_mean_1); %sample values I picked.
b1Guessed = 0;
c1Guessed = min(yy_mean_1);
f_1 = @(b1,x) b1(1).*exp(-b1(2).*xx)+b1(3);    % Objective Function
B_1 = fminsearch(@(b1) norm(yy_mean_1 - f_1(b1,xx)), [a1Guessed; b1Guessed; c1Guessed]); % estimation
coef_f_1=B_1(2,1) 



a2Guessed = mean(yy_mean_2)-min(yy_mean_2); %sample values I picked.
b2Guessed = 0;
c2Guessed = min(yy_mean_2);
f_2 = @(b2,x) b2(1).*exp(-b2(2).*xx)+b2(3);    % Objective Function
B_2 = fminsearch(@(b2) norm(yy_mean_2 - f_2(b2,xx)), [a2Guessed; b2Guessed; c2Guessed]); % estimation
coef_f_2=B_2(2,1)

%% coefs. of exponential fits

%%%% Coef METHOD 1
% coef_b_0=fk_0.b;
% coef_b_1=fk_1.b;
% coef_b_2=fk_2.b;

%%%% Coef METHOD 3
coef_b_0=coef_f_0;
coef_b_1=coef_f_1;
coef_b_2=coef_f_2;

%% Now extract diffusion coefficient from the above curve fits
num1_0=(4*Len_alpha*V)/(N0_alpha*d*d*pi);
D_Hla_0=coef_b_0*num1_0 % in m^2/sec, Diffusion coeffecient for 0 alpha-Hla pore

num1_1=(4*Len_alpha*V)/(N1_alpha*d*d*pi);
D_Hla_1=coef_b_1*num1_1 % in m^2/sec, Diffusion coeffecient for 1 alpha-Hla pore

num1_2=(4*Len_alpha*V)/(N2_alpha*d*d*pi);
D_Hla_2=coef_b_2*num1_2 % in m^2/sec, Diffusion coeffecient for 2 alpha-Hla pore

%% plot for the fitted curves
%%% https://kakearney.github.io/2016/06/10/boundedline.html
%%% Patch 1 for 0 alpha-Hla
% lo_0=yy_mean_0-err_0;
% hi_0=yy_mean_0+err_0;
% hp_0=patch([xx;xx(end:-1:1)],[lo_0;hi_0(end:-1:1)],'r');

%%%% Patch 2 for 1 alpha-Hla
lo_1=yy_mean_1-err_1;
hi_1=yy_mean_1+err_1;
hp_1=patch([xx;xx(end:-1:1)],[lo_1;hi_1(end:-1:1)],'g');

% %%%% Patch 3 for 2 alpha-Hla
% lo_2=yy_mean_2-err_2;
% hi_2=yy_mean_2+err_2;
% hp_2=patch([xx;xx(end:-1:1)],[lo_2;hi_2(end:-1:1)],'b');

%%% 0 alpha-Hla plots: errorbar + curve fits
% hold on;
% hl_0=line(xx,yy_mean_0);
% set(hp_0, 'facecolor', [1 0.8 0.8], 'edgecolor', 'none');
% set(hl_0, 'color', 'r');
% hold on; 
% errorbar(xx,yy_mean_0,err_0,'',"MarkerSize",5,...
%     "MarkerEdgeColor","red","MarkerFaceColor",[1 0 0]); hold on;
% p0 = plot(xx,f_0(B_0,xx),'k'); hold on; %=plot(fk_0,'k');
% set(p0,'lineWidth',1.5); 

%%% 1 alpha-Hla plots: errorbar + curve fits
hold on;
hl_1=line(xx,yy_mean_1);
set(hp_1, 'facecolor', [0.8 1 0.8], 'edgecolor', 'none');
%set(hl_1, 'color', 'r', 'marker', 'x');
hold on;
errorbar(xx,yy_mean_1,err_1,'',"MarkerSize",5,...
    "MarkerEdgeColor","green","MarkerFaceColor",[0 1 0]); hold on;
p1 = plot(xx,f_1(B_1,xx),'k'); hold on; %=plot(fk_1,'k'); 
set(p1,'lineWidth',1.5); 

%%% 2 alpha-Hla plots: errorbar + curve fits
% hold on;
% hl_2=line(xx,yy_mean_2);
% set(hp_2, 'facecolor', [0.8 0.8 1], 'edgecolor', 'none');
% %set(hl_2, 'color', 'r', 'marker', 'x');
% hold on;
% errorbar(xx,yy_mean_2,err_2,'',"MarkerSize",5,...
%     "MarkerEdgeColor","red","MarkerFaceColor",[0 0 1]); hold on;
% p2 = plot(xx,f_2(B_2,xx),'k'); hold on; %=plot(fk_1,'k'); 
% set(p2,'lineWidth',1.5); 
axis tight

%%%% Figure properties
ylim([0 1])
xlabel('Time (sec.)')
ylabel('Normalized Intensity (AU)')
title('Exponential Regression: Curve fits for Wells With 0,1,2 \alpha-Hla Pores')
txt = ['D_1{\alpha-Hla} ='  num2str(D_Hla_1),'m^2 s^{-1}.'];
text(50,0.3,txt)
% txt1 = ['D_2{\alpha-Hla} ='  num2str(D_Hla_2),'m^2 s^{-1}.'];
% text(50,0.2,txt1)
% txt1 = ['D_2{\alpha-Hla} ='  num2str(D_Hla_2),'m^2 s^{-1}'];
% text(50,0.2,txt1)
txt2 = ('Diffusion Coef. for 1 \alpha-Hla Pores');
text(50,0.4,txt2,'FontSize', 12, 'FontWeight', 'bold')
legend off; 
