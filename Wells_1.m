%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Cellpose for segmentation of well %%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies: (1) Medical Imaging Toolbox
%               (2) Cellpose Trained models (in the filepath)
%               (3) Deep Learning Toolbox
%               (4) Statistics and Machine Learning Toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathname1 = ('C:\Users\Syd_R\OneDrive\Desktop\RA position\hemolysin_assays\19_02_2024 - run1\'); %Folder where your images are saved
time_interval =2; % in minutes
total_frames =30;
coord_x =1900; 
coord_y =1700;  
ROI_width = 500;
ROI_height =500;
nBins=30;
number_of_gaussians_to_fit=2;
%% Single frame analysis
[fn, pn]=uigetfile('*.tiff','Load first frame');
vv1=[pn,fn];
x1=(imread(vv1));
x1=mat2gray(im2double(x1(:,:,1)));
I3=imcrop(x1,[coord_x coord_y ROI_width ROI_height]);
figure(1);imshow(I3); title('Grayscale ROI')
cp = cellpose(Model="cyto2");
averageCellDiameter = 15;
labelsDefault = segmentCells2D(cp,I3,ImageCellDiameter=averageCellDiameter);
loverlayDefault = labeloverlay(I3,labelsDefault);
figure(2);imshow(loverlayDefault); title('Cellpose Mask')

tic
L=logical(labelsDefault);
seg=L.*I3;
surf(seg); colormap jet; shading interp; axis tight; colorbar;
figure(3);imshow(seg); colormap jet; shading interp; axis tight; colorbar;  title('Cellpose Segmentation')
toc

s = regionprops(L, 'PixelIdxList','Centroid');
k = 1:numel(s);
len_k=length(k);
figure(4); imshow(loverlayDefault);title(['Total number of wells in ROI =', ' ',num2str(len_k)])
hold on
for k = 1:numel(s)
    c = s(k).Centroid;
    text(c(1), c(2), sprintf('%d', k), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','Color','w','FontSize',8);
end
hold off

%% Multi frame analysis (using the same segemented mask for multiple images)

for m=1:total_frames
    I=imread([pathname1,num2str(m),'.tiff']);  
    I2=(im2double(I(:,:,1)));
    I3=imcrop(I2,[coord_x coord_y ROI_width ROI_height]);
    seg=L.*I3;
    for ii=1:length(s)
        idx = s(ii).PixelIdxList;
        zz=seg(idx);
        int_vec(m,ii)=mean(zz(:));
    end
end

%% Intensity vs time plots for all wells

total_time =0:time_interval:(time_interval*m)-1;
total_time=total_time.*60;
cmap = jet(length(int_vec));
for kk=1:length(int_vec)
    figure(6);
    p=plot(total_time',int_vec(:,kk));
    p.Color=[cmap(kk,1) cmap(kk,2) cmap(kk,3)];
    hold on;
    xlabel('Time (sec)')
    ylabel('Normalized Intensity')
    title(['Intensity Vs. Time plot for ',num2str(len_k),' ','wells.'])
    fontsize(16,"points");
    axis tight;
end
hold off;

%% Curve fit and Plot all column vectors

total_time =0:time_interval:(time_interval*total_frames)-1;
total_time=total_time*60;
cmap = jet(length(int_vec));
for kk=1:length(int_vec)
    f= fit(total_time',int_vec(:,kk),'exp2');
    coef_c_1(:,kk)=f.a;
    coef_c_2(:,kk) =f.c*exp(f.d*total_time);
    coef_k(:,kk) = f.b;
    figure(7);
    p=plot(f);
    p.Color=[cmap(kk,1) cmap(kk,2) cmap(kk,3)];
    xlim([0 max(total_time)])
    ylim([0 max(int_vec(:))+0.1])
    hold on
    title('Exponential Curve Fits =')
    legend off;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Histogram of rate constant (k)
nbins=30;
figure(8); 
h=histogram(abs(coef_k),nbins);
xlabel('Rate constant K (sec.^{-1})')
ylabel('Frequency')
title(['Rate constant K for ',num2str(len_k),' ','wells.'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Histogram of rate constant + Sum of Gaussians
IK=abs(coef_k);
gmdist = fitgmdist(IK', number_of_gaussians_to_fit);
gmsigma = gmdist.Sigma;
gmmu = gmdist.mu;
gmwt = gmdist.ComponentProportion;

figure(8);
h=histogram(IK,nbins,'Normalization', 'pdf');
h.FaceColor = [1 0.5 0];
h.LineWidth =0.5;
h.EdgeColor = 'k';
x = 0:0.0001:round(max(IK)*3,2);
xlim([0 round(max(IK)*3,2)])
%ylim([0 200])
hold on;
plot(x, pdf(gmdist, x'), 'k')
xlabel('Rate constant K (sec.^{-1})')
ylabel('Frequency')
title(['Rate constant K for ',num2str(len_k),' ','wells.'])
hold off;

%% Drift Check
% ROI
[fn, pn]=uigetfile('*.tiff','Load First Frame for Drift Check');
vv1=[pn,fn];
d1=(imread(vv1));
d2=mat2gray(im2double(d1(:,:,1)));
d3_crop=imcrop(d2,[coord_x coord_y ROI_width ROI_height]);
figure (9), imshow(d3_crop); title('Grayscale ROI (Select ROI)')
h = imrect(gca,[10 10 100 100]); 
position = wait(h); % returns coordinates in "position" when user doubleclicks on rectangle
ROI_DC_XY=imcrop(d3_crop,position);
%figure(10);imshow(ROI_DC_XY); title('ROI of Grayscale ROI')

% ROI & ROI of Grayscale ROI
figure(10);
subplot(121); imshow(d3_crop); title('Grayscale ROI')
hold all;  % hold image and plot rectangle on top
rectangle('Position', position,'EdgeColor','y','LineWidth',2);
subplot(122); imshow(ROI_DC_XY); title('ROI of Grayscale ROI')


%%%% Load second image
[fn, pn]=uigetfile('*.tiff','Load second Frame for Drift Check');
vv1_2=[pn,fn];
d1_2=(imread(vv1_2));
d2_2=mat2gray(im2double(d1_2(:,:,1)));
d3_crop_2=imcrop(d2_2,[coord_x coord_y ROI_width ROI_height]);
ROI_DC_XY_2=imcrop(d3_crop_2,position);


% Using fft cross correlations to detect the moving distance[1] https://medium.com/@windless99/drift-correction-in-matlab-f9cc02860c09
% [sy,sx]=size(I1);
% fftim1=fft2(I1);
% fftim2=fft2(I2);
% cc=fftshift(ifft2(fftim1.*conj(fftim2)));
% [shiftY,shiftX]=find(cc==max(cc(:)));
% shiftY=shiftY-fix(sy/2)-1;
% shiftX=shiftX-fix(sx/2)-1;
% figure(2);
% imshow(mat2gray(cc)); hold on;
% plot(fix(sx/2),fix(sy/2),'r+'); hold off;
%%%%%
c = normxcorr2(ROI_DC_XY,ROI_DC_XY_2);
[max_c, imax] = max(abs(c(:))); %find the max value 
[ypeak, xpeak] = ind2sub(size(c),imax(1)); %Find peak in cross-correlation.
corr_offset = round([(xpeak-(size(c,2)+1)/2) (ypeak-(size(c,1)+1)/2)]);
offset = corr_offset;
xoffset = offset(1)
yoffset = offset(2)
Corr_coef=max_c

%%%%%%%%%%% Multiframe drift check
%%%%% Select a ROI from the first frame for drift check

for m=1:total_frames
    I_x=imread([pathname1,num2str(m),'.tiff']);  
    I2_x=(im2double(I_x(:,:,1)));
    I3_x=imcrop(I2_x,[coord_x coord_y ROI_width ROI_height]);
    if exist('position')==0
        imshow (I3_x);
        h = imrect(gca,[10 10 100 100]); 
        position = wait(h); % returns coordinates in "position" when user doubleclicks on rectangle
        Frame1=imcrop(I3_x,position);

        c = normxcorr2(Frame1,Frame1);
        [max_c, imax] = max(abs(c(:))); %find the max value
        [ypeak, xpeak] = ind2sub(size(c),imax(1)); %Find peak in cross-correlation.
        corr_offset = round([(xpeak-(size(c,2)+1)/2) (ypeak-(size(c,1)+1)/2)]);
        offset = corr_offset;
        xoffset(1,m) = offset(1);
        yoffset(1,m) = offset(2);
        corr_coef(1,m)=max_c;
        
    else
    imshow(I3_x)
    Frame2 = imcrop(I3_x,[position(:,1) position(:,2) position(:,3) position(:,4)]);
    c = normxcorr2(Frame1,Frame2);
    [max_c, imax] = max(abs(c(:))); %find the max value
    [ypeak, xpeak] = ind2sub(size(c),imax(1)); %Find peak in cross-correlation.
    corr_offset = round([(xpeak-(size(c,2)+1)/2) (ypeak-(size(c,1)+1)/2)]);
    offset = corr_offset;
    xoffset(1,m) = offset(1);
    yoffset(1,m) = offset(2);
    corr_coef(1,m)=max_c; 
    end  
     
end

figure(12)
plot(corr_coef,'LineWidth',2,'Color','r');
xlabel('Frames')
ylabel('Corr. coef.')
title('Corr. Coef. Vs. Frames')
fontsize(16,"points");
axis tight;

figure(13)
plot(xoffset,'LineWidth',2,'Color','g');hold on;
plot(yoffset, 'LineWidth',2,'Color','b');hold off;
xlabel('Frames')
ylabel('Xoffset, Yoffset (in Pixels)')
title('Xoffset, Yoffset Vs. Frames')
fontsize(16,"points");
legend ('Xoffset','Yoffset')
axis tight;

