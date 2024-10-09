%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Nanowell Project: Image Analysis  %%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies: (1) Medical Imaging Toolbox
%               (2) Cellpose Trained models (in the filepath) Go to--> Get Add-ons--> Medical Imaging Toolbox Interface for Cellpose Library
%               (3) Deep Learning Toolbox
%               (4) Statistics and Machine Learning Toolbox
% $Author: Siddharth Rawat $    $Date: 2024/10/02 10:33:26 pm $    $Revision: 0.1 $
%                      $ Email: siddharth.rawat@unsw.edu.au
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Automatically detect 52 OECTS using Normalized cross-correlation.
% Either use cellpose segmentation (step 4.1) or Gaussian fit (step 4.2) to find well intensity 
%% STEP 1: Parameters
pathname1 = ('C:\Users\Syd_R\OneDrive\Desktop\20240927\20240927\assay_488\New folder\'); %Folder where your images are saved
time_interval =1; % between image frames in mins.
total_frames =25; % total number of images in your folder 
%% STEP 2: Crop a single device from the 52 locations (Template)
% Step 1: Read and display the image
[fn, pn]=uigetfile('*.tiff','Load First Image');
vv1=[pn,fn];
X=(imread(vv1));
main_img=X(:,:,1); main_img2=imcrop(main_img, [1500 1200 1200 1200]);
figure(1);imshow(main_img2.*5)
title('Select two points to crop the image');
% Step 2: Select two points using ginput
[x, y] = ginput(2);  % User selects two points with the mouse

% Ensure points are integers
x = round(x);
y = round(y);

% Step 3: Calculate the cropping rectangle
% Top-left corner of the crop
x1 = min(x);
y1 = min(y);

% Width and height of the crop
width = abs(x(2) - x(1));
height = abs(y(2) - y(1));

% Crop the image using the calculated rectangle
croppedImage = imcrop(main_img2, [x1, y1, width, height]);

% Step 4: Display the cropped image
figure(2);
imshow(croppedImage.*5);
title('Cropped Template Image');
%% STEP 3: Detect 52 device locations automatically using normalized cross-correlation
% Read the main image and the template image
mainImage = main_img.*5;    % Replace with your main image file
templateImage = croppedImage.*5;  % Replace with your template image file

% Convert images to grayscale if they are RGB
if size(mainImage, 3) == 3
    mainImageGray = rgb2gray(mainImage);
else
    mainImageGray = mainImage;
end

if size(templateImage, 3) == 3
    templateImageGray = rgb2gray(templateImage);
else
    templateImageGray = templateImage;
end

% Get the dimensions of the template
[templateHeight, templateWidth] = size(templateImageGray);
[imageHeight, imageWidth] = size(mainImageGray);

% Perform normalized cross-correlation
crossCorrelationOutput = normxcorr2(templateImageGray, mainImageGray);

% Crop cross-correlation output to match original image size
crossCorrelationOutput = crossCorrelationOutput(templateHeight:end-templateHeight+1, ...
                                                templateWidth:end-templateWidth+1);

% Define a threshold for detection
threshold = 0.8;  % Adjust the threshold as needed for your data

% Find peaks above the threshold
detectionMask = crossCorrelationOutput > threshold;

% Apply non-maximum suppression
nonMaxSuppressed = imregionalmax(crossCorrelationOutput) & detectionMask;

% Get coordinates of detected peaks
[ypeaks, xpeaks] = find(nonMaxSuppressed);

% Display the result
figure(3);
imshow(mainImage);
hold on; 

% Draw rectangles around all detected regions
for i = 1:length(xpeaks)
    rectangle('Position', [xpeaks(i), ypeaks(i), templateWidth, templateHeight], ...
              'EdgeColor', 'r', 'LineWidth', 2);
    % % Display text labels with x and y coordinates
    % text(xpeaks(i), ypeaks(i) - 10, sprintf('(%d, %d)', xpeaks(i), ypeaks(i)), ...
    %      'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
    % Display sequential number labels
    text(xpeaks(i), ypeaks(i) - 10, sprintf('%d', i), ...
         'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
end

title('Detected Patterns Using Template Matching with Non-Maximum Suppression');
hold off;

% Display the result
figure(4);
imshow(mainImage);
hold on;

% Draw rectangles around all detected regions
for i = 1:length(xpeaks)
    rectangle('Position', [xpeaks(i), ypeaks(i), templateWidth, templateHeight], ...
              'EdgeColor', 'r', 'LineWidth', 2);
    % Display text labels with x and y coordinates
    text(xpeaks(i), ypeaks(i) - 10, sprintf('(%d, %d)', xpeaks(i), ypeaks(i)), ...
         'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
end

title('Detected Patterns Using Template Matching with Non-Maximum Suppression');
hold off;


%% STEP 4.1: Extract time series data from the 52 device locations (OPTION 1: Segmentation based)

f = waitbar(0,'Please wait...');
for m=1:total_frames

    I=imread([pathname1,num2str(m),'.tiff']);  
    I2=(im2double(I(:,:,1)));
    background = imopen(I2, strel('disk', 10));
    J1 = imsubtract(I2, background);
    for ii=1:length(xpeaks)
        I3=imcrop(J1,[xpeaks(ii)+15, ypeaks(ii)+35 40 40]);   
        
        back=imcrop(J1,[xpeaks(ii), ypeaks(ii)+15 15 50]);

       cp = cellpose(Model="nuclei"); % define the cellpose model of your choice e.g., cyto2 etc.
       averageCellDiameter = 20; % in pixels
       labelsDefault = segmentCells2D(cp,I3,ImageCellDiameter=averageCellDiameter); 
       loverlayDefault = labeloverlay(I3,labelsDefault);
       %figure(2);imagesc(loverlayDefault); title('Cellpose Mask')
       tic
       L=logical(labelsDefault);
       seg=L.*I3;
       figure(3);imshow(seg); colormap jet; shading interp; axis tight;  title('Cellpose Segmentation')
       toc 
        
       a=(mean(seg(:)));b=mean(back(:));
       %figure(5); imshow(mat2gray(I3));
       int_vec(m,ii)=(a-b);
    
    end
    
    waitbar(m/total_frames, f, sprintf('Time Series Analysis: %d %%', floor(m/total_frames*100)));
end
close(f)
%% STEP 4.2: Extract time series data from the 52 device locations (OPTION 2: Gaussian curve fits based)

f = waitbar(0,'Please wait...');
for m=1:total_frames

    I=imread([pathname1,num2str(m),'.tiff']);  
    I2=(im2double(I(:,:,1)));
    I2=adapthisteq(I2);
    %background = imopen(I2, strel('disk', 10));
    %J1 = imsubtract(I2, background);
    for ii=1:length(xpeaks)
        I3=imcrop(I2,[xpeaks(ii)+15, ypeaks(ii)+35 40 40]); 
        data=I3;
        [X, Y] = meshgrid(1:size(data, 2), 1:size(data, 1));

        % Reshape data to use for fitting
        xData = [X(:), Y(:)];
        zData = data(:);  % Your 2D spot data reshaped into a vector

        % Initial guess for the parameters [A, x0, y0, sigma_x, sigma_y, B]
        A_initial = max(data(:));         % Peak amplitude (initial guess)
        x0_initial = size(data, 2) / 2;   % Initial guess for x center
        y0_initial = size(data, 1) / 2;   % Initial guess for y center
        sigma_x_initial = 10;             % Initial guess for sigma_x
        sigma_y_initial = 10;             % Initial guess for sigma_y
        B_initial = min(data(:));         % Background offset

        initialGuess = [A_initial, x0_initial, y0_initial, sigma_x_initial, sigma_y_initial, B_initial];

        % Define the 2D Gaussian function
       gauss2D = @(params, xy) params(1) * exp( ...
        -((xy(:, 1) - params(2)).^2 / (2 * params(4)^2) + (xy(:, 2) - params(3)).^2 / (2 * params(5)^2)) ...
          ) + params(6);

       % Set optimization options
     options = optimset('Display', 'off');

     % Perform the fit using lsqcurvefit
    [paramFit, resnorm] = lsqcurvefit(gauss2D, initialGuess, xData, zData, [], [], options);

     % Extract the fitted parameters
    A_fit = paramFit(1);    % Amplitude (mean intensity)
    x0_fit = paramFit(2);   % x center
    y0_fit = paramFit(3);   % y center
    sigma_x_fit = paramFit(4); % Standard deviation in x
    sigma_y_fit = paramFit(5); % Standard deviation in y
    B_fit = paramFit(6);    % Background offset

    % Display the fitted parameters
   % disp('Fitted Parameters:');
   % fprintf('A (Mean Intensity) = %.3f, x0 = %.3f, y0 = %.3f, sigma_x = %.3f, sigma_y = %.3f, B = %.3f\n', ...
   %  A_fit, x0_fit, y0_fit, sigma_x_fit, sigma_y_fit, B_fit);

    % The mean intensity is the amplitude A
    meanIntensity = A_fit;
    disp(['Mean Intensity (from fit) = ', num2str(meanIntensity)]);
    
    % % Plot the original data and the fitted Gaussian
    % figure;
    % subplot(1, 2, 1);
    % imagesc(data);
    % title('Original Data');
    % axis image;
    % 
    % % Generate the fitted Gaussian for plotting
    % fittedData = reshape(gauss2D(paramFit, xData), size(data));
    % 
    % subplot(1, 2, 2);
    % imagesc(fittedData);
    % title('Fitted 2D Gaussian');
    % axis image;   
    int_vec(m,ii)=abs(A_fit);
    % pause(1);
    % close all;
    end
    
    waitbar(m/total_frames, f, sprintf('Time Series Analysis: %d %%', floor(m/total_frames*100)));
end
close(f)
%% STEP 5.1: Extract time series data from the 52 device locations (Normalized intensity (plot starts with 1 then decays))
f0 = waitbar(0,'Please wait...');
az=size(int_vec);
final_mat=zeros(total_frames,length(xpeaks));
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
    p=plot(total_time',final_mat(:,kk1),'-'); 
    p.Color=[cmap(kk1,1) cmap(kk1,2) cmap(kk1,3)];
    hold on;
    xlabel('Time (sec)')
    ylabel('Intensity (AU)')
    title(['Normalized Intensity Vs. Time plot for ',num2str(az(1,2)),' ','wells.'])
    fontsize(16,"points");
    axis tight;
    waitbar(kk1/jj1(1,2), f0, sprintf('Normalized Intensity vs. time plots: %d %%', floor(kk1/jj1(1,2)*100)));
end
hold off;
close(f0)
%% STEP 5.2: Normalized plot for a SIGNLE WELL
well_number=35; %% enter well number
xx=1:az(1,1);

plot(xx,final_mat(:,well_number),'--gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])

xlabel('Time (min.)')
ylabel('Normalized Intensity (au)')
title('Normalized Intensity vs. Time Plot')
%% STEP 6.1: Curve_fits1 (a*exp(-bx)) (OPTION 1)
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
    title(['Exponential Curve Fits for ',num2str(length(xpeaks)),' ','wells.'])
   legend off;
   fontsize(16,"points");
   waitbar(kk/jj(1,2), f2, sprintf('exp1 curve fits (a*exp(-bx)): %d %%', floor(kk/jj(1,2)*100)));
end
close(f2)
%% STEP 6.2: Curve_fits2 (a*exp(-bx)+c)  (OPTION 2)
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
    title(['Exponential Curve Fits for ',num2str(length(xpeaks)),' ','wells.'])
    legend off;
    fontsize(16,"points");
    waitbar(kk/jj(1,2), f6, sprintf('exp1 curve fits (a*exp(-bx)+c) : %d %%', floor(kk/jj(1,2)*100)));
end
close(f6)
coef_k1=coef_b;
%% STEP 6.3: Curve_fits3 (a exp(-bx)+c)  (OPTION 3)
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
   title(['Exponential Curve Fits for ',num2str(length(xpeaks)),' ','wells.'])
   fontsize(16,"points");
   waitbar(kk/jj(1,2), f7, sprintf('exp1 curve fits (a*exp(-bx)+c): %d %%', floor(kk/jj(1,2)*100)));
end
close(f7)
coef_k1=coef_b_star;

%% STEP 7: Check which wells shows change (reaction rates vs. well number)
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
%% STEP 8: Simple Histogram of rate constant (k)
nbins=10;
figure(8); 
h=histogram(abs(coef_k1),nbins); 
h(1).FaceColor = [1 0.5 0];
h(1).LineWidth =0.3;
h(1).EdgeColor = 'g';
%h(2).Color = [0 0 0];
xlabel('Rate constant K (sec.^{-1})')
ylabel('Frequency')
title(['Rate constant K for ',num2str(length(xpeaks)),' ','wells.'])
%set(gca,'YScale','log')

%% STEP 9: thresholded histogram
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

%% STEP 10: Sum of Gaussian of above Histogram for rate constant (k) calculation
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
title(['Rate constant K for ',num2str(length(xpeaks)),' ','wells.'])
%% OPTIONAL: Individual images of single well over time
img_no=13; % change the well number
n = 256; % Number of colors
fluorescent_green = [linspace(0, 0, n)' linspace(0, 1, n)' linspace(0, 0, n)'];
I_1=imread([pathname1,num2str(img_no),'.tiff']); 
I_2=double(I_1);
Ia=imcrop(I_2,[xpeaks(well_number)+15, ypeaks(well_number)+35 40 40]);
Ik1=Ia./15000;
imshow((Ik1)); colormap(fluorescent_green); 
% Construct the filename
filename = strcat(num2str(img_no), ".jpg");
% Write the image
imwrite(Ik1.*200,fluorescent_green, filename);




