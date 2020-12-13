clear; close all;

filename = '9MWTF3 central-1-Image Export-01_c1+2.tif';
image_name = string(filename(1:end-4));
I = imread(strcat('/Users/pierreboerkoel/Desktop/Glaucoma Project/images/Group 1 Wildtype/', filename));
outliers = readtable('/Users/pierreboerkoel/git/Glaucoma_Image_Analysis/Outliers.csv');
outfolder = '/Users/pierreboerkoel/Desktop/Glaucoma Project/images/Group 1 Wildtype/';

cy3threshold = 20;
curve_dilation_radius = 50;

Ib = I(:,:,3);

% 1. otsu threshold
Ibthre = Ib>30;

%remove outlier nuclei
outliersRemoved = '_no_outliers_removed_';
if ismember(filename, outliers.FileName)
    outlierMask = ones(size(Ibthre));
    outs = outliers(ismember(outliers.FileName, filename),:);
    for i = 1:size(outs, 1)
       outlierMask = outlierMask & ~bwselect(Ibthre, outs(i, :).x, outs(i, :).y);
    end
    Ibthre = Ibthre & outlierMask;
    Ib = double(Ib).*outlierMask;
    outliersRemoved = strcat('_', string(i), '_outliers_removed_');
end

% 2. close gaps to create clusters   
se = strel('disk',12);
Icl = imdilate(Ibthre, se);
Icl = imclose(Ibthre,se);

% 3. identify clusters (CC) 
CC = bwconncomp(Icl,26);

% 4. identify the DAPI band clusters
numPixels = cellfun(@numel,CC.PixelIdxList);
[max1 ind1] = max(numPixels);
numPixels(ind1) = 0; 
[max2 ind2] = max(numPixels);

% 4.1. identify the top DAPI band cluster 
if max2/max1 < .2
    numLayers = 1;
    topInd = ind1;
else
    numLayers = 2;
    topInd = ind2;
end

% 5. find the top of top DAPI band cluster 
topchunk = zeros(size(Icl));
topchunk(CC.PixelIdxList{topInd}) = 1;

[X, Y] = meshgrid(1:size(topchunk,1));
Ymasked = Y.*topchunk;

Ymask = ones(size(Icl));
for i = 1:size(Ymask,2)
    topY = find(Ymasked(:,i),1);
    maskind = Y(:,i)+150 > topY;
    Ymask(maskind,i) = 0;
end

%figure; subplot(1,2,1); imagesc(Ib); subplot(1,2,2); imagesc(double(Ib).*Ymask);

%%%%%%%%%%%%%%%%%%%%%%%%%%

Ibcut = double(Ib).*Ymask;
% 1. threshold
Ibthre = Ibcut>30;
%figure; imagesc(Ibthre);

% 2. fit curve
%[row, col] = find(Ibthre);
%options = fitoptions('Method', 'SmoothingSpline', 'SmoothingParam', 0.000001);
%curve  = fit(col, row, 'SmoothingSpline', 'SmoothingParam', 0.000001);
%curve  = fit(col, row);
%figure; subplot(1,2,1); imagesc(I(:,:,1:3)); hold on; plot(curve,'-g'); subplot(1,2,2); imagesc(Ibthre);
%figure; subplot(1,2,1); imagesc(I(:,:,1:3)); subplot(1,2,2); imagesc( I(:,:,1:3).*uint8(repmat(Ibthre,1,1,3)));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% better identification of RGC
%se = strel('disk',3);
Ithre_mf2 = medfilt2(Ibthre);
Ithre_cleaned = imfill(logical(abs(Ithre_mf2-1)),[631 321]);

% cy3 mask
Im_cy3 = I(:,:,1);
cy3_mask = Im_cy3 > cy3threshold;

% Find cy3 staining using curve
[row, col] = find(abs(Ithre_cleaned-1));
curve = fit(col, row, 'smoothingspline','smoothingParam',0.000001);
lower_curve = fit(col, row+50, 'smoothingspline','smoothingParam',0.000001);
upper_curve = fit(col, row-50, 'smoothingspline','smoothingParam',0.000001);

rgc_curve = zeros(size(cy3_mask));
[numRows,numCols] = size(rgc_curve);
current_slope = 0;
for i = 1:numCols
%     if i < numCols
%         x_current = i;
%         y_current = ppval(coeffvalues(upper_curve),i);
%         x_next = i + 1;
%         y_next = -(ppval(coeffvalues(upper_curve),x_next));
%
%         slope = y_next - y_current;
%         if slope ~= 0
%             x_upper = round(x_next + sqrt((curve_distance^(2))./(1+(1/slope^(2)))));
%             x_lower = round(x_next - sqrt((curve_distance^(2))./(1+(1/slope^(2)))));
%             y_upper = round(y_next - ((1/slope)*(x_upper-x_next)));
%             y_lower = round(y_next - ((1/slope)*(x_lower - x_next)));
%         else
%             x_upper = x_current;
%             x_lower = x_current;
%             y_upper = round(y_current + curve_distance/2);
%             y_lower = round(y_current - curve_distance/2);
%         end
%     end

%     y_upper = round(ppval(coeffvalues(curve),i) - curve_dilation_radius);
%     y_lower = round(ppval(coeffvalues(curve),i) + curve_dilation_radius);
%
%     x_upper = i - curve_dilation_radius;
%     x_lower = i + curve_dilation_radius;

    [xx,yy] = ndgrid((1:numRows)-ppval(coeffvalues(curve),i),(1:numCols)-i);
    rgc_mask = (xx.^2 + yy.^2)<curve_dilation_radius^2;
    rgc_curve(rgc_mask) = 1;

%     if y_upper < 1
%         y_upper = 1;
%     elseif y_upper > numRows
%         continue;
%     end
%
%     if y_lower < 1
%         continue;
%     elseif y_lower > numRows
%         y_lower = numRows;
%     end
%
%     if x_lower < 1
%         x_lower = 1;
%     elseif x_lower > numCols
%         x_lower = numCols;
%     end
%
%     if x_upper < 1
%         x_upper = 1;
%     elseif x_upper > numCols
%         x_upper = numCols;
%     end
    %rgc_curve(y_upper:y_lower,x_upper:x_lower) = 1;
end

rgc_cy3_curve = rgc_curve & cy3_mask;
rgc_percent_cy3_curve_method = nnz(rgc_cy3_curve)./nnz(rgc_curve);
max_intensity_curve_sample = bsxfun(@times, rgc_cy3_curve, 255);
max_intensity_curve_sample = double(cat(3, max_intensity_curve_sample, max_intensity_curve_sample, max_intensity_curve_sample));
max_intensity_curve_sample(:,:,2) = 0;
max_intensity_curve_sample(:,:,3) = 0;

% Find cy3 staining using dilation
se = strel('disk', curve_dilation_radius, 8);
rgc_dilation = imdilate(abs(Ithre_cleaned-1), se);
rgc_dilation = imclose(rgc_dilation, 50);

rgc_cy3_dilation = cy3_mask & rgc_dilation;
rgc_percent_cy3_dilation_method = nnz(rgc_cy3_dilation)./nnz(rgc_dilation);
max_intensity_dilation_sample = bsxfun(@times, rgc_cy3_dilation, 255);
max_intensity_dilation_sample = double(cat(3, max_intensity_dilation_sample, max_intensity_dilation_sample, max_intensity_dilation_sample));
max_intensity_dilation_sample(:,:,2) = 0;
max_intensity_dilation_sample(:,:,3) = 0;

% Write cy3 data to table
cy3_table = table(image_name, rgc_percent_cy3_dilation_method, rgc_percent_cy3_curve_method);
writetable(cy3_table, strcat(outfolder, image_name, '_threshold-', num2str(cy3threshold), '_cy3_table.csv'));

% Write figures
figure = subplot(3,2,1); imagesc(I(:,:,1:3)); hold on; plot(upper_curve,'c'); plot(curve,'g'); plot(lower_curve,'y'); title('Curve Method');
subplot(3,2,2); imagesc(abs(Ithre_cleaned-1)); title('Dilation Method');
subplot(3,2,3); imagesc(rgc_curve);
subplot(3,2,4); imagesc(rgc_dilation);
subplot(3,2,5); imagesc(max_intensity_curve_sample);
subplot(3,2,6); imagesc(max_intensity_dilation_sample);

saveas(figure, strcat('/Users/pierreboerkoel/Desktop/Glaucoma Project/images/Group 1 Wildtype/', image_name, outliersRemoved, 'dilation_sample_RGC.png'));
imwrite(max_intensity_curve_sample, strcat(outfolder, image_name, '_threshold-', num2str(cy3threshold), '_curve_image.tif'));
imwrite(max_intensity_dilation_sample, strcat(outfolder, image_name, '_threshold-', num2str(cy3threshold), '_dilation_image.tif'));
