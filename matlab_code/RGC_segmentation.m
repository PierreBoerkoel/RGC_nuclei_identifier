clear; close all;

infolder = '/Users/pierreboerkoel/Desktop/Glaucoma Project/images/Group 3 9M baseline/input_images/';
outfolder = '/Users/pierreboerkoel/Desktop/Glaucoma Project/images/Group 3 9M baseline/analysis_35/';

outliers = readtable(strcat(infolder, 'outlier_nuclei.csv'));
cy3threshold = 35;
curve_dilation_radius = 50;

imglist = dir(strcat(infolder, '**/*+2.tif'));
imgnames = extractfield(imglist, 'name');

cy3_table = table('Size', [length(imglist),4], 'VariableTypes', ...
   {'string', 'double', 'double', 'double'},'VariableNames', {'FileName', 'outliers_removed', 'cy3_threshold', 'cy3_percentage'});

for i = 1:length(imglist)

    filename = imgnames{i}; % string(filename(1:end-4));
    disp([num2str(i), '/', num2str(length(imglist)),' - ', filename]);
    I = imread(strcat(infolder, filename));

    Ib = I(:,:,3);

    % dapi threshold
    Ibthre = Ib>50;
    outThre = Ib>10; % to limit the connectedness of outlier nuclei in outlierMask

    % remove outlier nuclei
    outliersRemoved = '_no_outliers_removed_';
    if ismember(filename, outliers.FileName)
        outlierMask = ones(size(Ibthre));
        outs = outliers(ismember(outliers.FileName, filename),:);
        outsFound = nnz(outlierMask);
        for k = 1:size(outs, 1)
           outlierMask = outlierMask & ~bwselect(outThre, outs(k, :).x, outs(k, :).y);
        end
        Ibthre = Ibthre & outlierMask;
        Ib = double(Ib).*Ibthre;
        outliersRemoved = strcat('_', string(size(outs, 1)), '_outliers_removed_');
        cy3_table.outliers_removed(i) = size(outs, 1);
    end

    % close gaps to create clusters
    se = strel('disk',12);
    Icl = imdilate(Ibthre, se);
    Icl = imclose(Icl,se);

    % 3. identify clusters (CC)
    CC = bwconncomp(Icl,26);

    % 4. identify the DAPI band clusters
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [max1, ind1] = max(numPixels);
    numPixels(ind1) = 0;
    [max2, ind2] = max(numPixels);

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
    for k = 1:size(Ymask,2)
        topY = find(Ymasked(:,k),1);
        maskind = Y(:,k)+150 > topY;
        Ymask(maskind,k) = 0;
    end

    Ibcut = double(Ib).*Ymask;
    % 1. threshold
    Ibthre = Ibcut>30;
    Ithre_mf2 = medfilt2(Ibthre);
    Ithre_cleaned = imfill(logical(abs(Ithre_mf2-1)),[631 321]);

    % cy3 mask
    Im_cy3 = I(:,:,1);
    cy3_mask = Im_cy3 > cy3threshold;

    % Find cy3 staining using curve
    [row, col] = find(abs(Ithre_cleaned-1));
    curve = fit(col, row, 'smoothingspline','smoothingParam',0.000001);
    % lower_curve = fit(col, row+50, 'smoothingspline','smoothingParam',0.000001);
    % upper_curve = fit(col, row-50, 'smoothingspline','smoothingParam',0.000001);

    rgc_curve = zeros(size(cy3_mask));
    [numRows,numCols] = size(rgc_curve);

    for k = 1:numCols
        [xx,yy] = ndgrid((1:numRows)-ppval(coeffvalues(curve),k),(1:numCols)-k);
        rgc_mask = (xx.^2 + yy.^2)<curve_dilation_radius^2;
        rgc_curve(rgc_mask) = 1;
    end

    rgc_cy3_curve = rgc_curve & cy3_mask;
    rgc_percent_cy3_curve_method = nnz(rgc_cy3_curve)./nnz(rgc_curve);
    max_intensity_curve_sample = bsxfun(@times, rgc_cy3_curve, 255);
    max_intensity_curve_sample = double(cat(3, max_intensity_curve_sample, max_intensity_curve_sample, max_intensity_curve_sample));
    max_intensity_curve_sample(:,:,2) = 0;
    max_intensity_curve_sample(:,:,3) = 0;

    % Find cy3 staining using dilation
    % se = strel('disk', curve_dilation_radius, 8);
    % rgc_dilation = imdilate(abs(Ithre_cleaned-1), se);
    % rgc_dilation = imclose(rgc_dilation, 50);
    %
    % rgc_cy3_dilation = cy3_mask & rgc_dilation;
    % rgc_percent_cy3_dilation_method = nnz(rgc_cy3_dilation)./nnz(rgc_dilation);
    % max_intensity_dilation_sample = bsxfun(@times, rgc_cy3_dilation, 255);
    % max_intensity_dilation_sample = double(cat(3, max_intensity_dilation_sample, max_intensity_dilation_sample, max_intensity_dilation_sample));
    % max_intensity_dilation_sample(:,:,2) = 0;
    % max_intensity_dilation_sample(:,:,3) = 0;

    % update cy3 table
    cy3_table.FileName(i) = filename;
    cy3_table.cy3_threshold(i) = cy3threshold;
    cy3_table.cy3_percentage(i) = rgc_percent_cy3_curve_method;

    % Write figures
    f = figure('visible', 'off');
    subplot(2,1,1); imagesc(I(:,:,1:3)); hold on; plot(curve,'g'); title('Curve Method');
    % subplot(3,2,2); imagesc(abs(Ithre_cleaned-1)); title('Dilation Method');
    % subplot(3,1,2); imagesc(rgc_curve); axis square;
    % subplot(3,2,4); imagesc(rgc_dilation);
    subplot(2,1,2); imagesc(max_intensity_curve_sample);
    % subplot(3,2,6); imagesc(max_intensity_dilation_sample);

    saveas(f, strcat(outfolder, string(filename(1:end-4)), outliersRemoved, 'RGC.png'));
    % waitfor(figure);
    imwrite(max_intensity_curve_sample, strcat(outfolder, string(filename(1:end-4)), '_threshold-', num2str(cy3threshold), '_curve_image.tif'));
    % imwrite(max_intensity_dilation_sample, strcat(outfolder, string(filename(1:end-4)), '_threshold-', num2str(cy3threshold), '_dilation_image.tif'));
end

% Write cy3 data to table
writetable(cy3_table, strcat(outfolder, 'threshold-', num2str(cy3threshold), '_normalized_cy3_table.csv'));
