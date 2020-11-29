clear; close all;
%I = imread('Y-8 peripheral-1-Image Export-52_c1+2.tif');
% I = imread('Y-1 peripheral-1-Image Export-40_c1+2.tif');
 I = imread('T-2 central-1-Image Export-22_c1+2.tif');
%I = imread('9MWTF3 peripheral-1-Image Export-03_c1+2.tif');

Ib = I(:,:,3); 

% binarize 

% 1. otsu threshold
%Ibthre = imbinarize(Ib);
Ibthre = Ib>30;
%figure; imagesc(Ibthre); 

% 2. close gaps to create clusters   
se = strel('disk',12);
Icl = imdilate(Ibthre, se); 
Icl = imclose(Ibthre,se);
%figure; imagesc(Icl); 

% 3. idnentify clusters (CC) 
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
[row, col] = find(Ibthre); 
curve  = fit(col, row, 'smoothingspline','smoothingParam',0.000001);
figure; subplot(1,2,1); imagesc(I(:,:,1:3)); hold on; plot(curve,'-g'); subplot(1,2,2); imagesc(Ibthre); 
%figure; subplot(1,2,1); imagesc(I(:,:,1:3)); subplot(1,2,2); imagesc( I(:,:,1:3).*uint8(repmat(Ibthre,1,1,3)));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% better identification of RGC 
se = strel('disk',3);
Ithre_mf2 = medfilt2(Ibthre); 
Ithre_cleaned = imfill(logical(abs(Ithre_mf2-1)),[631 321]); 

[row, col] = find(abs(Ithre_cleaned-1)); 
curve  = fit(col, row, 'smoothingspline','smoothingParam',0.000001);
figure; subplot(1,2,1); imagesc(I(:,:,1:3)); hold on; plot(curve,'-g'); subplot(1,2,2); imagesc(Ithre_mf2); 


