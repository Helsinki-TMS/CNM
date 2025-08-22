function runClusterAnalysis(subject,sessionPrefix)

addpath('/m/nbe/scratch/istim/Tuukka/programs/NIfTI_20140122');

if isempty(sessionPrefix)
    session      = [sessionPrefix 'FunImgARCFB'];
    result       = [sessionPrefix 'Results'];
else
    session      = [sessionPrefix '_FunImgARCFB'];
    result       = [sessionPrefix '_Results'];
end

datapath         = ['/m/nbe/scratch/istim/3Rpilot/' subject];
outputPath       = [datapath,'/cluster/',result,'/FC_FunImgARCFB'];

fMRI_prefix      = [subject '_4DVolume.nii'];
fMRI_path        = [datapath '/connectivity/' session '/' subject '/' fMRI_prefix];

seed_prefix      = ['FCROI_9_' subject '.nii'];
seed_path        = [datapath '/connectivity/Masks/WarpedMasks/' seed_prefix];

brainMask_prefix = [subject '_BrainMask_05_91x109x91.nii'];
brainMask_path   = [datapath '/connectivity/Masks/WarpedMasks/' brainMask_prefix]; 

numberOfClusters = 8;


%% Create output directory

mkdir(outputPath);

%% Read fMRI image
nii_fMRI    = load_untouch_nii(fMRI_path);

%% Smooth image
smoothImg   = zeros(size(nii_fMRI.img));
for i=1:size(nii_fMRI.img,4)
    tmpImg = double(nii_fMRI.img(:,:,:,i));
    smoothImg(:,:,:,i) = imgaussfilt3(tmpImg,0.6667); %approx. 4 mm smoothing
end

nii_fMRI.img               = single(smoothImg);
nii_fMRI.hdr.dime.datatype = 16;
sfMRI_prefix               = ['s', fMRI_prefix];
save_untouch_nii(nii_fMRI,[outputPath,'/',sfMRI_prefix,'.gz']);

%% Read seed image
nii_seed    = load_untouch_nii(seed_path);

%% Find indices of seed points in the image
indices     = find(nii_seed.img(:));
[x,y,z]     = ind2sub(size(nii_seed.img),indices);

%% Get connectivity matrix
N           = size(x,1);

conn        = zeros(N,N);

for i = 1:N
    for j = 1:N
        
        fMRI_p1   = smoothImg(x(i),y(i),z(i),:);
        fMRI_p1   = double(squeeze(fMRI_p1));
        
        fMRI_p2   = smoothImg(x(j),y(j),z(j),:);
        fMRI_p2   = double(squeeze(fMRI_p2));
         
        conn(i,j) = corr(fMRI_p1,fMRI_p2, 'type', 'Spearman');
        
    end
end

%% Plot connectivity matrix
threshold   = 0.5;
conn_thresh = conn.*(conn > threshold);
imagesc(conn_thresh);
axis square;

%% Check time series visually
% i = 1;
% j = 10;
% 
% fMRI_p1   = smoothImg(x(i),y(i),z(i),:);
% fMRI_p1   = squeeze(fMRI_p1);
% fMRI_p2   = smoothImg(x(j),y(j),z(j),:);
% fMRI_p2   = squeeze(fMRI_p2);
% 
% corr(fMRI_p1,fMRI_p2, 'type', 'Spearman');
% 
% figure;
% % plot(fMRI_p1,'k'); hold on;
% % plot(fMRI_p2,'r');
% scatter(fMRI_p1,fMRI_p2,'or');

%% Clustering

Z = linkage(conn_thresh);
T = cluster(Z,'maxclust',numberOfClusters);

numberOfVoxelsInClusters = zeros(numberOfClusters,1);   
for i=1:numberOfClusters
    numberOfVoxelsInClusters(i) = nnz(T==i);
end

%% Write clusters in image

cluster_img             = int16(zeros(size(nii_seed.img)));

for i=1:numberOfClusters
    cluster_img(indices(T==i))  = int16(i);
end

nii_seed.img = cluster_img;
save_untouch_nii(nii_seed,[outputPath,'/','clusterLabels_', seed_prefix,'.gz']);

%% Time series for each cluster, saved as "clusterTimeSeries"

voxelClusterTimeSeries       = cell(numberOfClusters,1);

for i=1:length(indices)
    voxelTimeSeries         = squeeze(smoothImg(x(i),y(i),z(i),:));
    voxelClusterTimeSeries{T(i)} = [voxelClusterTimeSeries{T(i)};voxelTimeSeries'];
end

medianClusterTimeSeries = zeros(numberOfClusters,size(smoothImg,4));
clusterTimeSeries       = zeros(numberOfClusters,size(smoothImg,4));
for i=1:numberOfClusters
    medianClusterTimeSeries(i,:)  = median(voxelClusterTimeSeries{i},1)';
    maxCorr = 0;
    for j=1:size(voxelClusterTimeSeries{i},1)
        curCorr = corr(voxelClusterTimeSeries{i}(j,:)',medianClusterTimeSeries(i,:)', 'type', 'Spearman');
        if (curCorr>maxCorr)
            clusterTimeSeries(i,:) = voxelClusterTimeSeries{i}(j,:);
        end
    end
end

%% Connectivity with the rest of the brain

nii_brainMask    = load_untouch_nii(brainMask_path);

indices          = find(nii_brainMask.img(:));
[x,y,z]          = ind2sub(size(nii_brainMask.img),indices);

cluster2brainConn= zeros(numberOfClusters,length(indices));

parfor i=1:length(indices)
    for c=1:numberOfClusters

        fMRI_b   = smoothImg(x(i),y(i),z(i),:);
        fMRI_b   = double(squeeze(fMRI_b));
        
        cluster2brainConn(c,i) = corr(fMRI_b,clusterTimeSeries(c,:)', 'type', 'Spearman');
    end
end

%% Scale cluster connectivities with cluster volume and save as nifti images
for c=1:numberOfClusters
    img                        = single(zeros(size(smoothImg)));
    img(indices)               = numberOfVoxelsInClusters(c)*cluster2brainConn(c,:);
    
    nii_seed.img               = single(img);
    nii_seed.hdr.dime.datatype = 16;
    
    cfMRI_prefix               = ['ROI', num2str(c), 'FCMap_', subject '.nii'];
    save_untouch_nii(nii_seed,[outputPath,'/',cfMRI_prefix]);
    
    % Save also z-scored image
    nii_seed.img = (nii_seed.img(indices)-mean(nii_seed.img(indices)))/std(nii_seed.img(indices));
    
    cfMRI_prefix               = ['zROI', num2str(c), 'FCMap_', subject '.nii'];
    save_untouch_nii(nii_seed,[outputPath,'/',cfMRI_prefix]);
    
end
