%%


% Load image
path_img = uigetdir; % Select path for files to be taken from
filePattern_img = fullfile(path_img, '*.tif'); % Only allow tif files from path to be selected
Din_img = dir(filePattern_img);
Ns_img = length(Din_img);
SelFiles_img = listdlg(...
            'PromptString', 'Choose Specific Files to Analyze', ...
            'SelectionMode', 'Multiple', ...
            'Name', 'File List', ...
            'InitialValue', 1:Ns_img, ...
            'ListString', {Din_img.name},...
            'Listsize',[300 400]);
       
        numImg = length(SelFiles_img); % Total number of images selected
        
                
% Read Image Data and Info
Images_projections = cell(numImg,1);
Images_projections{1} = imread(fullfile(path_img,Din_img(SelFiles_img(1)).name));

% Define strain conditions, change lines 25+26 to match desired strain
applied_strain = .1; % --> strain, this would be 10% strain
poissons_ratio = 1; 
fixed = Images_projections{1};
tform = imregcorr(fixed,Images_projections{1});

% apply digital strain transformation matrix
tform.T = [1+applied_strain 0 0; 0 1-poissons_ratio*applied_strain 0; 0 0 1];
warped = imwarp(fixed,tform);

%show warped image and fixed image
imshowpair(fixed,warped)
title('Warped, uncropped') 

% crop to equivalent size of fixed
padsize = size(warped);
for k = 1:(1024-padsize(1))
    warped(padsize(1)+k,:) = 0;
end
warped = imcrop(warped,[0 0 1024 1024]);

% show cropped, warped image and fixed image
figure
imshowpair(fixed,warped)

% Name as desired
imwrite(warped,'10_dsi.tif')
title('Warped,Cropped')

