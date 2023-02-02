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

%{
% Define strain conditions, change lines 25+26 to match desired strain
applied_strain = .1; % --> strain, this would be 10% strain
poissons_ratio = 1; 
%}

prompt = {'Enter Applied Strain Increment','Enter Max Applied Strain:', 'Enter Poisson''s'' Ratio'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0.02','0.1','1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
si = str2double(answer(1));
mas = str2double(answer(2));
pr = str2double(answer(3));

for k = 0:(mas/si)
    fixed = Images_projections{1};
    tform = affinetform2d;
    if k == 0
        titlek=sprintf('DigitallyTransformed%d%%Strain.tif',k*si*100);
        imwrite(fixed,titlek)
    else
        % apply digital strain transformation matrix
        tform.A = [1+k*si 0 0; 0 1-pr*k*si 0; 0 0 1];
        warped = imwarp(fixed,tform);
    
        % crop to equivalent size of fixed
        padsize = size(warped);
        for m = 1:(1024-padsize(1))
            warped(padsize(1)+m,:) = 0;
        end
        warped = imcrop(warped,[0 0 1024 1024]);
    
        % show cropped, warped image and fixed image
        figure
        imshowpair(fixed,warped)
    
        %Name as desired
        titlek=sprintf('DigitallyTransformed%d%%Strain.tif',k*si*100);
        imwrite(warped,titlek)
    end
end
