 
%% Process TIF Images from Microscope
clearvars -except ini_path Nlines
close all
scrsz = get(0,'screensize');

% New algortihm identifies sample surface and takes pixel values from that
% plane. 
startImg=1;

% Thresholds
edgeT= 1.5; %1.4;      %1.4 fordtissue, 1.2 for gels
areaT=200;      %200

distT=80; %100;      %150 25x, 100 for 10xareaT

% Number of pixels over which to smooth lines (must be odd)
avg_L=35;

% Run left sides of Notched samples
Left=false;

%% Read Data
if ~exist('ini_path','var') || max(ini_path==0)
    ini_path='C:\Users\bep15\Box\b-szczesny-lab Shared\Student Folders\Ben\Data\Grip_assesment_multi-PBL\Grip_strain_test_PBL\0%';
end
[file,ini_path,type] = uigetfile({'*.lsm;*.tif;*.nd2;*.mat;*.fig','LSM(*.lsm), TIFF(*.tif),ND2(*.nd2), MAT-file(*.mat), or Figure(*.fig)'},'Choose LSM File, Image File from TIFF Set, or .mat File Containing Images or Figure',ini_path);
[path, name, ext] = fileparts(fullfile(ini_path,file));
cd(ini_path);

switch ext
    case '.nd2' 
        
   file_name_ = name(1:end);
   file_name = [file_name_, '.nd2'];
   
 %Import Metadata (Gives Stage Positions, resolultion (um/px), and z-stack
 %slice thickness. 
 %Uses imreadBFmeta function. Download link found in line 46
 meta = imreadBFmeta(file);
 
 z_1 = getfield(meta,'zsize'); %Outpus n-stack size 
 z_2 = getfield(meta,'nframes');  % Outputs z-stack size
 
 numI = max(z_1,z_2)    ;  
 
 parm = getfield(meta,'parameterNames');
 parm_ = char(parm);
 parm__ = cellstr(parm_)
 y_name_pos= find(strcmp('Series 0 Y position 1',parm__),1);
 x_name_pos= find(strcmp('Series 0 X position 1',parm__),1);
 
 d_cal_name = find(strcmp('dCalibration',parm__),1);
 d_z_name = find(strcmp('dZStep',parm__),1);
 
 
%Index value for given meta data location is fixed & input below.
values = getfield(meta,'parameterValues');
x_stg_pos = values(x_name_pos,1);  %Stage x position 
y_stg_pos = values(y_name_pos,1);  %Stage y position 
resolution = values(d_cal_name,1) ;    %um/px value
dz = values(d_z_name,1);             % z-stack thickness 
 
Res = [resolution, dz];

  % Read Images
 colorI = cell(numI,1);
 
 %import data using imreadBF function. Download found below
 %(https://www.mathworks.com/matlabcentral/fileexchange/32920-imread-for-multiple-life-science-image-file-formats)
 
 colorI = imreadBF(file_name,[1:z_1],[1:z_2],1); 

  % Convert to gray scale
  dimensions = size(colorI);
  pix_width = dimensions(1);
  pix_length = dimensions(2);
 
  I = uint16(colorI);
  I = rot90(I,-1);
  I = I*2^4 ;
    
  Isize = pix_width  
  
  Ns = numI
  
  prefix = file_name
  parent=fullfile(path,'Exported',[file_name_ 'data']);
  if ~exist(parent,'dir')
      mkdir(parent)
  end
        
    case '.lsm'
%         ind=find(name==' ',1,'first');
%         file_name=name(ind+1:end);        % When using "z" to denote z-stack in file name
        file_name=name(1:end);
        
        ind=regexp(file_name,' \d');
        max_ind=ind;%min(ind)-1;
        prefix=file_name(1:max_ind);
        parent=fullfile(path,'Exported',[prefix 'data']);
        
        if ~exist(parent,'dir')
            mkdir(parent);
        end
        
        d=dialog('Visible','off','WindowStyle', 'modal','Position',[50 50 200 100]);
        movegui(d,'center')
        uicontrol(d,'Style','text','String','Please Wait...Loading Data','Position',[1,35,200,20],'FontSize',10)
        set(d,'Visible','on')
        pause(.5)
        
        stack = tiffread(fullfile(path,[name ext]));
        numI=length(stack);
        Isize=size(stack(1).data,1);
        I=zeros([Isize,Isize,numI],'uint16');
        for i=startImg:numI
            I(:,:,i)=stack(i).data;
%             I(:,:,i)=fliplr(stack(i).data);
        end
        
%         Res=nan(1,2);    %Xres, Yres
%         % Read scaling
%         Res(1)=10^6*stack(1).lsm.VoxelSizeX;        %um/pix (x-dir)
%         Res(2)=10^6*stack(1).lsm.VoxelSizeZ;        %um/pix (z step size)
        
        close(d)
        
    
    case '.tif'
        [parent, file_name] = fileparts(path);
        ind=regexp(file_name,'\d');
        max_ind=min(ind)-1;
        prefix=file_name(1:max_ind);
        parent=fullfile(parent,[prefix 'data']);
        if ~exist(parent,'dir')
            mkdir(parent);
        end
        
        Din = dir(fullfile(path,'*.tif'));
        Ns = length(Din);
        
       % Res=[1.23 2.175] ; UPDATE FOR EACH INDIVDUAL TIFF FILE. DOESN NOT
       % AUTOMATICALLY UPDATE 
        
        % %Prompt user to choose which files to analyze
        % SelFiles = listdlg(...
        %     'PromptString', 'Choose Specific Files to Analyze', ...
        %     'SelectionMode', 'Multiple', ...
        %     'Name', 'File List', ...
        %     'InitialValue', 1:Ns, ...
        %     'ListString', {Din.name},...
        %     'Listsize',[300 400]);
        %
        % numI = length(SelFiles);
        
        % Select all files
        numI = Ns;
        
        % Read Images
        colorI = cell(numI,1);
        w = waitbar(0,'Reading Image Files');
        for i=1:numI
            colorI{i} = imread(fullfile(path,Din(i).name));
            waitbar(i/numI)
        end
        close(w)
        
        %Convert to grayscale
        Isize=size(colorI{1},1);
        I = zeros([Isize,Isize,numI],'uint16');
        for i=1:numI
            I(:,:,i) = colorI{i};
        end
        I=I*2^4;
        
%         % Save Images Array
%         if ~exist([fullfile(parent,file_name) ' I.mat'],'file')
%             save([fullfile(parent,file_name) ' I'], 'I')
%         end

    case '.mat'
        parent=path;
        file_in=name;
        file_name=name(1:find(file_in==' ',1,'last')-1);
        
        d=dialog('Visible','off','WindowStyle', 'modal','Position',[50 50 200 100]);
        movegui(d,'center')
        uicontrol(d,'Style','text','String','Please Wait...Loading Data','Position',[1,35,200,20],'FontSize',10)
        set(d,'Visible','on')
        pause(.5)
        load(fullfile(parent,[file_in ext]))
        Isize=size(I,1);
        close(d)
        
    case '.fig'
        parent=path;
        file_in=name;
        file_name=name(1:find(file_in==' ',1,'last')-1);
        
        open(fullfile(parent,[file_in ext]))
        line_data=get(gcf,'UserData');
        template=line_data.template;
        Isize=size(template,1);
        
        sides=line_data.sides;
        xL=min(sides{1}(:,1));
        xR=max(sides{2}(:,1));
        
        X=line_data.X;
        Y=line_data.Y;
    otherwise
        return
end

%% Profile Template
if ~exist('line_data','var')
% Find sides of sample
temp_value = round(Ns/2);
imshow(I(:,:,temp_value),'InitialMagnification',65)      %Changed 1/2/18 to use middle image to get better edge boundaries % imshow(I(:,:,(Ns/2)),'InitialMagnification',65)     Use final image
hold on
if strcmp(prefix,'cut ')    % Use zen for cut samples
    answer = inputdlg({'x-min:'},'Left Bound',1,{''});
    temp=str2double(answer{1});
    msgbox('Define Right Side. Right-click to End Each Line.','modal');
else
    msgbox('Define Left and then Right Sides. Right-click to End Each Line.','modal');
end

sides=cell(2,1);
for i=1:2
    if strcmp(prefix,'cut ') && i==1
        sides{i}=[temp 1;temp Isize];
    else
        [x,y]=getline;
        if atand(abs(y(2)-y(1))/abs(x(2)-x(1)))<45
            if i==1
                x=ones(2,1);
            elseif i==2
                x=Isize*ones(2,1);
            end
            y=[1;Isize];
        end
        sides{i}=[x,y];
    end
    plot(sides{i}(:,1),sides{i}(:,2),'g')
end
pause(.5)
close(gcf)
pause(.5)

xL=max(sides{1}(:,1));
xR=min(sides{2}(:,1));

end
%% Create composite image

if ~exist('line_data','var')
    
button='No';
while(strcmp(button,'No'))
template=im2uint16(zeros(Isize));
X=cell(Isize,1);    Y=(1:Isize)';

Wh = waitbar(0,'Processing Images...');
for x_pos=1:Isize
waitbar(x_pos/Isize)
    
% Longitudinal profile
Iprofile=im2uint16(zeros(Isize,numI));
for i=1:numI
    Iprofile(:,numI-(i-1))=I(:,x_pos,i);
end
Iratio=[Res(1) Res(2) 1];    %y-dir is first index(rows)

% if x_pos>=Bmin && x_pos<=Bmax
%     template(:,x_pos)=Iprofile(:,end);
%     continue
% end

% Determine Edges
[~,thresh]=edge(Iprofile(1:512,:),'sobel','vertical'); 
[iniEDGEStop,thresh]=edge(Iprofile(1:512,:),'sobel',edgeT*thresh,'vertical');           %1.6

[~,thresh]=edge(Iprofile(513:end,:),'sobel','vertical'); 
[iniEDGESbot,thresh]=edge(Iprofile(513:end,:),'sobel',edgeT*thresh,'vertical');

iniEDGES=cat(1,iniEDGEStop,iniEDGESbot);

% Remove noise
SE=strel('rectangle',[15,5]);  % 15,10 (drop 15 down if centroids are too big)
tempI=imdilate(iniEDGES,SE);
SE=strel('rectangle',[30,1]);  
EDGES=imerode(tempI,SE);

[rowI,colI]=find(iniEDGES);
iniCoords=[colI,rowI];     % y-dir is row
[rowE,colE]=find(EDGES);
Coords=[colE,rowE];

% Find centroids of edges that are above a certain size and are vertical
stats = regionprops(EDGES,'Centroid','Area','BoundingBox','Eccentricity','Orientation');
centroids=nan(length(stats),2);
for i=1:length(stats)
    if stats(i).Area>areaT        % 200
        if stats(i).Eccentricity<.92
            centroids(i,:)=stats(i).Centroid;
        else
            if abs(stats(i).Orientation)>85
                centroids(i,:)=stats(i).Centroid;
            end
        end
    end
end
centroids=centroids(~isnan(centroids(:,1)),:);

if isempty(centroids)
    if x_pos<xL || x_pos>xR
        template(:,x_pos)=Iprofile(:,1);      % last image
        X{x_pos}=ones(Isize,1);
        continue
    end
end

% Only keep centroids closest to surface (first clean up)
temp=nan(size(centroids));
for i=1:size(centroids,1)
    diff=abs(centroids(i,2) - centroids(:,2));      %y-dir
    inds1=find(diff<distT);           %150
    [~,inds2]=max(centroids(inds1,1));
    temp(i,:)=centroids(inds1(inds2),:);
end
centroids=unique(temp,'rows');

% Break up large regions
nbreak=1;
if size(centroids,1)==1
    if x_pos>=xL && x_pos<=xR
        for i=1:nbreak
            EDGES(i*Isize/(nbreak+1),:)=zeros(1,numI);
        end
        [rowE,colE]=find(EDGES);
        Coords=[colE,rowE];
        
        stats = regionprops(EDGES,'Centroid','Area','BoundingBox','Eccentricity','Orientation');
        centroids=nan(length(stats),2);
        for i=1:length(stats)
            if stats(i).Area>areaT        % 200
                if stats(i).Eccentricity<.92
                    centroids(i,:)=stats(i).Centroid;
                else
                    if abs(stats(i).Orientation)>85
                        centroids(i,:)=stats(i).Centroid;
                    end
                end
            end
        end
        centroids=centroids(~isnan(centroids(:,1)),:);
    end
end

% Only keep centroids closest to surface (second clean up)
temp=nan(size(centroids));
for i=1:size(centroids,1)
    diff=abs(centroids(i,2) - centroids(:,2));      %y-dir
    inds1=find(diff<distT);           %150
    [~,inds2]=max(centroids(inds1,1));
    temp(i,:)=centroids(inds1(inds2),:);
end
centroids=unique(temp,'rows');

if size(centroids,1)<2      % Note empty condition for image outside fascicle bounds is already covered above
    if x_pos<xL || x_pos>xR
        temp=round(centroids(1));
        template(:,x_pos)=Iprofile(:,temp);      % fixed depth
        X{x_pos}=temp*ones(Isize,1);
        continue
    else

    % If there still is only one centroid within sample edges
    % Use previous line if there is only one region    
    if length(stats)==1
        pF=pFold;
        pF(2)=centroids(1)-pF(1)*centroids(2);
        
        p=[1/pF(1) -pF(2)/pF(1)];
        xtop=round((1-p(2))/p(1));
        xbot=round((Isize-p(2))/p(1));
        
        centroids=cat(1,centroids,[xtop 1],[xbot Isize]);
        
        imshow(Iprofile)
        set(gca,'DataAspectRatio',Iratio)    %y-dir is first index(rows)
        hold on
        plot(iniCoords(:,1),iniCoords(:,2),'r.')
        plot(Coords(:,1),Coords(:,2),'w.')
        plot(centroids(:,1),centroids(:,2),'b.')
        set(gcf,'Position',[365 85 1188 881])
        
        plot([xtop xbot],[1 1024],'g-')
        hold off
        keyboard
        close(gcf)
        
    else
        % Select additional centroids
        imshow(Iprofile)
        set(gca,'DataAspectRatio',Iratio)    %y-dir is first index(rows)
        hold on
        plot(iniCoords(:,1),iniCoords(:,2),'r.')
        plot(Coords(:,1),Coords(:,2),'w.')
        plot(centroids(:,1),centroids(:,2),'b.')
        set(gcf,'Position',[365 85 1188 881])
        
%         plot([xtop xbot],[1 1024],'g-')
        hold off
        
        msgbox('Select Additional Centroids. Press Enter When Done.','modal');
        [tempx, tempy]=getpts;
        close(gcf)
        
        for i=1:length(tempx)
            
            for j=1:length(stats)
                temp=stats(j).BoundingBox;  % [xmin ymin width height]
                if tempx(i)>=temp(1) && tempx(i)<=temp(1)+temp(3) && tempy(i)>=temp(2) && tempy(i)<=temp(2)+temp(4)     % Inside bounding box
                    centroids=cat(1,centroids,stats(j).Centroid);
                end
            end
            
        end
        
        if size(centroids,1)==1
            temp=round(centroids(1));
            template(:,x_pos)=Iprofile(:,temp);      % fixed depth
            X{x_pos}=temp*ones(Isize,1);
            continue
        elseif size(centroids,1)==0
            template(:,x_pos)=Iprofile(:,1);      % fixed depth
            X{x_pos}=ones(Isize,1);
            continue
        end
        
    end
        
    end
end

% Fit line through centroids to find surface
% Fit needs to be flipped since fitting a vertical line is inaccurate
pF = polyfit(centroids(:,2),centroids(:,1),1);
X{x_pos} = polyval(pF,1:Isize);
pFold=pF;

p=[1/pF(1) -pF(2)/pF(1)];
xtop=round((1-p(2))/p(1));
xbot=round((Isize-p(2))/p(1));
numimgs=abs(xbot-xtop)+1;

% Noisy centroids that create extreme angles
if xtop<1 || xbot<1 || xtop>numI || xbot>numI
    temp=round(mean(centroids(:,1)));
    template(:,x_pos)=Iprofile(:,temp);      % avg depth
    X{x_pos}=temp*ones(Isize,1);
    continue
end

% Find pixel locations
% [X{x_pos}, ~]=bresenham(xtop,Y(1),xbot,Y(end));

% Capture pixel intensities
% True profiles
% if numimgs>2    % 2
% line=im2uint16(zeros(Isize,1));
% for i=1:Isize
%     line(i)=Iprofile(Y(i),X{x_pos}(i));
% end
% 
% else
% % Averaged profile to remove aliasing at center line
% line_sum=zeros(Isize,1);
% for i=min([xtop xbot]):max([xtop xbot])
%     line_sum=line_sum+im2double(Iprofile(:,i));
% end
% line=im2uint16(line_sum/numimgs);
% end

% Gradient Averaging
tempI=im2double(Iprofile);
line=zeros(Isize,1);
for i=1:Isize
    tempX=X{x_pos}(i);
    botI=floor(tempX);
    topI=ceil(tempX);
    if botI==0
        botI=1;
    end
    if topI>numI
        topI=numI;
    end
    line(i)=(ceil(tempX)-tempX)*(tempI(i,botI)) + (tempX-floor(tempX))*(tempI(i,topI));
end
line=im2uint16(line);

template(:,x_pos)=line;
end
close(Wh)

figure
 imshow(template,'InitialMagnification',65)
hold on
plot(sides{1}(:,1),sides{1}(:,2),'g')
plot(sides{2}(:,1),sides{2}(:,2),'g')

button = questdlg('Are sample boundaries acceptable?');

if strcmp(button,'No')
    old_sides=sides;
    sides=cell(2,1);
    for i=1:2
        if strcmp(prefix,'cut ') && i==1
            sides{i}=old_sides{i};
        else
            [x,y]=getline;
            if atand(abs(y(2)-y(1))/abs(x(2)-x(1)))<45
                if i==1
                    x=ones(2,1);
                elseif i==2
                    x=Isize*ones(2,1);
                end
                y=[1;Isize];
            end
            sides{i}=[x,y];
        end
        plot(sides{i}(:,1),sides{i}(:,2),'r')
    end
    pause(.5)
    
    xL=max(sides{1}(:,1));
    xR=min(sides{2}(:,1));
end
    close(gcf)
end

end
%}
%% Sketch Initial Guess Lines and Find Pixel Values of Lines
if ~exist('line_data','var') || Left
% Half-interval for search of minimum pixel intensity (next section)
range=15; %pixels (no more than 20)  

line_fig=figure;
% If guess lines already exist, ask if they should be redrawn
if exist('x_pts','var')
    Nlines=length(x_pts);
    figure(gcf)
     imshow(template,'InitialMagnification',65)
    hold on
    for i=1:Nlines
        plot(x_pts{i},y_pts{i},'r')
        plot(x_pts{i},y_pts{i}+range,'g',x_pts{i},y_pts{i}-range,'g')
    end
    hold off
    
    redo_lines = questdlg('Redo Guess Lines?');
else
    redo_lines='Yes';
end

switch redo_lines
    case 'Yes'
        figure(gcf)
         imshow(template,'InitialMagnification',65)
         zoom(1.0)
        hold on
        
        % Sketch a polyline for each bleach line
        if ~exist('Nlines','var')
            Nlines=4;
        end
        answer = inputdlg({'Number of Lines:'},'',1,{num2str(Nlines)});
        Nlines=str2double(answer{1});
        msgbox('Draw Multi-Point Line for each Bleach Line. Right-click to End Each Line.','modal');
        clear x_pts y_pts
        for i=1:Nlines
            if i>1
                plot(x_pts{i-1},y_pts{i-1},'r')
                plot(x_pts{i-1},y_pts{i-1}+range,'y',x_pts{i-1},y_pts{i-1}-range,'y')
            end
            [x_pts{i}, y_pts{i}] = getline(gca);
        end
        plot(x_pts{end},y_pts{end},'r')
        plot(x_pts{end},y_pts{end}+range,'y',x_pts{end},y_pts{end}-range,'y')
        hold off
        
        % Add guess lines to data file
%         d=dialog('Visible','off','WindowStyle', 'modal','Position',[50 50 200 100]);
%         movegui(d,'center')
%         uicontrol(d,'Style','text','String','Please Wait...Saving Data','Position',[1,35,200,20],'FontSize',10)
%         set(d,'Visible','on')
        pause(.5)
%         save([fullfile(parent,file_name) ' I'], 'I', 'x_pts', 'y_pts')
%         close(d)
        
    case 'Cancel'
        return

end

% Create a dummy imline element for creation of masks
 imshow(template,'InitialMagnification',65)
 zoom(1.0)
msgbox('Draw Dummy Line.','modal');
h = imline(gca);

% Create masks for each bleach line by adding masks of each line segment
% Find x,y positions of line masks
mask_lines=false(size(template,1),size(template,2),Nlines);
row=cell(1,Nlines);  col=row;
for i=1:Nlines
    num_pts=length(x_pts{i});
    for j=1:num_pts-1
        setPosition(h,x_pts{i}(j:j+1),y_pts{i}(j:j+1))
        temp_mask = createMask(h);
        % Mask for each line
        mask_lines(:,:,i) = logical(imadd(mask_lines(:,:,i),temp_mask));
    end
    % x,y pixel locations of each initial guess line
    [row{i},col{i}] = find(mask_lines(:,:,i));
end
set(line_fig,'Visible','off')

end
%% Find Line Locations based on minimum local pixel intensity
if ~exist('line_data','var') || Left
COL=cell(size(col));
ROW=cell(size(row));
% Find images to use for identifying pixel locations of bleach lines at
% every x-pos of lines
for i=1:Nlines   % for each line
    num_pts=length(col{i});
    x_pos = col{i};
    guess = row{i};
    temp_col=x_pos;      temp_row=nan(size(temp_col));
    for j=1:num_pts     % for each x-pos of each line
        if x_pos(j) < min(sides{1}(:,1)) || x_pos(j) > max(sides{2}(:,1))
            temp_col(j)=nan;
        end
               
        if ~isnan(temp_col(j))
            line=template(:,x_pos(j));
            
            % Find min intensity point about guess pixel
            [min_int, min_local] = min(line(guess(j)-range:guess(j)+range));
            contrast = double((max(line(guess(j)-range:guess(j)+range))-min_int))/mean(line(guess(j)-range:guess(j)+range));

            if contrast >= 0.1 && min([line(guess(j)-range), line(guess(j)+range)])>min_int % good contrast and true minimum
                temp_row(j)=guess(j)-range + min_local-1;
            end
        
            % Debugging
%             if i==4
%             if contrast < 0.1
%             if j==140
%                 figure(gcf)
%                 set(gcf,'Position',[100 330 1750 650],'Name',['Line ' num2str(i) ', x-pos ' num2str(x_pos(j)) ', contrast ' num2str(contrast)])
%                 subplot(1,2,1)
%                  imshow(template,'InitialMagnification',65)
%                 hold on
%                 plot([temp_col(j) temp_col(j)],[1 Isize])
%                 plot(temp_col(j),guess(j)-range,'y.',temp_col(j),guess(j)+range,'y.')
%                 plot(temp_col(j),guess(j),'r.')
%                 plot(temp_col(j),temp_row(j),'g.')
%                 hold off
%                 
%                 subplot(1,2,2)
%                 plot(1:Isize,line)
%                 hold on
%                 plot(guess(j)-range,line(guess(j)-range),'k.',guess(j)+range,line(guess(j)+range),'k.','MarkerSize',20)
%                 plot(guess(j),line(guess(j)),'r.')
%                 if ~isnan(temp_row(j))
%                     plot(temp_row(j),line(temp_row(j)),'g.')
%                 end
%                 hold off
%                 keyboard
%             end

        end
    end
    
    COL{i}=temp_col;
    ROW{i}=temp_row;
end
end
%% Refine Lines
if exist('line_data','var') && ~Left
    COL=line_data.COL;
    ROW=line_data.ROW;
    ROW_avg=line_data.ROW_avg;
    line_fig=line_data.line_fig;
    Nlines=length(COL);
end

if line_fig==2  % data_fig will also be 2
    figure
end

% Smooth lines and correct lines
data_fig=figure;
set(data_fig,'Name',num2str(file_name))
button='Correct';
while strcmp(button,'Correct') || strcmp(button,'Crop')
    
%     % Make all lines equal length
%     Llimit=0;
%     Rlimit=Inf;
%     for i=1:Nlines
%         temp_min=min(COL{i});
%         temp_max=max(COL{i});
%         if temp_min>Llimit
%             Llimit=temp_min;
%         end
%         if temp_max<Rlimit
%             Rlimit=temp_max;
%         end
%     end
%     
%     for i=1:Nlines
%         temp_col=COL{i};
%         temp_row=ROW{i};
%         COL{i}=temp_col(temp_col>=Llimit & temp_col<=Rlimit);
%         ROW{i}=temp_row(temp_col>=Llimit & temp_col<=Rlimit);
%     end

    % Remove NaN values
    for i=1:Nlines
        min_ind=find(~isnan(ROW{i}),1,'first');
        max_ind=find(~isnan(ROW{i}),1,'last');
        COL{i}=COL{i}(min_ind:max_ind);
        ROW{i}=ROW{i}(min_ind:max_ind);
    end
    
    % Remove duplicate x-positions (occurs for very steep lines)
    for i=1:Nlines
        num_pts=length(COL{i});
        for j=1:num_pts-1
            if j>=num_pts
                break
            elseif COL{i}(j+1)==COL{i}(j)
                temp_inds=find(COL{i}==COL{i}(j));      % Retain mid index
                temp_ind=floor(mean(temp_inds));
%                 if length(temp_inds)>2
%                     disp('%%%%%%Duplicates')
%                     keyboard
%                 end
                ROW{i}(temp_ind)=mean(ROW{i}(temp_inds));       % Set to average value
                COL{i}(temp_inds(temp_inds~=temp_ind))=[];              % Remove duplicates
                ROW{i}(temp_inds(temp_inds~=temp_ind))=[];
                num_pts=length(COL{i});
%                 if length(temp_inds)>2
%                     keyboard
%                 end
            end
        end
    end
    
    % Smooth lines via moving average
    ROW_avg=cell(size(ROW));
    for i=1:Nlines
        num_pts=length(COL{i});
        vals=ROW{i};    avg_val=nan(length(ROW{i}),1);
        for j=1:num_pts
            
            if j>avg_L/2 && j<length(vals)-avg_L/2+1
                avg_val(j)=nanmean(vals(j-(avg_L-1)/2:j+(avg_L-1)/2));
            end
            
        end
        ROW_avg{i}=avg_val;
        
    end
    
    % Plot smoothed lines
    figure(data_fig)
     imshow(template,'InitialMagnification',65)
     zoom(1.0)
    hold on
    for i=1:Nlines
        plot(COL{i},ROW{i},'r')
        plot(COL{i},ROW_avg{i},'y')
    end
    hold off

    % Correct lines
    %button = questdlg3_custom('Correct Bad Portions of Lines.','Correct Lines','Correct','Crop','Done');
    
    button = questdlg('Correct Bad Portions of Lines.','Correct Lines','Correct','Crop','Done');
    
    switch button
        case 'Correct'
            [corr_x, corr_y] = getline(gcf);
            % Eliminate portion of line in correction region
            temp=[];
            diff=[];
            max_ind=nan(Nlines,1);
            min_ind=nan(Nlines,1);
            for i=1:Nlines
                a = find(COL{i}>=min(corr_x),1,'first');
                b = find(COL{i}<=max(corr_x),1,'last');
                if isempty(a) || isempty(b)
                    continue
                end
                min_ind(i) = a;
                max_ind(i) = b;
                if sum(ROW{i}(min_ind(i):max_ind(i)) >= min(corr_y)-15 & ROW{i}(min_ind(i):max_ind(i)) <= max(corr_y)+15) > 0      % Find line(s) that is inside crop region
                    temp=cat(1,temp,i);
                    [~, temp1]=unique(corr_x);
                    temp2=abs(ROW{i}(min_ind(i):max_ind(i))-interp1(corr_x(temp1),corr_y(temp1),COL{i}(min_ind(i):max_ind(i))));
                    diff=cat(1,diff,sum(temp2(~isnan(temp2))));
                end
            end
            [~,ind]=min(diff);
            line_num=temp(ind);
            ROW{line_num}(min_ind(line_num):max_ind(line_num))=NaN;            
            
            % Create new mask for corrected portion
            % Input corrected x,y positions into data
            figure(line_fig)
            if ~exist('h','var')
                 imshow(template,'InitialMagnification',65)
                msgbox('Draw Dummy Line.','modal');
                h = imline(gca);
            end
                
            mask_lines=false(size(template));
            num_pts=length(corr_x);
            for i=1:num_pts-1
                setPosition(h,corr_x(i:i+1),corr_y(i:i+1))
                temp_mask = createMask(h);
                % Mask for each line
                mask_lines = logical(imadd(mask_lines,temp_mask));
            end
            % x,y pixel locations
            [temp_row,temp_col] = find(mask_lines);
            for i=1:length(temp_col)
                ind=find(COL{line_num}==temp_col(i));
                ROW{line_num}(ind)=temp_row(i);
            end
            
        case 'Crop'
            crop = getrect(gcf);
            % Eliminate portion of line in crop region
            for i=1:Nlines
                min_ind = find(COL{i}>=crop(1),1,'first');
                max_ind = find(COL{i}<=crop(1)+crop(3),1,'last');
                if sum(ROW{i}(min_ind:max_ind) >= crop(2) & ROW{i}(min_ind:max_ind) <= crop(2)+crop(4)) > 0      % Find line that is inside crop region
                    COL{i}(min_ind:max_ind)=NaN;
                    ROW{i}(min_ind:max_ind)=NaN;
                end
            end
    end

end

% Prepare data for storage
clear line_data
line_data.template=template;
line_data.COL=COL;
line_data.ROW=ROW;
line_data.ROW_avg=ROW_avg;
line_data.line_fig=line_fig;
line_data.X=X;
line_data.Y=Y;
line_data.sides=sides;

close(line_fig)

% Save Data with Figure
set(data_fig,'UserData',line_data)

if ~Left
    saveas(data_fig, [fullfile(parent,file_name_) ' lines_' num2str(avg_L) '.fig'])
else
    parent=fullfile(parent,'Left Side');
    if ~exist(parent,'dir')
        mkdir(parent);
    end
    saveas(data_fig, [fullfile(parent,file_name_) ' lines_' num2str(avg_L) '.fig'])
end


table(x_stg_pos,y_stg_pos)
