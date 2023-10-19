clc, close all
clear
%% Loading Images

% NAMING CONVENTION: NAME IMAGE FILES IN NUMERICAL ORDER
% NAME FIXED/UNDEFORMED/ZEmovingRO STRAIN IMAGE BEGINNING WITH 00
% NAME IMAGES IN NUMERICAL ORDER AFTER THAT (01,02,...) IN ORDER FROM
% LEAST DEFORMED TO MOST DEFORMED
% Example: 00_undeformed
        %  01_1 percent strain
        %  02_2 percent strain


path_img = uigetdir; % Select path for files to be taken from
filePattern_img = fullfile(path_img, '*.tif'); % Only allow tif files from path to be selected
Din_img = dir(filePattern_img);
        Ns_img = length(Din_img);
        
        %Prompt user to choose which files to analyze 
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
        w = waitbar(0,'Reading Image Files');
        for i=1:numImg
            Images_projections{i} = imread(fullfile(path_img,Din_img(SelFiles_img(i)).name));
            waitbar(i/numImg)
        end
        close(w)
        


%% Input Parameters

% Define correlation coefficient threshold to identify "bad" regions:
corr_threshold = .5;

% Define subsize for search window (for correlation coeff.):
subSize = 32;

% Define number of neighboring points included in strain calculations:
numP = 12; 

%% Displacement Field Calculations

% Load reference (undeformed) image 
fixed =  Images_projections{1}; % Choose first image from loaded images to be fixed image

for l = 1:numImg-1
    
    % Call DIC Algohrithm for nodal positions and displacements at said
    % nodes. input is the current iteration "l" and the path for images
    if l == 1
        global storage
        storage = [4 4 4 4 4];
    end
    [DICmesh,ResultDisp,RTxA,RTyA]=main_ALDIC(l,path_img);
    
    % Define x and y positions from ALDIC output
    if l == 1
        xposi = DICmesh.coordinatesFEM(:,1);
        yposi = DICmesh.coordinatesFEM(:,2);
    end
    xpos = DICmesh.coordinatesFEM(:,1);
    ypos = DICmesh.coordinatesFEM(:,2);
    dif = xpos(2) - xpos(1);

    % Define x and y displacements from ALDIC output
    xdisp = zeros(length(ResultDisp{1,1}.U)/2,1);
    ydisp = zeros(length(ResultDisp{1,1}.U)/2,1);
    countx = 0;
    county = 0;
    for i = 1 : length(ResultDisp{1,1}.U)
        if round(i/2) == i/2
           county = county + 1;
           ydisp(county,1) = ResultDisp{1,1}.U(i);
        else
           countx = countx + 1;
           xdisp(countx,1) = ResultDisp{1,1}.U(i);
        end
    end

    % Allocate space and define variables
    moving=uint16(nan(max(ypos),max(xpos),numImg-1)); % Allocate space for deformed images
    movingReg=moving; 
    clear moving
    [row, col] = size(movingReg(:,:,1)); % Size of deformed images
    numX = floor((col)/subSize); % Number of windows along x-direction
    numY = floor((row)/subSize); % Number of windows along y-direction
    corrcoeff=nan(numY,numX,numImg-1); % Allocate space for correlation coefficient for each window
    BW = zeros(max(yposi),max(xposi),numImg-1); % Allocate space for mask
    
    % Load deformed image
    moving(:,:,l)=Images_projections{l+1}; 

    % Create meshgrid of x and y vectors
    [X,Y] = meshgrid(min(xpos):dif:max(xpos),min(ypos):dif:max(ypos));
    [Xi,Yi] = meshgrid(min(xposi):dif:max(xposi),min(yposi):dif:max(yposi));
    
    ref = Images_projections{l};
    def = Images_projections{l+1};
    tform = imregcorr(ref,def,'translation');
    RTx = tform.T(3,1);
    RTy = tform.T(3,2);
    % Reshape Dtemp to match meshgrid(temporary displacement field)
    clear Dtemp
    Dtemp(:,:,1)=reshape(xdisp,size(X,2),size(X,1))'+RTx; % x displacements
    Dtemp(:,:,2)=reshape(ydisp,size(X,2),size(X,1))'+RTy; % y displacements

    % When beyond the first iteration, sum the current dataframe with the
    % previous ones.
    if l == 1
        DtempF = Dtemp; 
    else
        DtempP(:,:,1) = interp2(X,Y,Dtemp(:,:,1),Xprime,Yprime,'spline',NaN);
        DtempP(:,:,2) = interp2(X,Y,Dtemp(:,:,2),Xprime,Yprime,'spline',NaN);
        DtempF(:,:,1) = DtempP(:,:,1) + DtempF(:,:,1);
        DtempF(:,:,2) = DtempP(:,:,2) + DtempF(:,:,2);
    end
    
    % Remove NaN Logic. This is done by locating the NaNs in DtempF and
    % summing them along their respective rows and columns. The logic then
    % loops through all of the NaN locations and determines if there is a
    % greater presence (by percentage) in its respective row or column. If
    % far more NaNs are present in the column, it will eliminate the
    % column. If the row and column are within 1% of eachother, it removes
    % the row or column that is farthest from the center of the dataframe
    [rowNaN,colNaN] = find(isnan(DtempF(:,:,1)));
    DtempDF = DtempF;
    [rDtempDF,cDtempDF,zDtempDF] = size(DtempDF);
    rowrem = [];
    colrem = [];
    rcount = 0;
    ccount = 0;
    for jj = 1:length(rowNaN)
        if any(rowrem == rowNaN(jj)) || any(colrem == colNaN(jj))
            continue
        else
            rcount = 0;
            ccount = 0;
            for ii = 1:sum(colNaN(jj)==colNaN)
               clist = find(colNaN(jj) == colNaN); 
               if sum(rowNaN(clist(ii)) == rowrem) == 1
                   ccount = ccount + 0;
               else
                   ccount = ccount + 1;
               end
            end
            colNaNp = ccount / rDtempDF;
            
            for ii = 1:sum(rowNaN(jj)==rowNaN)
               rlist = find(rowNaN(jj) == rowNaN); 
               if sum(colNaN(rlist(ii)) == (colrem)) == 1
                   rcount = rcount + 0;
               else
                   rcount = rcount + 1;
               end
            end
            rowNaNp = rcount / cDtempDF;            
            
            if rowNaNp < .03 && colNaNp < .03
                rowdis = min(rDtempDF - rowNaN(jj), rowNaN(jj));
                coldis = min(cDtempDF - colNaN(jj), colNaN(jj));
                if rowdis < coldis
                    rowrem = cat(2, rowrem, rowNaN(jj));
                elseif coldis < rowdis
                    colrem = cat(2, colrem, colNaN(jj));
                else
                    rowrem = cat(2, rowrem, rowNaN(jj));
                end
            elseif rowNaNp > colNaNp
                rowrem = cat(2, rowrem, rowNaN(jj));
            else
                colrem = cat(2, colrem, colNaN(jj));
            end
        end
    end
    
    % Eliminate spaces in displacements and meshgrid according to NaN
    % removal
    DtempDF(rowrem,:,:) = [];
    DtempDF(:,colrem,:) = [];
    Xn = Xi;
    Yn = Yi;
    Xn(rowrem,:,:) = [];
    Xn(:,colrem,:) = [];
    Yn(rowrem,:,:) = [];
    Yn(:,colrem,:) = [];
    
    vsize = size(DtempDF);
    vrhigh = max(ypos);
    vchigh = max(xpos);
    clear Xq; clear Yq
    
    % Create meshgrid of query points for interpolation
    [Xq,Yq] = meshgrid(min(xposi):max(xposi),min(yposi):max(yposi));
    
    clear DF
    
    % Calculate final displacement field for every pixel (ex: 1:512)
    DF(:,:,1) = interp2(Xn,Yn,DtempDF(:,:,1),Xq,Yq,'spline');
    DF(:,:,2) = interp2(Xn,Yn,DtempDF(:,:,2),Xq,Yq,'spline');  
    
    % Calculate new positions of original meshgrid points in most recent 
    % deformed image(X and Y are original points)
    Xprime=Xi+DtempF(:,:,1);
    Yprime=Yi+DtempF(:,:,2);

    % Show unwarped image based on calculated displacments
    movingReg(min(yposi):max(yposi),min(xposi):max(xposi),l) = imwarp(moving(min(yposi):max(yposi),min(xposi):max(xposi),l),DF); 

    % Show true fixed image
    f2 = figure;
    imshowpair(fixed,movingReg(:,:,l),'Scaling','joint')
    set(gcf,'Name','Comparison')
    %{
    rect=getrect;
    rect = ceil(rect);
    rect(1) = rect(1)+subSize;
    rect(2) = rect(2)+subSize;
    rect(3) = rect(3)-2*subSize;
    rect(4) = rect(4)-2*subSize;
    mask1=false(size(movingReg(:,:,l)));
    mask1(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3))=true;
    %}

    f1 = figure;
    set(gcf,'Name','Bad Regions')
    imshow(movingReg(:,:,l)*2^4) % Nikon images 12 bit, MATLAB forces 16 bit
    title('Define Region of Interest')
    subtitle('Use a 4-point polygon. First point is top left, define clockwise')
    roi = drawpolygon;

    % Create mask to eliminate "bad" regions
    % Compare subwindows using correlation coefficient between original and
    % unwarped images
    for ii = 1:numX
        for j = 1:numY
             windowY = 1+(j-1)*subSize:j*subSize; 
             windowX = 1+(ii-1)*subSize:ii*subSize;
             temp = normxcorr2(fixed(windowY,windowX), movingReg(windowY,windowX,l));
             corrcoeff(j,ii,l) = temp(subSize,subSize);
             % If correlation coefficcient is less than threshold, create
             % box around "bad" region
             if  corrcoeff(j,ii,l) >= corr_threshold
                 % if correlation coefficient is greater than or equal to set threshold, mark as "good"
                 BW(windowY,windowX,l) = 1; 
             else
                  % if correlation coefficient is not greater than or equal to set threshold, mark as "bad"
                  % and create box around this region
                 hold on
                 mask = false(size(movingReg(:,:,l)));
                 mask(windowY,windowX,l) = true;
                 visboundaries(mask(:,:,l),'Color','b');
                 hold off
             end
        end     
    end
    
    mask1=false(size(movingReg(:,:,l)));
    good = 0;
    total = 0;
    m12 = (roi.Position(2,2)-roi.Position(1,2))/(roi.Position(2,1)-roi.Position(1,1));
    b12 = roi.Position(1,2)-m12*roi.Position(1,1);
    m23 = (roi.Position(3,2)-roi.Position(2,2))/(roi.Position(3,1)-roi.Position(2,1));
    b23 = roi.Position(2,2)-m23*roi.Position(2,1);
    m34 = (roi.Position(4,2)-roi.Position(3,2))/(roi.Position(4,1)-roi.Position(3,1));
    b34 = roi.Position(3,2)-m34*roi.Position(3,1);
    m41 = (roi.Position(1,2)-roi.Position(4,2))/(roi.Position(1,1)-roi.Position(4,1));
    b41 = roi.Position(4,2)-m41*roi.Position(4,1);
    for q = 1:max(xposi)
        yb = m12*q+b12;
        yt = m34*q+b34;
        for p = 1:max(yposi)
            if abs(m41) > 1000
                xl = (roi.Position(1,1)+roi.Position(4,1))/2;
            else
                xl = (p-b41)/m41;
            end
            if abs(m23) > 1000
                xr = (roi.Position(2,1)+roi.Position(3,1))/2;
            else
                xr = (p-b23)/m23;
            end
            if xl < q && xr > q
                if yb < p && yt > p
                    accepted = 1;
                else
                    accepted = 0;
                end
            else
                accepted = 0;
            end
            if BW(p,q,l) == 1
                if accepted == 1
                    good = good + 1;
                    total = total + 1;
                    mask1(p,q,l) = 1;
                else %(accepted == 0)
                    BW(p,q,l) = 0;
                    mask1(p,q,l) = 0;
                end
            else %(BW(p,q,l) == 0)
                if accepted == 1
                    total = total + 1;
                    mask1(p,q,l) = 1;
                else %(accepted == 0)
                    mask1(p,q,l) = 0;
                end
            end
        end
    end
    clear roi
    hold on
    visboundaries(mask1(:,:,l),'Color','w')
    hold off
    
    

    %%  Strain Calculations
    
    % Calculating strain only using "good" data points (Bxy_trunc and
    % Cxy_trunc)
    
    % Define positions of search windows (original and temporary)
    if isempty(rowrem) == 0 && isempty(colrem) == 1
        Bxy = [yposi, xposi]; % original x and y positions
        del = find(sum(Bxy(:,1) == rowrem.*dif+(min(yposi)-dif),2));
        del = unique(del);
        Bxy(del,:) = [];
    elseif isempty(rowrem) == 1 && isempty(colrem) == 0
        Bxy = [yposi, xposi]; % original x and y positions
        del = find(sum(Bxy(:,2) == colrem.*dif+(min(xposi)-dif),2));
        del = unique(del);
        Bxy(del,:) = [];
    elseif isempty(rowrem) == 0 && isempty(colrem) == 0
        Bxy = [yposi, xposi]; % original x and y positions
        del = find(sum(Bxy(:,1) == rowrem.*dif+(min(yposi)-dif),2));
        del = cat(1,del,find(sum(Bxy(:,2) == colrem.*dif+(min(xposi)-dif),2)));
        del = unique(del);
        Bxy(del,:) = [];
    else
        Bxy = [yposi, xposi];
    end
    Bxy_trunc = Bxy; % truncated original x and y positions (bad points removed)
    Cxy= [Bxy(:,1) + reshape(DtempDF(:,:,2)',vsize(1) * vsize(2) ,1), Bxy(:,2) + reshape(DtempDF(:,:,1)',vsize(1) * vsize(2),1)]; % new x and y positions
    Cxy_trunc = Cxy; % truncated new x and y positions (bad points removed)


    % Determine if point is within "bad" region; if so, delete it
    % Since bad_inds will change everytime, it must be cleared
    clear bad_inds
    c = 1;
    for m = 1:size(Bxy,1)
        if BW(Bxy(m,1),Bxy(m,2),l) == 0
            bad_inds(c) = m;
            c= c+1;
        end
    end
    

    % Deleting  
    Bxy_trunc(bad_inds,:)=[];
    Cxy_trunc(bad_inds,:)=[];


    % Allocate space for strain variables
    lambdaxx = nan(vsize(1) * vsize(2),1);
    lambdayy = nan(vsize(1) * vsize(2),1);
    epsilonxx = nan(vsize(1) * vsize(2),1);
    epsilonyy = nan(vsize(1) * vsize(2),1);
    gamma = nan(vsize(1)*vsize(2),1);
    maxP = nan(vsize(1) * vsize(2),1);
    minP = nan(vsize(1) * vsize(2),1);
    maxVecU = nan(vsize(1) * vsize(2),1);
    maxVecV = nan(vsize(1) * vsize(2),1);
    minVecU = nan(vsize(1) * vsize(2),1);
    minVecV = nan(vsize(1) * vsize(2),1);
    
    for k=1:size(Bxy,1)
        
        
        % If bad index, strain value is NaN
        if any(bad_inds==k)
            continue
        end
        % Find good points closest to current point for analysis
        distance = sqrt((Bxy(k,1) - Bxy_trunc(:,1)).^2 + (Bxy(k,2) - Bxy_trunc(:,2)).^2);
        [sortD,indsD] = sort(distance);
    
        % Pick numP closest neighboring points for strain calculation
        % If current point is on an edge, leave strain as NaN
        if Bxy(k,1) == min(Bxy_trunc(:,1)) ||  Bxy(k,2) == min(Bxy_trunc(:,2)) || Bxy(k,1) == max(Bxy_trunc(:,1)) || Bxy(k,2) == max(Bxy_trunc(:,2))
            % Top edge
            % inds = [k-1, k+1, i+numX];
            continue
        else
            % Selecting closest point (2) through the numP closest point
            % (numP+1)
            inds=indsD(2:numP+1);
        end
        dX=nan(2,length(inds));
        dx=nan(2,length(inds));
        for ii=1:length(inds)
            % Calculating relative distances between points in original image
            dX(1,ii) = Bxy_trunc(inds(ii),2) - Bxy(k,2);
            dX(2,ii) = Bxy_trunc(inds(ii),1) - Bxy(k,1);
            % Calculating relative distances between points in deformed image
            dx(1,ii) = Cxy_trunc(inds(ii),2) - Cxy(k,2);
            dx(2,ii) = Cxy_trunc(inds(ii),1) - Cxy(k,1);
        end
            % Calculating deformation gradient
            F=dx/dX;                %x = FX + c -> derivative
            C=(transpose(F))*F;
            U=sqrtm(C);             % Removes rigid body rotation
        
            [EigV,~] = eig(U);
            
            Nx=[1; 0];
            Ny=[0; 1];
            % Calculating stretches and strains in x and y directions
            lambdaxx(k)=sqrt((transpose(Nx))*C*Nx);
            epsilonxx(k)=lambdaxx(k)-1;
            lambdayy(k)=sqrt((transpose(Ny))*C*Ny);
            epsilonyy(k)=lambdayy(k)-1;
            gamma(k)=pi/2-acos(C(1,2)/(sqrt(C(1,1))*sqrt(C(2,2))));
            % Calculating principal strains and principal directions
            princLam(1,1) = sqrt((transpose(EigV(:,1)))*C*EigV(:,1));
            princLam(1,2) = sqrt((transpose(EigV(:,2)))*C*EigV(:,2));

            [maxP(k),ind] = max(princLam(1,:));
            maxVecU(k) = EigV(1,ind);
            maxVecV(k) = EigV(2,ind);
            [minP(k),ind] = min(princLam(1,:));
            minVecU(k) = EigV(1,ind);
            minVecV(k) = EigV(2,ind);
            
    end

    % Plotting x strains
    f3 = figure;
    set(gcf,'Name','X Strain')
    temp=reshape(epsilonxx,vsize(2),vsize(1))';
    contourf(Xn,Yn,temp)
    colorbar
    set(gca, 'ydir','reverse')
    
    % Plotting y strains
    f4 = figure;
    set(gcf,'Name','Y Strain')
    temp=reshape(epsilonyy,vsize(2),vsize(1))';
    contourf(Xn,Yn,temp)
    colorbar
    set(gca, 'ydir','reverse')

    % Plotting shear strains
    f5 = figure;
    set(gcf,'Name','Shear Strain')
    temp=reshape(gamma,vsize(2),vsize(1))';   
    contourf(Xn,Yn,temp)
    colorbar
    set(gca, 'ydir','reverse')

    % Plotting max principal strains
    f6 = figure;
    set(gcf,'Name','Max Principal Strain')
    temp=reshape(maxP,vsize(2),vsize(1))'-1;
   
    contourf(Xn,Yn,temp)
    hold on
    XnS = Xn(1:5:vsize(1),1:5:vsize(2));
    YnS = Yn(1:5:vsize(1),1:5:vsize(2)); 
    temp1 = reshape(maxVecU,vsize(2),vsize(1))';
    temp1S = temp1(1:5:vsize(1),1:5:vsize(2));
    temp2 = reshape(maxVecV,vsize(2),vsize(1))';
    temp2S = temp2(1:5:vsize(1),1:5:vsize(2));
    quiver(XnS,YnS,temp1S,temp2S,1,"w",'ShowArrowHead',"off");
    colorbar
    set(gca, 'ydir','reverse')
    hold off

    % Plotting minimum principal strains
    f7 = figure;
    set(gcf,'Name','Min Principal Strain')
    temp=reshape(minP,vsize(2),vsize(1))'-1;

    contourf(Xn,Yn,temp)
    hold on
    temp1 = reshape(minVecU,vsize(2),vsize(1))';
    temp1S = temp1(1:5:vsize(1),1:5:vsize(2));
    temp2 = reshape(minVecV,vsize(2),vsize(1))';
    temp2S = temp2(1:5:vsize(1),1:5:vsize(2));
    quiver(XnS,YnS,temp1S,temp2S,1,"w",'ShowArrowHead',"off");
    colorbar
    set(gca, 'ydir','reverse')
    hold off
    

%%


    % Save plots and values
    global Storage
    if l == 1
        Storage=cell(8,numImg);
        Storage{1,1} = '//////';
        Storage{2,1} = 'Av X Strain';
        Storage{3,1} = 'X Std';
        Storage{4,1} = 'Av Y Strain';
        Storage{5,1} = 'Y Std';
        Storage{6,1} = 'Av Shear Strain';
        Storage{7,1} = 'Shear Std';
        Storage{8,1} = 'Bad Region %';
        for zzz = 1:(numImg-1)
            Storage{1,1+zzz} = sprintf('0-%d',2*zzz);
        end
    end
    Storage{2,l+1} = mean(mean(epsilonxx,'omitnan'),'omitnan');
    Storage{3,l+1} = std(epsilonxx,"omitnan");
    Storage{4,l+1} = mean(mean(epsilonyy,'omitnan'),'omitnan');
    Storage{5,l+1} = std(epsilonyy,"omitnan");
    Storage{6,l+1} = mean(mean(gamma,'omitnan'),'omitnan');
    Storage{7,l+1} = std(gamma,"omitnan");
    Storage{8,l+1} = 100*(total-good)/total;
    
    count = 2*l; % Counter used to name images according to 2% increments
    oldFolder = pwd;
    mkdir NuclearTrackingResults
    cd NuclearTrackingResults
    location = [pwd, '\NuclearTrackingResults'];
    f1name = sprintf('Bad Regions 0-%d',count);
    f2name = sprintf('Comparison 0-%d',count);
    f3name = sprintf('X Strain 0-%d',count);
    f4name = sprintf('Y Strain 0-%d',count);
    f5name = sprintf('Shear Strain 0-%d',count);
    f6name = sprintf('Max Principal Strain 0-%d',count);
    f7name = sprintf('Min Principal Strain 0-%d',count);
    saveas(f1,f1name)
    saveas(f2,f2name)
    saveas(f3,f3name)
    saveas(f4,f4name)
    saveas(f5,f5name)
    saveas(f6,f6name)
    saveas(f7,f7name)
    if l == (numImg - 1)
        writecell(Storage,'Results.xlsx')
    end
    cd(oldFolder)


    keyboard
end
