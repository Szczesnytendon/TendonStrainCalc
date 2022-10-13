%% Calculate Strains from Bleach Lines
clearvars -except ini_path
%close all
scrsz = get(0,'screensize');
labs=gcp('nocreate');
if labs==0
    parpool('MATLAb Parallel Cloud') %parpool('local') 
end

% Pixels over which to average strain or find slope
strain_L=35;
ang_L=35;

% Add spaces for extra plots in contour figs?
addPlots=0;

% Code for ajusting strain contours
% a=get(gcf,'Position');
% a(2)=a(2)+40;
% set(gcf,'Position',a)

% Size for avg strain plot
% a=[ 590   300   230   430];

%% Select Figures for Analysis
if ~exist('ini_path','var') || max(ini_path==0)
    ini_path='C:\Users\bep15\Box\b-szczesny-lab Shared\Student Folders\Ben\Data\Abstract\Test7_E20\Tif images\0%\Center';
end
[file,ini_path,type] = uigetfile({'*.fig','Figure(*.fig)'},'Choose Figure',ini_path);
path = fileparts(fullfile(ini_path,file));
parent=path;

% Find files at same location
ind=regexp(file,'\d');
max_ind=min(ind)-1;
prefixO=file(1:max_ind);

% Find files with same avgeraging length
max_ind=find(file=='_',1,'last');
if ~isempty(max_ind)
    suffix=file(max_ind:end);
    avg_L=suffix(2:end-4);
else
    suffix='.fig';
end

Din = dir(fullfile(parent,[prefixO '*' suffix]));
Ns = length(Din);

%Prompt user to choose which files to analyze
SelFiles = listdlg(...
    'PromptString', 'Choose Specific Files to Analyze', ...
    'SelectionMode', 'Multiple', ...
    'Name', 'File List', ...
    'InitialValue', 1:Ns, ...
    'ListString', {Din.name},...
    'Listsize',[300 400]);

numF = length(SelFiles);

%% Read and Initialize Data

% Sort Figures by Applied Strain Magnitudes
Nom_strain=nan(numF,1);
for i=1:numF
    temp=Din(SelFiles(i)).name;
    max_ind=find(temp=='_',1,'first');
    if ~isempty(max_ind)
        ind=regexp(temp(1:max_ind),'\d');
    else
        ind=regexp(temp,'\d');
    end
    Nom_strain(i)=str2double(temp(ind));
end
[dum, inds]=sort(Nom_strain);

names=cell(numF,1);
COL_orig=cell(numF,1);
ROW_orig=cell(numF,1);
COL_rot=cell(numF,1);
ROW_rot=cell(numF,1);
center_x=nan(numF,1);
center_y=nan(numF,1);

for fig=1:numF
    % Assign names
    temp=Din(SelFiles(inds(fig))).name;
    names{fig}=temp(1:find(temp==' ',1,'last')-1);

    % Open each figure
    open(fullfile(parent,Din(SelFiles(inds(fig))).name))
    line_data=get(gcf,'UserData');
    close(gcf)
    if fig==1
        Isize=size(line_data.template,1);
        template=uint16(zeros(Isize,Isize,numF));
    end
    
    % Initialize variables
    template(:,:,fig)=line_data.template;
    COL=line_data.COL;
    ROW=line_data.ROW;
    ROW_avg=line_data.ROW_avg;
    
    temp1=regexp(names{fig},'ini', 'once');
    temp2=regexp(names{fig},'end', 'once');
    if ~isempty(temp1)
        prefix=['INI___' prefixO];
    elseif ~isempty(temp2)
        prefix=['END___' prefixO];
    else
        prefix=prefixO;
    end
    
%   Remove certain lines if desired
%     remL=[3 4];
% %     if exist('remL','var')
% %         prefix=['DEL___' prefixO];
% %     end
%     COL(remL)=[];
%     ROW(remL)=[];
%     ROW_avg(remL)=[];

    Nlines=length(COL);
    % Remove NaN values
    COL_avg=COL;
    for i=1:Nlines
        min_ind=find(~isnan(ROW_avg{i}),1,'first');
        max_ind=find(~isnan(ROW_avg{i}),1,'last');
        COL_avg{i}=COL_avg{i}(min_ind:max_ind);
        ROW_avg{i}=ROW_avg{i}(min_ind:max_ind);
    end
    
    %%% Rotate FOV to make Reference Lines Horizontal
    if fig==1
        % Calculate rotation tensor
        slope0=nan(Nlines,1);
        ang0=slope0;
        for i=1:Nlines
            [coeffs, gof] = fit(COL_avg{i},ROW_avg{i},'poly1');
            slope0(i)=coeffs.p1;
            ang0(i)=atand(slope0(i));
            
            % figure(gcf)
            % hold on
            % plot(COL_avg0{i},feval(coeffs,COL_avg{i}))
            % hold off
        end
        Q=[cosd(mean(ang0)) sind(mean(ang0)); -sind(mean(ang0)) cosd(mean(ang0))];
    end
    
    % Perform rotation
    temp_colrot=COL_avg;
    temp_rowrot=ROW_avg;
    temp_colorig=COL;
    temp_roworig=ROW;
    center_x(fig)=mean([max(COL_avg{1}),min(COL_avg{1})]);
    center_y(fig)=mean([mean(ROW_avg{Nlines}),mean(ROW_avg{1})]);
    
    for i=1:Nlines
        X=cat(2,COL_avg{i}-center_x(fig),ROW_avg{i}-center_y(fig))';    % Rotate about image center (just for visualization purposes; center of rotation doesn't affect strains or angles)
        temp=Q*X;
        temp_pt=temp(1,1)+center_x(fig);
        if (round(temp_pt) - temp_pt) > 0 % Use rounding direction of left most pixel
            temp_colrot{i}=ceil(temp(1,:)'+center_x(fig));  % Align x-pos to nearest pixel for simplicity; not concerned with lateral deformations
        else
            temp_colrot{i}=floor(temp(1,:)'+center_x(fig));
        end
        temp_rowrot{i}=temp(2,:)'+center_y(fig);
        
        % Raw line points
        Xorig=cat(2,COL{i}-center_x(fig),ROW{i}-center_y(fig))';    % Rotate about image center (just for visualization purposes; center of rotation doesn't affect strains or angles)
        temp_orig=Q*Xorig;
        temp_colorig{i}=temp_orig(1,:)'+center_x(fig);
        temp_roworig{i}=temp_orig(2,:)'+center_y(fig);

        %     [coeffs, gof] = fit(temp_colrot{i},temp_rowrot{i},'poly1');
        %     figure(gcf)
        %     hold on
        %     plot(temp_colrot{i},feval(coeffs,temp_colrot{i}),'r')
        %     hold off
    end
    COL_orig{fig}=temp_colorig;
    ROW_orig{fig}=temp_roworig;
    
    % Remove duplicate x-positions (occurs for very steep lines)
    for i=1:Nlines
        num_pts=length(temp_colrot{i});
        for j=1:num_pts-1
            if j>=num_pts
                break
            elseif temp_colrot{i}(j+1)==temp_colrot{i}(j)
                temp_inds=find(temp_colrot{i}==temp_colrot{i}(j));      % Retain mid index
                temp_ind=floor(mean(temp_inds));
%                 if length(temp_inds)>2
%                     disp('%%%%%%Duplicates')
%                     keyboard
%                 end
                temp_rowrot{i}(temp_ind)=mean(temp_rowrot{i}(temp_inds));       % Set to average value
                temp_colrot{i}(temp_inds(temp_inds~=temp_ind))=[];              % Remove duplicates
                temp_rowrot{i}(temp_inds(temp_inds~=temp_ind))=[];
                num_pts=length(temp_colrot{i});
%                 if length(temp_inds)>2
%                     keyboard
%                 end
            end
        end
    end
    
    % Add missing x-positions
    for i=1:Nlines
        num_pts=length(temp_colrot{i});
        j=1;
        while j<num_pts
            if j>=num_pts
                break
            elseif temp_colrot{i}(j+1)~=temp_colrot{i}(j)+1
                ins_col=(temp_colrot{i}(j)+1:temp_colrot{i}(j+1)-1)';
                ins_row=interp1(temp_colrot{i}(j:j+1),temp_rowrot{i}(j:j+1),ins_col);       % Interpolate from bounding points
%                 disp('%%%%%%Missing')
%                 keyboard
                temp_colrot{i}=cat(1,temp_colrot{i}(1:j),ins_col,temp_colrot{i}(j+1:end));
                temp_rowrot{i}=cat(1,temp_rowrot{i}(1:j),ins_row,temp_rowrot{i}(j+1:end));
                num_pts=length(temp_colrot{i});
%                 keyboard
            end
            j=j+1;
        end
    end
    
    % Make all lines equal length
    Llimit=0;
    Rlimit=Inf;
    for i=1:Nlines
        temp_min=min(temp_colrot{i});
        temp_max=max(temp_colrot{i});
        if temp_min>Llimit
            Llimit=temp_min;
        end
        if temp_max<Rlimit
            Rlimit=temp_max;
        end
    end
    
    TEMP_col=temp_colrot;
    TEMP_row=temp_rowrot;
    temp_colrot=nan(Nlines,Rlimit-Llimit+1);
    temp_rowrot=nan(Nlines,Rlimit-Llimit+1);
    sum_error=0;
    for i=1:Nlines
        temp_col=TEMP_col{i};
        temp_row=TEMP_row{i};
        temp_colrot(i,:)=temp_col(temp_col>=Llimit & temp_col<=Rlimit);
        temp_rowrot(i,:)=temp_row(temp_col>=Llimit & temp_col<=Rlimit);
        sum_error=sum_error + sum(temp_colrot(1,:)-temp_colrot(i,:));
    end
    
    if sum_error~=0 || sum(temp_colrot(1,:)-(Llimit:Rlimit))~=0
        errordlg('COL_rot has missing or duplicated points')
        keyboard
    end
    COL_rot{fig}=Llimit:1:Rlimit;
    ROW_rot{fig}=temp_rowrot;
    
%     keyboard

%         figure(gcf)
%         hold on
%         plot(COL_rot{fig},ROW_rot{fig},'y')
%         plot(center_x,center_y,'r+')
%         hold off
%         keyboard
%         close(gcf)
end

%% Calculate Variables

% Smooth lines via moving average
ang=cell(numF,1);
DIST0=nan(Nlines-1,1);
DIST=cell(numF,1);
strain=cell(numF,1);
strain_max=-Inf;
strain_min=Inf;
ang_max=-Inf;
ang_min=Inf;
max_pts=0;

for fig=1:numF
    num_pts=length(COL_rot{fig});
    if num_pts>max_pts
        max_pts=num_pts;
    end
    
    % Calculate distance between each pair of bleach lines
    temp=nan(Nlines-1,num_pts);
    for i=1:Nlines-1
        y_posT=ROW_rot{fig}(i,:);       % Top row
        y_posB=ROW_rot{fig}(i+1,:);     % Bottom row
    
        for j=1:num_pts
%             if j>strain_L/2 && j<num_pts-strain_L/2+1
%                 temp(i,j)=mean(y_posB(j-(strain_L-1)/2:j+(strain_L-1)/2)) - mean(y_posT(j-(strain_L-1)/2:j+(strain_L-1)/2));
                temp(i,j)=y_posB(j) - y_posT(j);
%             end
        end
        
        % Initial distance between each line pair
        if fig==1
            DIST0(i)=nanmean(temp(i,:));
        end
        
    end
    DIST{fig}=temp;         % Raw unaveraged pixel-wise distances
    
    temp=nan(Nlines-1,num_pts);
    % Calculate strains
    for i=1:Nlines-1
        for j=1:num_pts
            if j>strain_L/2 && j<num_pts-strain_L/2+1
%                 temp(i,j)=(DIST{fig}(i,j)-DIST0(i))/DIST0(i);
                temp(i,j)=(mean(DIST{fig}(i,j-(strain_L-1)/2:j+(strain_L-1)/2))-DIST0(i))/DIST0(i);     % Strain smoothed by averaging distances in vicintiy of pixel
            end
        end
    end
    strain{fig}=temp;
    temp=max(max(strain{fig}));
    if temp > strain_max
        strain_max=temp;
    end
    temp=min(min(strain{fig}));
    if temp < strain_min
        strain_min=temp;
    end
end

w = waitbar(0,['Calculating Angles for Fig 1 of ' num2str(numF)]);
for fig=1:numF
    waitbar(0,w,['Calculating Angles for Fig ' num2str(fig) ' of ' num2str(numF)])
    num_pts=length(COL_rot{fig});
    % Calculate angles
    temp_ang=nan(Nlines,num_pts);
    for i=1:Nlines
        x_pos=COL_rot{fig};
        y_pos=ROW_rot{fig}(i,:);
        
        parfor j=1:num_pts
            if j>ang_L/2 && j<num_pts-ang_L/2+1
                [coeffs, gof] = fit(x_pos(j-(ang_L-1)/2:j+(ang_L-1)/2)',y_pos(j-(ang_L-1)/2:j+(ang_L-1)/2)','poly1');
            temp=coeffs.p1;
            temp_ang(i,j)=-atand(temp);     % Define pos angle as CCW
            end
        end
        waitbar(i/Nlines,w)
    end
    ang{fig}=temp_ang;
    temp=max(max(ang{fig}));
    if temp > ang_max
        ang_max=temp;
    end
    temp=min(min(ang{fig}));
    if temp < ang_min
        ang_min=temp;
    end
    
end
close(w)
%% Summary Variables

WinSize=floor(Isize/1000)/2*1000;

Nom_strain=sort(Nom_strain,'descend');
appS=cell(numF,1);
for fig=1:numF
    appS{fig}=[num2str(Nom_strain(fig)) '%'];
end

Slabel=cell(Nlines-1,1);
Slabel{1}='Top';    Slabel{Nlines-1}='Bottom';
for i=2:length(Slabel)-1
   Slabel{i}=''; 
end

%%%%%%%%%%%%%%%%%%% Average strains%%%%%%%%%%%%%%%%%%%%
avgS_horz=nan(Nlines-1,numF);
avgS_vert=nan(numF,max_pts);
x_pos=nan(numF,max_pts);
for fig=1:numF
    avgS_horz(:,fig)=nanmean(strain{fig},2);       % Average across width between pairs of bleach lines (i=1 is top-most pair)
    avgS_vert(fig,1:size(strain{fig},2))=nanmean(strain{fig},1);   % Average across tendon length
    x_pos(fig,1:length(COL_rot{fig}))=COL_rot{fig};%-center_x(fig);
end

figure(1)
set(gcf,'Position',[100 100 .5*scrsz(3) .8*scrsz(4)],'Name',['Avg Strains, Avg_L ' avg_L])
movegui(gcf,'center')
c_lines=colormap(lines);
% Average across width
subplot(2,1,1)
set(gca,'ColorOrder',flipud(c_lines(1:numF,:)),'NextPlot','replacechildren')
plot(fliplr(avgS_horz),(1:Nlines-1)')
set(gca,'YLim',[0 Nlines],'YTick',1:Nlines-1,'YTickLabel',Slabel,'YDir','reverse')
title('Strain Averaged Across Tendon Width')
xlabel('Average Strain')
ylabel('Tendon Position')

for fig=1:numF
    if rem(fig,2)   % for odd numbers
        text(avgS_horz(1,fig),.6,appS(end-fig+1))
    else
        text(avgS_horz(1,fig),.3,appS(end-fig+1))
    end
end

% Average across length
subplot(2,1,2)
set(gca,'ColorOrder',flipud(c_lines(1:numF,:)),'NextPlot','replacechildren')
plot(fliplr(x_pos'),fliplr(avgS_vert'))
% set(gca,'XLim',[-WinSize WinSize],'XTick',[-WinSize WinSize]*.8,'XTickLabel',{'Left';'Right'})
set(gca,'XLim',[1 Isize],'XTick',[Isize*(1/.8-1) Isize]*.8,'XTickLabel',{'Left';'Right'})
title('Strain Averaged Across Tendon Length')
ylabel('Average Strain')
xlabel('Tendon Position')
legend(appS,'location','NorthEast')

Data.avgH=avgS_horz;
Data.avgV=avgS_vert;
Data.xpos=x_pos;
set(gcf,'UserData',Data)
clear Data
saveas(gcf,fullfile(parent,[prefix 'Avg Strains, Avg_L ' avg_L ', Strain_L ' num2str(strain_L) '.fig']))

% All Strains Grouped by Bleach Line Pairs
figure(2)
set(gcf,'Position',[100 100 .5*scrsz(3) .8*scrsz(4)],'Name',['Srain Plots, Avg_L ' avg_L])
movegui(gcf,'center')
c_lines=colormap(lines);

for i=1:Nlines-1
    subplot(2,2,i)
    set(gca,'ColorOrder',flipud(c_lines(1:numF,:)),'NextPlot','replacechildren')
    temp=nan(numF,max_pts);
    for fig=1:numF
        temp(fig,1:size(strain{fig},2))=strain{fig}(i,:);
        x_pos(fig,1:length(COL_rot{fig}))=COL_rot{fig};%-center_x(fig);
    end
    plot(fliplr(x_pos'),fliplr(temp'))
%     set(gca,'XLim',[-WinSize WinSize],'XTick',[-WinSize WinSize]*.8,'XTickLabel',{'Left';'Right'})
    set(gca,'XLim',[1 Isize],'XTick',[Isize*(1/.8-1) Isize]*.8,'XTickLabel',{'Left';'Right'})
    set(gca,'YLim',[floor(strain_min*100)/100 ceil(strain_max*100)/100])
    title(['Raw Strain for Bleach Line Pair ' num2str(i)])
    ylabel('Strain')
    xlabel('Tendon Position')
    if i==2
        legend(appS,'location','NorthEast','Position',[0.8970 0.7843 0.0875 0.1350]);
    end
end
saveas(gcf,fullfile(parent,[prefix 'Strain Plots, Avg_L ' avg_L ', Strain_L ' num2str(strain_L) '.fig']))

% Find strain errors (non-zeros vals) for reference image
temp=strain{1}.^2;
temp=nanmean(temp,2);
RMSerrorS=sqrt(temp);

temp=ang{1}.^2;
temp=nanmean(temp,2);
RMSerrorA=sqrt(temp);

%%%%%%%%%%%%%%%%%%%%%%%%%% Representative Angles%%%%%%%%%%%%%%%%%%%%%%%%%%%
RMSang=nan(Nlines,numF);
ang_gen=nan(Nlines,numF);
avgA_vert=nan(numF,max_pts);

Alabel=cell(Nlines,1);
Alabel{1}='Top';    Alabel{Nlines}='Bottom';
for i=2:length(Alabel)-1
   Alabel{i}=''; 
end

for fig=1:numF
    temp_ang=ang{fig}.^2;
    temp_ang=nanmean(temp_ang,2);
    RMSang(:,fig)=sqrt(temp_ang);       % RMS angle of each bleach line
    temp_max=max(max(RMSang));
    
    for i=1:Nlines
        [coeffs, gof] = fit(COL_rot{fig}',ROW_rot{fig}(i,:)','poly1');
        temp=coeffs.p1;
        ang_gen(i,fig)=-atand(temp);     % Overall angle of each bleach line (again, pos is CCW)
        if max(max(ang_gen))>temp_max
            temp_max=max(max(ang_gen));
        end
        
    end
    
    avgA_vert(fig,1:size(ang{fig},2))=mean(ang{fig},1);   % Average across tendon length
end

figure(3)
set(gcf,'Position',[100 100 .5*scrsz(3) .8*scrsz(4)],'Name',['Rep Angles, Avg_L ' avg_L ', Ang_L ' num2str(ang_L)])
movegui(gcf,'center')
c_lines=colormap(lines);
subplot(2,1,1)
% Overall angles
set(gca,'ColorOrder',flipud(c_lines(1:numF,:)),'NextPlot','replacechildren')
plot(fliplr(ang_gen),(1:Nlines)')
set(gca,'YLim',[0 Nlines+1],'YTick',1:Nlines,'YTickLabel',Alabel,'YDir','reverse','XLim',[floor(min(min(ang_gen))) ceil(temp_max)])
title('Overall Angle of Each Bleach Line')
xlabel('Angle (deg)')
ylabel('Tendon Position')

for fig=1:numF
    if rem(fig,2)   % for odd numbers
        text(ang_gen(1,fig),.6,appS(end-fig+1))
    else
        text(ang_gen(1,fig),.3,appS(end-fig+1))
    end
end

% RMS Angles
subplot(2,1,2)
set(gca,'ColorOrder',flipud(c_lines(1:numF,:)),'NextPlot','replacechildren')
plot(fliplr(RMSang),(1:Nlines)')
set(gca,'YLim',[0 Nlines+1],'YTick',1:Nlines,'YTickLabel',Alabel,'YDir','reverse','XLim',[floor(min(min(ang_gen))) ceil(temp_max)])
title('RMS Angle of Each Bleach Line')
xlabel('RMS Angle (deg)')
ylabel('Tendon Position')

for fig=1:numF
    if rem(fig,2)   % for odd numbers
        text(RMSang(1,fig),.6,appS(end-fig+1))
    else
        text(RMSang(1,fig),.3,appS(end-fig+1))
    end
end

Data.ang_gen=ang_gen;
Data.RMSang=RMSang;
set(gcf,'UserData',Data)
clear Data
saveas(gcf,fullfile(parent,[prefix 'Rep Angles, Avg_L ' avg_L ', Ang_L ' num2str(ang_L) '.fig']))

%%%%%%%%%%%%% Average Strain and Angle across length%%%%%%%%%%%%%%%%%%%%
figure(4)
set(gcf,'Position',[100 100 .5*scrsz(3) .8*scrsz(4)],'Name',['Both, Avg_L ' avg_L ', Ang_L ' num2str(ang_L)])
movegui(gcf,'center')
c_lines=colormap(lines);
% Strain
subplot(2,1,1)
set(gca,'ColorOrder',flipud(c_lines(1:numF,:)),'NextPlot','replacechildren')
plot(fliplr(x_pos'),fliplr(avgS_vert'))
% set(gca,'XLim',[-WinSize WinSize],'XTick',[-WinSize WinSize]*.8,'XTickLabel',{'Left';'Right'})
set(gca,'XLim',[1 Isize],'XTick',[Isize*(1/.8-1) Isize]*.8,'XTickLabel',{'Left';'Right'})
title('Strain Averaged Across Tendon Length')
ylabel('Average Strain')
xlabel('Tendon Position')
legend(appS,'location','NorthEast')

% Angles
subplot(2,1,2)
set(gca,'ColorOrder',flipud(c_lines(1:numF,:)),'NextPlot','replacechildren')
plot(fliplr(x_pos'),fliplr(avgA_vert'))
% set(gca,'XLim',[-WinSize WinSize],'XTick',[-WinSize WinSize]*.8,'XTickLabel',{'Left';'Right'})
set(gca,'XLim',[1 Isize],'XTick',[Isize*(1/.8-1) Isize]*.8,'XTickLabel',{'Left';'Right'})
title('Angles Averaged Across Tendon Length')
ylabel('Average Angle (deg)')
xlabel('Tendon Position')
legend(appS,'location','NorthEast')

Data.avg_angV=avgA_vert;
set(gcf,'UserData',Data)
clear Data
saveas(gcf,fullfile(parent,[prefix 'Avg Strains and Angles, Avg_L ' avg_L ', Strain_L ' num2str(strain_L) ', Ang_L ' num2str(ang_L) '.fig']))

%%%%%%%% Angles for Each Line and Average for Each Figure %%%%%%%
LegText=cell(Nlines+1,1);
for i=1:Nlines
    LegText{i}=['Line ' num2str(i)];
end
LegText{Nlines+1}='Average';

% All Angles and Average for each Figure
figure(5)
set(gcf,'Position',[100 100 .85*scrsz(3) .8*scrsz(4)],'Name','Angles for Each Bleach Line')
movegui(gcf,'center')

RMSang2=nan(1,numF);
SDevAng=nan(1,numF);
for fig=1:numF
    temp_ang=ang{fig};
    temp_avg=mean(temp_ang,1);
    RMSang2(fig)=sqrt(nanmean(temp_avg.^2));
    SDevAng(fig)=std(temp_avg(~isnan(temp_avg)),1);
    subplot(3,ceil(numF/3),fig)
%     plot((COL_rot{fig}-center_x(fig))',temp_ang')
    plot((COL_rot{fig})',temp_ang')
    hold on
%     plot((COL_rot{fig}-center_x(fig))',temp_avg','k','LineWidth',3)
    plot((COL_rot{fig})',temp_avg','k','LineWidth',3)
    hold off
%     set(gca,'XLim',[-WinSize WinSize],'XTick',[-WinSize WinSize]*.8,'XTickLabel',{'Left';'Right'})
    set(gca,'XLim',[1 Isize],'XTick',[Isize*(1/.8-1) Isize]*.8,'XTickLabel',{'Left';'Right'})
    set(gca,'YLim',[floor(ang_min/10)*10 ceil(ang_max/10)*10])
    title(['Angles for Each Bleach Line and Avg Line - ' appS{end-fig+1}])
    ylabel('Angle (deg)')
    xlabel('Tendon Position')
    text(.05,.1,['Sdev = ' sprintf('%.1f', SDevAng(fig)) ' deg'],'Units','normalized','Color',[0 0 0])
    if fig==2
        legend(LegText,'location','NorthEast','Position',[0.8970 0.7843 0.0875 0.1350]);
    end
end
Data.RMSang2=RMSang2;
Data.SDevAng=SDevAng;
set(gcf,'UserData',Data)
clear Data
saveas(gcf,fullfile(parent,[prefix 'Angles per Line, Avg_L ' avg_L ', Ang_L ' num2str(ang_L) '.fig']))

%% Plot Bleach Lines from All Figures Together
% figure
% c=colormap(lines);
% hLines=cell(numF,1);
% hGroups=cell(numF,1);
% imshow(template(:,:,1))
% hold on
% for fig=1:numF
%     offset_x=center_x(fig)-center_x(1);
%     offset_y=center_y(fig)-center_y(1);
%     hLines{fig}=plot(COL_rot{fig}-offset_x,ROW_rot{fig}-offset_y,'Color',c(fig,:));
%     hGroups{fig}=hggroup;
%     set(hLines{fig},'Parent',hGroups{fig})
%     set(get(get(hGroups{fig},'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','on'); 
% end
% legend(names)
% hold off

%% Plot Strain and Angle Contours
figS=6;
numPlots=numF+addPlots;

for k=figS:figS+1     %Strain=figS and Angle=figS+1
    figure(k)
    if k==figS
        CLim=[strain_min strain_max];
        FigName=['Strain Contours, Smoothing Length=' avg_L ', Strain Bin L=' num2str(strain_L) ', Avg RMS Error=' num2str(mean(RMSerrorS))];
    else
        CLim=[ang_min ang_max];
        FigName=['Angle Contours, Smoothing Length=' avg_L ', Slope Length=' num2str(ang_L) ', Avg RMS Error=' num2str(RMSang2(1))];
    end
    
    % Set figure and subplot axis sizes
    fig_w=scrsz(3)-20;
    margin=10; %.5;
    plot_w=(fig_w - 100 - margin*numPlots)/(fig_w*numPlots);
    fig_h=1.5*plot_w*scrsz(3);
    set(gcf,'Position',[10   150   fig_w   fig_h],'Name',FigName)
    movegui(gcf,'center')
    margin= margin/fig_w;
    
    for fig=1:numF
        % Assign strain/angle values to all points between bleach lines at each
        % x-pos
        background=cat(3,template(:,:,fig),template(:,:,fig),template(:,:,fig));
        x_pos=COL_rot{fig};
        num_pts=length(COL_rot{fig});
        
        if k==figS
            
            strainC=nan(Isize);
            for i=1:Nlines-1
                y_posT=round(ROW_rot{fig}(i,:));       % Top row
                y_posB=round(ROW_rot{fig}(i+1,:));     % Bottom row
                
                for j=1:num_pts
                    strainC(y_posT(j):y_posB(j),x_pos(j))=strain{fig}(i,j);     % Note x-dir is 2nd dimension
                end
            end
            
        else
            
            angC=nan(Isize);
            for i=1:Nlines
                if i==1
                    y_posT=round(ROW_rot{fig}(i,:));       % Top row
                    y_posB=round(ROW_rot{fig}(i,:) + DIST{fig}(1,:)/2);     % Bottom row
                elseif i==Nlines
                    y_posT=round(ROW_rot{fig}(i,:) - DIST{fig}(Nlines-1,:)/2);       % Top row
                    y_posB=round(ROW_rot{fig}(i,:));     % Bottom row
                else
                    y_posT=round(ROW_rot{fig}(i,:) - DIST{fig}(i-1,:)/2);       % Top row
                    y_posB=round(ROW_rot{fig}(i,:) + DIST{fig}(i,:)/2);     % Bottom row
                end
                
                for j=1:num_pts
                    angC(y_posT(j):y_posB(j),x_pos(j))=ang{fig}(i,j);     % Note x-dir is 2nd dimension
                end
            end
            
        end
        
        subplot(1,numPlots,fig)
        set(gca,'Position',[(fig*margin+(fig-1)*plot_w) .1 .95*plot_w .8])
        imagesc(background,CLim)
        axis image
        hold on
        if k==figS
            contourf(strainC,12)
            text(.05,.05,['Avg Strain = ' sprintf('%.4f', mean(avgS_horz(:,fig)))],'Units','normalized','Color',[0 0 0],'BackgroundColor',[1 1 1])
        else
             contourf(angC,12)
             text(.05,.05,sprintf('Sdev = %.1f deg', SDevAng(fig)),'Units','normalized','Color',[0 0 0],'BackgroundColor',[1 1 1])
%              text(.05,.05,sprintf('Avg RMS = %.1f deg\nSdev = %.1f deg', [RMSang2(fig) SDevAng(fig)]),'Units','normalized','Color',[0 0 0],'BackgroundColor',[1 1 1])
        end
        shading flat
        plot(COL_rot{fig},ROW_rot{fig},'w','LineWidth',4);
        col_orig=COL_orig{fig};
        row_orig=ROW_orig{fig};
        for i=1:Nlines
            plot(col_orig{i},row_orig{i},'k','LineWidth',2)
        end
        hold off
        set(gca,'XTick',[],'YTick',[])
        title(names{fig})
    end

    % Create separate colorbar
    axes('Position', [0.95 0.05 0.1 0.9], 'Visible', 'off','CLim',CLim);
    c=colorbar('West');
    
    if k==figS
        ylabel(c,'Strain')
        Data.strain=strain;
        Data.xpos=x_pos;
        Data.RMS=RMSerrorS;
        set(gcf,'UserData',Data)
        clear Data
        saveas(gcf,fullfile(parent,[prefix 'Strain Contours, Avg_L ' avg_L ', Strain_L ' num2str(strain_L) '.fig']))
    else
        ylabel(c,'Angle (deg)')
        Data.ang=ang;
        Data.xpos=x_pos;
        Data.RMS=RMSang2;
        set(gcf,'UserData',Data)
        clear Data
        saveas(gcf,fullfile(parent,[prefix 'Angle Contours, Avg_L ' avg_L ', Ang_L ' num2str(ang_L) '.fig']))
    end
end

%% Add subplots

for fig=numF+1:numPlots
    
    [file,path,type] = uigetfile({'*.fig','Figure(*.fig)'},'Choose Figure',ini_path);
    open(fullfile(path,file))
    temp=get(gcf,'UserData');
    close(gcf)
    pic=temp.template;
    background=cat(3,pic,pic,pic);
    name=file(1:find(file=='l')-1);
    
    for k=figS:figS+1     %Strain=figS and Angle=figS+1
        figure(k)
        subplot(1,numPlots,fig)
        set(gca,'Position',[(fig*margin+(fig-1)*plot_w) .1 .95*plot_w .8])
        
        if k==figS
            CLim=[strain_min strain_max];
        else
            CLim=[ang_min ang_max];
        end

        imagesc(background,CLim)
        axis image
        set(gca,'XTick',[],'YTick',[])
        title(name)
        
        if k==figS
            saveas(gcf,fullfile(parent,[prefix 'Strain Contours, Avg_L ' avg_L ', Strain_L ' num2str(strain_L) '.fig']))
        else
            saveas(gcf,fullfile(parent,[prefix 'Angle Contours, Avg_L ' avg_L ', Ang_L ' num2str(ang_L) '.fig']))
        end
        
    end

end