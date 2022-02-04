%% Dutta.D(1), Misra.S.(1), and Karmakar.S.(2)
% Deformation mechanisms and characteristics of the meta-BIFs from an early
% Proterozoic shear system of the Southern Granulite Terrane (SGT), India
% Submitted to Journal of Structural Geology
% (1) Experimental Rock Deformation Laboratory, Department of Earth
%     Sciences, IIT Kanpur
% (2) Department of Geological Sciences, Jadavpur University

% Please contact the first author at dripta.dutta@gmail.com for queries

%% NOTE:

% The following script aims to automate the subgrain boundary trace
% analysis (e.g. Prior et al. 2002, Reddy et al. 2007, Barrie et al. 2008)
% and plot the orientations of the subgrain boundary wall and slip plane 
% along with their poles as equal area, lower hemisphere stereoplots.
% It also plots the dispersion axis for a selected subset of the
% EBSD/orientation data using the principal geodesic analysis of Michels
% et al. (2015).

% Users are advised to clean the EBSD data and reconstruct the grains
% beforehand. Please also download the PGA code from the following link:
% https://github.com/zmichels/CVA/blob/master/PGA.m

% This script identifies the nature of the subgrain boundaries i.e., tilt,
% tilt/twist, and calculates the orientations of the subgrain boundary wall
% and the slip plane (and their poles) accordingly using the orientation
% of the misorientation axes maximum and the mean trend of the selected
% subgrain boundary traces. Please see the text in Sections 3.2.2 and
% 3.2.3, and Fig. 2 of the manuscript for more details.

% Steps in brief:
% -> Select the EBSD data and grains of the phase of interest and set
%    the lower angular threshold for the subgrain boundaries.
% -> Plot the misorientation to mean orientation (mis2mean) map.
% -> Select the grain of interest.
% -> Plot the orientation map of the grain selected in the previous step
%    along with the subgrain boundary traces. The traces are colorized
%    according to the angle (sgb_ang) between them and the projection
%    of the corresponding misorientation axes on the plane of observation
%    i.e., the XZ-section.
% -> Select similarly oriented subgrain boundary traces using the polygon
%    feature. Save the polygon variable.
% -> Derive the misorientation axes maximum.
% -> Choose the miller indices and calculate their orientations
%    corresponding to the grain mean orientation.
% -> Categorize boundaries as tilt ONLY if sgb_ang<15 degrees, and as
%    tilt/twist boundaries if sgb_ang>75 degrees.
% -> Plot and save the resulting stereoplots

% The MTEX documentation on tilt and twist boundaries and discussions
% on the MTEX forum were helpful.

% Full references of the articles cited here:
% Barrie, C. D., Boyle, A. P., Cox, S. F., & Prior, D. J. (2008). 
%   Slip systems and critical resolved shear stress in pyrite: an electron
%   backscatter diffraction (EBSD) investigation. Mineralogical Magazine,72

% Michels, Z. D., Kruckenberg, S. C., Davis, J. R., & Tikoff, B. (2015).
%    Determining vorticity axes from grain-scale dispersion of
%    crystallographic orientations. Geology,43

% Prior, D. J., Wheeler, J., Peruzzo, L., Spiess, R., & Storey, C. (2002).
%    Some garnet microstructures: an illustration of the potential of
%    orientation maps and misorientation analysis in microstructural
%    studies. Journal of Structural Geology,24

% Reddy, S. M., Timms, N. E., Pantleon, W., & Trimby, P. (2007).
%    Quantitative characterization of plastic deformation of zircon and
%    geological implications. Contributions to Mineralogy and Petrology,153

%% Clear the workspace and command window. Close all figures.
clear all
close all
clc

% plotting convention (Check this with the RAW ebsd file before executing)
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% Subgrain boundary trace analysis

% set the threshold angles to distinguish between tilt and twist boundaries
thresh_tilt = 15; 
thresh_tilt_twist = 75;

% set the range of misorientation angles to select particular subgrain boundaries
sgb_low = 4; % we have chosen 4, 5, and 7 degrees for different grains 
sgb_high = 10;

% Plot the selected grain and the subgrain boundaries
load('oslo.mat') % colormap by Crameri et al. (2020)
clc
ebsd = ebsd('Magnetite'); % specify the phase of interest
gmag = grains('Magnetite'); % choose the corresponding grains

% set the corresponding crystal symmetry
cs = crystalSymmetry('m-3m',[8.4 8.4 8.4],'mineral','Magnetite');

% plot the mis2mean map
m2mean = ebsd.mis2mean.angle./degree;
plot(ebsd,m2mean,'figsize','medium');
hold on
plot(gmag.boundary,'figsize','medium');
colormap(flipud(oslo))
mtexColorbar('Title','misorientation to mean (\circ)','Location','southoutside')
mtexTitle('Select the grain of interest - ONE mouse click','FontSize',12)
hold off

% manually select the grain of interest
[x_select,y_select] = ginput(1);
selected_grain_ID = findByLocation(gmag,[x_select,y_select]);
gmag = gmag(selected_grain_ID); % assign ID of the selected grain

% select all the low-angle boundaries
sgb_mag = gmag.innerBoundary('Magnetite','Magnetite');
mis_ang = sgb_mag.misorientation.angle./degree;
cond = mis_ang>=sgb_low & mis_ang<sgb_high;
subGB = sgb_mag(cond);

% plot the orientation map along with all the subgrain boundaries
if (isempty(subGB)==0)
    oM = ipfHSVKey(cs.Laue);
    oM.inversePoleFigureDirection = vector3d.X;
    oM_colors = oM.orientation2color(ebsd(gmag).orientations);
    plot(ebsd(gmag),oM_colors); hold on
    plot(gmag.boundary,'linewidth',1.5);
    
    % colorize the subgrain boundaries according to the angle between the 
    % traces of the boundaries and the projection of the 
    % corresponding misorientation axes
    oriGB = ebsd('id',subGB.ebsdId).orientations;
    axS = axis(oriGB(:,1),oriGB(:,2),'antipodal');
    sgb_ang = angle(cross(cross(axS,vector3d.Z),vector3d.Z),...
        subGB.direction)./degree;
    plot(subGB,sgb_ang,'linewidth',2)
    mtexColorMap blue2red
    mtexColorbar('Title','\alpha (\circ)','Location','southoutside')
    mtexTitle('Orientation map (IPF-X) with ','FontSize',12)
end

% select the region of interest
prompt = {'How many polygons do you wish to draw?'};
prompt_title = 'Polygon selection';
lines = 1;
npoly = str2double(inputdlg(prompt,prompt_title,lines));
poly = [];
for i = 1:npoly
    p = selectPolygon;
    poly = [poly;p];
end
save('poly.mat','poly') % save the polygon variable
close all

% derive the misorientation axes corresponding to all the 
% subgrain boundaries enclosed within the polygon
subGB = subGB(inpolygon(subGB,poly));
oriGB = ebsd('id',subGB.ebsdId).orientations;
axS = axis(oriGB(:,1),oriGB(:,2),'antipodal');
sgb_ang = angle(cross(cross(axS,vector3d.Z),vector3d.Z),...
    subGB.direction)./degree; 
ebsdg = ebsd(gmag);
ebsd_poly = ebsdg(inpolygon(ebsdg,poly));

% derive the misorientation axes maximum for the subgrain boundaries
% enclosed within the polygon
if (isempty(subGB)==0)    
    density = calcDensity(axS,'halfwidth',5*degree);
    [~,MisAx] = max(density);    
end
hold off

% marker colors
mcolors = {[50 100 100]/255;...
           [250 150 0]/255;...
           [150 50 50]/255}; % for <100>,<110>, and <111>
sgb_col = [86 166 219]/255; % subgrain bounadry wall and its pole
slip_plane_col = [204 51 102]/255; % slip plane and its pole
disp_axis_col = [255 255 153]/255; % dispersion axis

% marker type
mtype1 = 's'; % for <100>,<110>, and <111>
mtype2 = 'o'; % for other

% display names
disp_name = {'<100>';'<110>';'<111>'};

% define the miller indices
mag_mInd = {Miller(1,0,0,ebsd.CS,'uvw');...
            Miller(1,1,0,ebsd.CS,'uvw');...
            Miller(1,1,1,ebsd.CS,'uvw')};

% calculate the pole figure components correpsonding to the selected grain
vec_axes = {gmag.meanOrientation*symmetrise(mag_mInd{1},'unique','antipodal');...
            gmag.meanOrientation*symmetrise(mag_mInd{2},'unique','antipodal');...
            gmag.meanOrientation*symmetrise(mag_mInd{3},'unique','antipodal')...
            };
        
% separate tilt and tilt/twist boundaries
subGB_tilt = subGB(sgb_ang<thresh_tilt);
subGB_tilt_twist = subGB(sgb_ang>thresh_tilt_twist);
fprintf('Identifying subgrain boundary types. Please wait... \n\n')

if (length(subGB_tilt_twist)>=10) % analyze at least 10 boundary traces
    fprintf('Misorientation axes at high angle to the subgrain boundary traces! \n\n')
    fprintf('Subgrain boundaries could either be tilt OR twist! \n\n')   
    mtexfig = mtexFigure('figsize','huge');
    str_title = sprintf('Twist(?) (n_b = %i)\n\n',...
            length(subGB_tilt_twist));
    sgtitle(str_title,'FontSize',12);
    nextAxis(mtexfig)

    % Plot the axes and their symmetrically equivalent ones
    for i = 1:length(vec_axes)
        plot(vec_axes{i},'antipodal','Marker',mtype1,...
            'MarkerSize',8,'MarkerFaceColor',mcolors{i},...
            'DisplayName',disp_name{i});
        hold on
    end
    % Label the axes
    for i = 1:length(mag_mInd)
        uvw_var{i} = symmetrise(mag_mInd{i},'unique','antipodal');
        for j = 1:length(uvw_var{i})
            uvw_label = num2str(round(uvw_var{i}.uvw));
            plot(vec_axes{i}(j),'label',uvw_label(j,:),'antipodal','Marker','none');
            hold on
        end
    end 

    % ********************************************
    fprintf('Plot assuming TWIST character! \n\n')
    % ********************************************
    
    % Calculate misorientation axes and the average
    [~,mis_Ax] = max(calcDensity(subGB_tilt_twist.misorientation.axis,...
        'antipodal','halfwidth',10*degree));
    oriGB = ebsd('id',subGB_tilt_twist.ebsdId).orientations;
    axS = axis(oriGB(:,1),oriGB(:,2),'antipodal');    
    density = calcDensity(axS,'halfwidth',5*degree);
    [~,MisAx] = max(density);
    
    % Derive the mean trend of the subgrain boundary segments
    sgb_trace = mean(subGB_tilt_twist.direction,'omitnan','antipodal');
    sgb_pole = MisAx; % subgrain boundary normal
    slip_p_n = MisAx; % slip plane normal

    % plot all the components in the stereoplot
    % subgrain boundary wall
    plot(sgb_pole,'antipodal','plane','linestyle','-','linecolor',sgb_col,...
        'linewidth',1.5,'DisplayName',' subgrain boundary wall');
    hold on
    % misorientation axes
    plot(axS,'Marker','o','MarkerSize',3,'MarkerFaceColor','k',...
        'DisplayName',' misorientation axes (LAX)'); 
    % slip plane 
    plot(slip_p_n,'plane','antipodal','linestyle','-','linecolor',...
        slip_plane_col,'linewidth',1.5,'DisplayName',' slip plane');
    % misorientation axes maximum
    plot(MisAx,'antipodal','Marker','s','MarkerSize',8,...
        'MarkerFaceColor','y','MarkerEdgeColor','k',...
        'DisplayName',' LAX maximum');
    % subgrain boundary pole
    plot(sgb_pole,'antipodal','Marker',mtype2,'MarkerSize',10,...
        'MarkerFaceColor',sgb_col,'MarkerEdgeColor','k',...
        'DisplayName',' subgrain boundary pole'); 
    % slip plane normal
    plot(slip_p_n,'antipodal','Marker',mtype2,'MarkerSize',8,...
        'MarkerFaceColor',slip_plane_col,'MarkerEdgeColor','k',...
        'DisplayName',' slip plane normal');
    annotate([xvector,yvector,zvector],'MarkerSize',15)    
    % subgrain boundary trend on stereonet periphery
    plot(sgb_trace,'antipodal','Marker','x','MarkerSize',8,...
        'MarkerFaceColor',sgb_col,'LineWidth',2,'DisplayName',' subgrain boundary trend');
    
    % derive the dispersion axis using PGA and plot it
    [ev,mags] = PGA(ebsd_poly.orientations);    
    plot(ev(1),'antipodal','Marker',mtype2,'MarkerSize',8,...
        'MarkerFaceColor',disp_axis_col,'MarkerEdgeColor','k',...
        'DisplayName',' dipsersion axis')
    
    % add a legend
    lgnd = legend('location','southoutside','NumColumns',2,'FontSize',10);
    lgnd.Title.String = 'Key';
    lgd.Title.FontSize = 10;
    hold off
    
    % Save the plot
    saveas(gcf,'stereoplot_twist.png')
    close all

    % ********************************************
    fprintf('Plot assuming TILT character! \n\n')
    % ********************************************  
    mtexfig = mtexFigure('figsize','huge');
    str_title = sprintf('Tilt(?) (n_b = %i)\n\n',...
            length(subGB_tilt_twist));
    sgtitle(str_title,'FontSize',12);
    nextAxis(mtexfig)
    % Plot the axes and their symmetrically equivalent ones
    for i = 1:length(vec_axes)
        plot(vec_axes{i},'antipodal','Marker',mtype1,...
            'MarkerSize',8,'MarkerFaceColor',mcolors{i},...
            'DisplayName',disp_name{i});
        hold on
    end
    % Label the axes
    for i = 1:length(mag_mInd)
        uvw_var{i} = symmetrise(mag_mInd{i},'unique','antipodal');
        for j = 1:length(uvw_var{i})
            uvw_label = num2str(round(uvw_var{i}.uvw));
            plot(vec_axes{i}(j),'label',uvw_label(j,:),'antipodal','Marker','none');
            hold on
        end
    end 
    
    % Derive the mean trend of the boundary segments
    sgb_trace = mean(subGB_tilt_twist.direction,'omitnan','antipodal');
    sgb_pole = cross(sgb_trace,MisAx); % subgrain boundary normal
    slip_p_n = cross(sgb_pole,MisAx); % slip plane normal

    % plot all the components in the stereoplot
    % subgrain boundary wall
    plot(sgb_pole,'antipodal','plane','linestyle','-','linecolor',sgb_col,...
        'linewidth',1.5,'DisplayName',' subgrain boundary wall');
    hold on
    % misorientation axes
    plot(axS,'Marker','o','MarkerSize',3,'MarkerFaceColor','k',...
        'DisplayName',' misorientation axes (LAX)'); 
    % slip plane 
    plot(slip_p_n,'plane','antipodal','linestyle','-','linecolor',...
        slip_plane_col,'linewidth',1.5,'DisplayName',' slip plane');
    % misorientation axes maximum
    plot(MisAx,'antipodal','Marker','s','MarkerSize',8,...
        'MarkerFaceColor','y','MarkerEdgeColor','k',...
        'DisplayName',' LAX maximum');
    % subgrain boundary pole
    plot(sgb_pole,'antipodal','Marker',mtype2,'MarkerSize',10,...
        'MarkerFaceColor',sgb_col,'MarkerEdgeColor','k',...
        'DisplayName',' subgrain boundary pole'); 
    % slip plane normal
    plot(slip_p_n,'antipodal','Marker',mtype2,'MarkerSize',8,...
        'MarkerFaceColor',slip_plane_col,'MarkerEdgeColor','k',...
        'DisplayName',' slip plane normal');
    annotate([xvector,yvector,zvector],'MarkerSize',15)    
    % subgrain boundary trend on stereonet periphery
    plot(sgb_trace,'antipodal','Marker','x','MarkerSize',8,...
        'MarkerFaceColor',sgb_col,'LineWidth',2,'DisplayName',' subgrain boundary trend');
    
    % derive the dispersion axis using PGA and plot it
    [ev,mags] = PGA(ebsd_poly.orientations);    
    plot(ev(1),'antipodal','Marker',mtype2,'MarkerSize',8,...
        'MarkerFaceColor',disp_axis_col,'MarkerEdgeColor','k',...
        'DisplayName',' dipsersion axis') 

    % add a legend
    lgnd = legend('location','southoutside','NumColumns',2,'FontSize',10);
    lgnd.Title.String = 'Key';
    lgd.Title.FontSize = 10;
    % Save the plot
    saveas(gcf,'stereoplot_tilt.png')
    close all

else
    fprintf('Fewer than 10 subgrain boundary traces at high angle to misorientation axes!\n')
    fprintf('Plots skipped...\n\n')
end

if (length(subGB_tilt)>=10)
    % *********************************************************
    fprintf('Plotting for TILT ONLY subgrain boundaries! \n\n')
    % *********************************************************
    mtexfig = mtexFigure('figsize','huge');
    str_title = sprintf('Tilt (n_b = %i)\n\n',...
            length(subGB_tilt));
    sgtitle(str_title,'FontSize',12);
    nextAxis(mtexfig)
    % Plot the axes and their symmetrically equivalent ones
    for i = 1:length(vec_axes)
        plot(vec_axes{i},'antipodal','Marker',mtype1,...
            'MarkerSize',8,'MarkerFaceColor',mcolors{i},...
            'DisplayName',disp_name{i});
        hold on
    end
    % Label the axes
    for i = 1:length(mag_mInd)
        uvw_var{i} = symmetrise(mag_mInd{i},'unique','antipodal');
        for j = 1:length(uvw_var{i})
            uvw_label = num2str(round(uvw_var{i}.uvw));
            plot(vec_axes{i}(j),'label',uvw_label(j,:),'antipodal','Marker','none');
            hold on
        end
    end

    % Calculate misorientation axes and the average    
    oriGB = ebsd('id',subGB_tilt.ebsdId).orientations;
    axS = axis(oriGB(:,1),oriGB(:,2),'antipodal');    
    density = calcDensity(axS,'halfwidth',5*degree);
    [~,MisAx] = max(density);

    % Derive the mean trend of the TILT boundary segments
    sgb_trace = mean(subGB_tilt.direction,'antipodal','omitnan');
    sgb_pole = cross(sgb_trace,MisAx); % subgrain boundary normal
    slip_p_n = cross(sgb_pole,MisAx); % slip plane normal

    % plot all the components in the stereoplot
    % subgrain boundary wall
    plot(sgb_pole,'antipodal','plane','linestyle','-','linecolor',sgb_col,...
        'linewidth',1.5,'DisplayName',' subgrain boundary wall');
    hold on
    % misorientation axes
    plot(axS,'Marker','o','MarkerSize',3,'MarkerFaceColor','k',...
        'DisplayName',' misorientation axes (LAX)'); 
    % slip plane 
    plot(slip_p_n,'plane','antipodal','linestyle','-','linecolor',...
        slip_plane_col,'linewidth',1.5,'DisplayName',' slip plane');
    % misorientation axes maximum
    plot(MisAx,'antipodal','Marker','s','MarkerSize',8,...
        'MarkerFaceColor','y','MarkerEdgeColor','k',...
        'DisplayName',' LAX maximum');
    % subgrain boundary pole
    plot(sgb_pole,'antipodal','Marker',mtype2,'MarkerSize',10,...
        'MarkerFaceColor',sgb_col,'MarkerEdgeColor','k',...
        'DisplayName',' subgrain boundary pole'); 
    % slip plane normal
    plot(slip_p_n,'antipodal','Marker',mtype2,'MarkerSize',8,...
        'MarkerFaceColor',slip_plane_col,'MarkerEdgeColor','k',...
        'DisplayName',' slip plane normal');
    annotate([xvector,yvector,zvector],'MarkerSize',15)    
    % subgrain boundary trend on stereonet periphery
    plot(sgb_trace,'antipodal','Marker','x','MarkerSize',8,...
        'MarkerFaceColor',sgb_col,'LineWidth',2,'DisplayName',' subgrain boundary trend');
    
    % derive the dispersion axis using PGA and plot it
    [ev,mags] = PGA(ebsd_poly.orientations);    
    plot(ev(1),'antipodal','Marker',mtype2,'MarkerSize',8,...
        'MarkerFaceColor',disp_axis_col,'MarkerEdgeColor','k',...
        'DisplayName',' dipsersion axis')
    
    % add a legend
    lgnd = legend('location','southoutside','NumColumns',2,'FontSize',10);
    lgnd.Title.String = 'Key';
    lgd.Title.FontSize = 10;
    hold off
    
    % Save the plot 
    saveas(gcf,'stereoplot_tilt_only.png')
    close all

else
    fprintf('Fewer than 10 TILT ONLY subgrain boundary traces!\n')
    fprintf('Plots skipped...\n\n')
end    
