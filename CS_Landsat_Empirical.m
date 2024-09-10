function [Carbon_Stock1,Carbon_Stock_ton,Carbon_Stock_ha,C02_year]=CS_Landsat_Empirical(agb_baseline,ref_baseline,NDVI_baseline,EVI_baseline,L8_dummy,ref_dummy,info_dummy,predicted_year);
% This function needs m_map toolbox to be successfully run
% If no m_map toolbox available please add '%' or delete line 138-194. The
% generated map is still available to be opened in other GIS/RS softwares



NDVIImage = double(NDVI_baseline);
NDVIImage = NDVIImage./10000;
EVIImage = double(EVI_baseline);
EVIImage = EVIImage./10000;
agbImage = double(agb_baseline);
agbVector = agbImage(:);

% Ensure images are the same size
if ~isequal(size(NDVIImage), size(agbImage))
    error('Images must be of the same size.');
else if ~isequal(size(EVIImage), size(agbImage))
    error('Images must be of the same size.');
end
end


% Min-Max Normalization because the indices used here have different value
% range
minNDVI = min(min(NDVIImage));
maxNDVI = max(max(NDVIImage));
normalizedNDVI = (NDVIImage - minNDVI) ./ (maxNDVI - minNDVI);
NDVIVector=normalizedNDVI(:);
minEVI = min(min(EVIImage));
maxEVI = max(max(EVIImage));
normalizedEVI = (EVIImage - minEVI) ./ (maxEVI - minEVI);
EVIVector=normalizedEVI(:);
% Perform linear regression
% Model: AGB = a + b * Index
mdl_NDVI = fitlm(NDVIVector, agbVector);
mdl_EVI = fitlm(EVIVector, agbVector);

% Extract coefficients
slope_NDVI = mdl_NDVI.Coefficients.Estimate(2);
intercept_NDVI = mdl_NDVI.Coefficients.Estimate(1);
slope_EVI = mdl_EVI.Coefficients.Estimate(2);
intercept_EVI = mdl_EVI.Coefficients.Estimate(1);


% Predict AGB using the regression model
predictedAGB_NDVI = intercept_NDVI + (slope_NDVI .* normalizedNDVI);
predictedAGB1_NDVI=predictedAGB_NDVI(:);
predictedAGB_EVI = intercept_EVI + (slope_EVI .* normalizedEVI);
predictedAGB1_EVI=predictedAGB_EVI(:);

% Accuracy Assessment
corrNDVI=corr(predictedAGB1_NDVI,agbVector);
rmseNDVI=rmse(predictedAGB1_NDVI,agbVector);
corrEVI=corr(predictedAGB1_EVI,agbVector);
rmseEVI=rmse(predictedAGB1_EVI,agbVector);

corr_all=cat(2,corrNDVI,corrEVI);
rmse_all=cat(2,rmseNDVI,rmseEVI);

[max_corr,idx_corr]=max(corr_all);
[min_rmse,idx_rmse]=min(rmse_all);

% Predicted Above Ground Biomass L8
blueBand = double(L8_dummy(:,:,1));  % Band 2 (Blue)
redBand = double(L8_dummy(:,:,3));   % Band 4 (Red)
nirBand = double(L8_dummy(:,:,4));   % Band 5 (NIR)

% Constants for EVI calculation
G = 2.5;  % Gain factor
C1 = 6;   % Coefficient for red
C2 = 7.5; % Coefficient for blue
L = 1;    % Canopy background adjustment

% Calculate EVI using the formula
NDVI_dummy1=(nirBand-redBand)./(nirBand+redBand);
EVI_dummy1 = G * ((nirBand - redBand) ./ (nirBand + C1 * redBand - C2 * blueBand + L));


  if idx_corr == 1 && idx_rmse == 1
    disp('Using NDVI for Above Biomass Ground Estimation');
    %minNDVI1 = min(min(NDVI_dummy1));
    %maxNDVI1 = max(max(NDVI_dummy1));
    normalizedNDVI1 = (NDVI_dummy1 - minNDVI) ./ (maxNDVI - minNDVI);
    NDVI_dummy2=normalizedNDVI1;
    predictedAGB= intercept_NDVI + (slope_NDVI .* NDVI_dummy2);
else if idx_corr == 2 && idx_rmse == 2
    disp('Using EVI for Above Biomass Ground Estimation');
    normalizedEVI1 = (EVI_dummy1 - minEVI) ./ (maxEVI - minEVI);
    EVI_dummy2=normalizedEVI1;
    predictedAGB= intercept_EVI + (slope_EVI .* EVI_dummy2);
else disp('Prediction Error')
end
end



%
predictedBGB1 = 0.2 .* predictedAGB;
Total_Biomass1 = predictedAGB + predictedBGB1;
Carbon_Stock1 = Total_Biomass1.*0.5;
Carbon_Stock1(Carbon_Stock1<0)=-9999;
Carbon_Stock1(Carbon_Stock1==-9999)=NaN;
C02_1=Carbon_Stock1.*3.67;

%
pixel_size=(ref_dummy.CellExtentInLatitude*1000)*(ref_dummy.CellExtentInLongitude*1000);
pixel_ha=1/pixel_size;
valid_area=sum((sum(~isnan(Carbon_Stock1))));
Area_size=(pixel_size*valid_area);
Carbon_Stock_ton=sum(sum(Carbon_Stock1,'omitnan'));
Carbon_Stock_ha=(Carbon_Stock_ton/Area_size)/pixel_ha;
C02_year=Carbon_Stock_ton*3.67;


year_used=num2str(predicted_year);
pop_year=(['Carbon Stock Estimation in ',year_used,' (Landsat 8)']);
tcbs= num2str(Carbon_Stock_ton,'%.2f');
tcbs1=(['Total Carbon Stock per year within the region:  ',tcbs,' ton/year']);
cbsh= num2str(Carbon_Stock_ha,'%.2f');
cbsh1=(['Estimated Carbon Stock in hectare per year:  ',cbsh,' ton/ha']);
cosh= num2str(C02_year,'%.2f');
cosh1=(['Total CO2 Equivalent per year within the region:  ',cosh,' tCO2e/year']);
uiwait(msgbox({tcbs1, cbsh1, cosh1},pop_year));


disp('-------------------------');
disp(['Carbon Stock Estimation in ',year_used,' (Landsat 8)']);
disp('-------------------------');
disp(['Total Carbon Stock per year within the region:  ',tcbs,' ton/year']);
disp(['Estimated Carbon Stock in hectare per year:  ',cbsh,' ton/ha']);
disp(['Total CO2 Equivalent per year within the region:  ',cosh,' tCO2e/year']);
disp('-------------------------');



%Map display
x1=ref_dummy.LongitudeLimits(1);
x2=ref_dummy.LongitudeLimits(2);
dx=ref_dummy.CellExtentInLongitude;
X=x1:dx:x2;
y1=ref_dummy.LatitudeLimits(1);
y2=ref_dummy.LatitudeLimits(2);
dy=ref_dummy.CellExtentInLatitude;
Y=y2:-dy:y1;
X=X-dx/2;
X=X(:,[1:end-1]);
Y=Y-dy/2;
Y=Y(:,[1:end-1]);
longi1=X(1,[end-2]);
lati1=Y(2);

[lon,lat]=meshgrid(X,Y);

x12=ref_baseline.LongitudeLimits(1);
x22=ref_baseline.LongitudeLimits(2);
dx1=ref_baseline.CellExtentInLongitude;
X1=x12:dx1:x22;
y12=ref_baseline.LatitudeLimits(1);
y22=ref_baseline.LatitudeLimits(2);
dy12=ref_baseline.CellExtentInLatitude;
Y1=y22:-dy12:y12;
X1=X1(:,[1:end-1]);
Y1=Y1(:,[1:end-1]);


[lon1,lat1]=meshgrid(X1,Y1);

t=tiledlayout(1,2,'TileSpacing','tight');
nexttile
m_proj('mercator','long',[x12 longi1],'lat',[y1 lati1]); 
m_pcolor(lon1,lat1,agb_baseline);shading flat;
hold on
m_gshhs_i('color','k');
m_grid('linestyle','none','tickdir','in','xtick',100:0.05:120,'ytick',0:0.05:50,'fontsize',10); 
caxis([0 160])
cmap = cbrewer2('RdYlGn');
colormap(cmap);
title('Baseline AGB from WHRC (2010)')
nexttile
m_proj('mercator','long',[x12 longi1],'lat',[y1 lati1]); 
m_pcolor(lon,lat,Carbon_Stock1);shading flat;
hold on
m_gshhs_i('color','k');
m_grid('linestyle','none','tickdir','in','xtick',100:0.05:120,'ytick',0:0.05:50,'fontsize',10); 
caxis([0 160])
cmap = cbrewer2('RdYlGn');
colormap(cmap);
tit=(['Predicted Carbon Stock from Landsat 8 in ',year_used]);
title(tit)

title(t,'Generated Map For Carbon Stock Estimation','color','k','FontWeight','bold');
h=colorbar;
h.Label.String='ton/year';

% Saving report
fname=(['Result_Report for Carbon Estimation in ',year_used,'.txt']);
file=fopen(fname,'wt');
fprintf(file,['Carbon Stock Estimation in ',year_used,' (Landsat 8)\n']);
fprintf(file,'-------------------------\n');
fprintf(file,['Total Carbon Stock per year within the region:  ',tcbs,' ton/year\n']);
fprintf(file,['Estimated Carbon Stock in hectare per year:  ',cbsh,' ton/ha\n']);
fprintf(file,['Total CO2 Equivalent per year within the region:  ',cosh,' tCO2e/year\n']);
fprintf(file,'-------------------------\n');

%Saving map
fname_map1=(['Result_CarbonStock(Landsat8)_',year_used,'.tif']);
geotiffwrite(fname_map1, Carbon_Stock1, ref_dummy, 'GeoKeyDirectoryTag',...
info_dummy.GeoTIFFTags.GeoKeyDirectoryTag);



end
