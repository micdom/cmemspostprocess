% M. De Dominicis (2020)
% monthly averages from CMEMS data of SST and BT for ORE Supergen

clear all
path(path,'/login/micdom/matlab/m_map')
basedir1=('/scratch/micdom/CMEMS_data_releaseDec2020/Temperature/')
basedir2=('/scratch/micdom/CMEMS_data_releaseDec2020/BottomTemperature/')

root_name1=('metoffice_foam1_amm7_NWS_TEM_dm')
root_name2=('metoffice_foam1_amm7_NWS_BED_dm')
Year=['1993';'1994';'1995';'1996';'1997';'1998';'1999';...
    '2000';'2001';'2002';'2003';'2004';'2005';'2006';'2007';'2008';'2009';'2010';'2011';'2012';'2013';'2014';'2015';'2016';'2017';'2018';'2019'];

YY=length(Year)%-1

% fmask='/scratch/micdom/CMEMS_data/NWS-MFC_004_001_b_mask_bathy.nc'
% tmask=double(ncread(fmask,'mask'));
% tmask(tmask==0)=nan;

%myvar={'BT'};
myvar={'SST','BT'};
var=2
%for var=2%:length(myvar)
dummyfile='/scratch/micdom/CMEMS_data_releaseDec2020/Temperature/1993/01/metoffice_foam1_amm7_NWS_TEM_dm19930101.nc'
longitude=ncread(dummyfile,'longitude');
latitude=ncread(dummyfile,'latitude');
I=find(longitude>180);longitude(I)=longitude(I)-360;
[latitude,longitude]=meshgrid(latitude,longitude);

for year=1:YY
  y1=Year(year,:);
  for Month=1:12
      if Month<10;mm=['0' num2str(Month) ];else;mm=num2str(Month);end

      if Month==1;dmonth=30; elseif Month==2;dmonth=28; elseif Month==3;dmonth=31; elseif Month==4;dmonth=30;
         elseif Month==5;dmonth=31; elseif Month==6;dmonth=30; elseif Month==7;dmonth=31; elseif Month==8;dmonth=31;
         elseif Month==9;dmonth=30; elseif Month==10;dmonth=31; elseif Month==11;dmonth=30; elseif Month==12;dmonth=31;  
      end
      if strcmp(y1,'1996') || strcmp(y1,'2000') || strcmp(y1,'2004')|| strcmp(y1,'2008') || strcmp(y1,'2012') || strcmp(y1,'2016')
         if Month==2; dmonth=29;end
      end
      
      
      for Day=1:dmonth
          if Day<10;dd=['0' num2str(Day) ];else;dd=num2str(Day);end
          
          if var==1
           fname = [basedir1 y1 '/' mm '/' root_name1 y1 mm dd '.nc']  
           T=ncread(fname,'thetao');
           SST(:,:,Day)=squeeze(T(:,:,1));
           VAR=SST;
%            I=isnan(tmask(:,:,1)); VAR(I)=nan;
         elseif var==2
            fname = [basedir2 y1 '/' mm '/' root_name2 y1 mm dd '.nc']  
            BT(:,:,Day)=ncread(fname,'bottomT');
            VAR=BT;
%            VAR=squeeze(BT(:,:,50));
%            I=isnan(tmask(:,:,50)); VAR(I)=nan;
          end
      end %day loop
      VAR_mean_month(:,:,Month)=mean(VAR,3);
      VAR_var_month(:,:,Month)=std(VAR,0,3);
      VAR_median_month(:,:,Month)=median(VAR,3);
      VAR_max_month(:,:,Month)=max(VAR,[],3);
      VAR_min_month(:,:,Month)=max(VAR,[],3);
     
  end %month loop

%WRITE NETCDF OUTPUT FILE

%Open the file
ncid = netcdf.create(['/scratch/micdom/CMEMS_data_releaseDec2020/postproc/' myvar{var} '_' y1 '.nc'],'NC_WRITE')
 
%Define the dimensions
dimidt = netcdf.defDim(ncid,'time',12);
dimidlat = netcdf.defDim(ncid,'y',375);
dimidlon = netcdf.defDim(ncid,'x',297);

%Define IDs for the dimension variables (pressure,time,latitude,...)
%time_ID=netcdf.defVar(ncid,'time','double',[dimidt]);
lon_ID = netcdf.defVar(ncid,'longitude','double',[dimidlon dimidlat]);
lat_ID = netcdf.defVar(ncid,'latitude','double',[dimidlon dimidlat]);

used_varids= cell(0)
used_vnames= cell(0)
var_name=[ myvar{var} '_mean'];used_vnames = [used_vnames, var_name];
var_id=[myvar{var} '_avg_ID'];used_varids = [used_varids, var_id];
var_name=[ myvar{var} '_std'];used_vnames = [used_vnames, var_name];
var_id=[myvar{var} '_std_ID'];used_varids = [used_varids, var_id];
var_name=[ myvar{var} '_median'];used_vnames = [used_vnames, var_name];
var_id=[myvar{var} '_med_ID'];used_varids = [used_varids, var_id];
var_name=[ myvar{var} '_max'];used_vnames = [used_vnames, var_name];
var_id=[myvar{var} '_max_ID'];used_varids = [used_varids, var_id];
var_name=[ myvar{var} '_min'];used_vnames = [used_vnames, var_name];
var_id=[myvar{var} '_min_ID'];used_varids = [used_varids, var_id];
   
for ff=1:5
    eval(['' used_varids{ff} '=netcdf.defVar(ncid,''' used_vnames{ff} ''',''double'',[dimidlon dimidlat dimidt])'])
end

%We are done defining the NetCdf
netcdf.endDef(ncid);


%Then store the dimension variables in
%netcdf.putVar(ncid,time_ID,time);
netcdf.putVar(ncid,lat_ID,latitude);
netcdf.putVar(ncid,lon_ID,longitude);

%Then store my main variables
eval(['netcdf.putVar(ncid,' used_varids{1} ',VAR_mean_month)']);
eval(['netcdf.putVar(ncid,' used_varids{2} ',VAR_var_month)']);
eval(['netcdf.putVar(ncid,' used_varids{3} ',VAR_median_month)']);
eval(['netcdf.putVar(ncid,' used_varids{4} ',VAR_max_month)']);
eval(['netcdf.putVar(ncid,' used_varids{5} ',VAR_min_month)']);

%We're done, close the netcdf
netcdf.close(ncid)
%end
end %year loop


clear all
path(path,'/login/micdom/matlab/cmocean_v1.4/cmocean')
Year=['1993';'1994';'1995';'1996';'1997';'1998';'1999';...
    '2000';'2001';'2002';'2003';'2004';'2005';'2006';'2007';'2008';'2009';'2010';'2011';'2012';'2013';'2014';'2015';'2016';'2017';'2018';'2019'];
YY=length(Year)%-1


for year=1:YY
y1=Year(year,:);
SST_mean(:,:,:,year)=ncread(['/scratch/micdom/CMEMS_data_releaseDec2020/postproc/SST_' y1 '.nc'],'SST_mean');
SST_max(:,:,:,year)=ncread(['/scratch/micdom/CMEMS_data_releaseDec2020/postproc/SST_' y1 '.nc'],'SST_max');
end

figure
for year=1:YY
         year,
         pippo=year;
    for Month=1:12
         Month,

h=subplot_tight(12,YY,pippo,[0.0001])
%pcolor(squeeze(SST_mean(:,:,Month,year))'); shading flat; caxis([-2 20]);colormap(jet);%cmocean('thermal',20);
pcolor(squeeze(SST_max(:,:,Month,year))'); shading flat; caxis([-2 20]);colormap(jet);%cmocean('thermal',20);

axis off
pippo=pippo+YY
   end
 end


for year=1:YY
y1=Year(year,:);

BT_mean(:,:,:,year)=ncread(['/scratch/micdom/CMEMS_data_releaseDec2020/postproc/BT_' y1 '.nc'],'BT_mean');
BT_max(:,:,:,year)=ncread(['/scratch/micdom/CMEMS_data_releaseDec2020/postproc/BT_' y1 '.nc'],'BT_max');

end


figure
for year=1:YY
         year,
         pippo=year;
    for Month=1:12
         Month,

h=subplot_tight(12,YY,pippo,[0.0001])
%pcolor(squeeze(BT_mean(:,:,Month,year))'); shading flat; caxis([-2 20]); colormap(jet); 
pcolor(squeeze(BT_max(:,:,Month,year))'); shading flat; caxis([-2 20]); colormap(jet); 

axis off
pippo=pippo+YY
    end
end
 
