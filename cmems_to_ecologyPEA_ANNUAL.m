% M. De Dominicis (2020)
% monthly averages from CMEMS data of PEA for ORE Supergen

clear all, close all
path(path,'/login/micdom/matlab/m_map')
addpath '/vkamino/work/micdom/EcoWatt2050/scripts_ecowatt'
basedir1=('/scratch/micdom/CMEMS_data_releaseDec2020/Temperature/')
basedir2=('/scratch/micdom/CMEMS_data_releaseDec2020/Salinity/')

root_name1=('metoffice_foam1_amm7_NWS_TEM_dm')
root_name2=('metoffice_foam1_amm7_NWS_SAL_dm')
Year=['1993';'1994';'1995';'1996';'1997';'1998';'1999';...
    '2000';'2001';'2002';'2003';'2004';'2005';'2006';'2007';'2008';'2009';'2010';'2011';'2012';'2013';'2014';'2015';'2016';'2017';'2018';'2019'];

YY=length(Year)%-1

%fmask='/scratch/micdom/CMEMS_data_releaseDec2020/NWS-MFC_004_001_b_mask_bathy.nc'
%bathy=ncread(fmask,'deptho');bathy(bathy==-999)=nan;
%tmask=double(ncread(fmask,'mask'));%this mask has 33 layers instead then 24....
%tmask(tmask==0)=nan;

myvar={'PEA'};
var=1
%for var=2%:length(myvar)
dummyfile='/scratch/micdom/CMEMS_data_releaseDec2020/Temperature/1998/01/metoffice_foam1_amm7_NWS_TEM_dm19980101.nc'
longitude=ncread(dummyfile,'longitude');
latitude=ncread(dummyfile,'latitude');
depth=ncread(dummyfile,'depth');
I=find(longitude>180);longitude(I)=longitude(I)-360;
[latitude,longitude]=meshgrid(latitude,longitude);
temp=ncread(dummyfile,'thetao');

%adjust the last thickness
for lev=1:23
thickness(lev)=depth(lev+1)-depth(lev)
thick_mtx(:,:,lev)=ones(297,375)*thickness(lev); 
H(:,:,lev)=sum(thick_mtx(:,:,1:lev),3);               
end

bathymod=nan(297,375);
for i=1:297
 for j=1:375
 last=find(isnan(temp(i,j,:)),1,'first');                     
 if length(last)>=1 & last<24
 bathymod(i,j)=depth(last);
 else
 bathymod(i,j)=depth(24);
 end

 end
end



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
          
          
           fname1 = [basedir1 y1 '/' mm '/' root_name1 y1 mm dd '.nc']  
           fname2 = [basedir2 y1 '/' mm '/' root_name2 y1 mm dd '.nc']
           temp=ncread(fname1,'thetao');
           salt=ncread(fname2,'so');
           salt(salt<0)=0;
           rho = density_nopress(salt,temp);
           
             for lev=1:23
             thickness=depth(lev+1)-depth(lev)
             rhoa(:,:,lev)=rho(:,:,lev)*thickness;
             end
             rhoa_mean=(nansum(rhoa,3)./bathymod);
           
             PEA_level=zeros(297,375,23);
             for lev=1:23
             thickness=depth(lev+1)-depth(lev)  
             %PEA_level(:,:,i)=(((rho(:,:,i)-rhobar).*g).*H(:,:,i)).*dz_3d_weight(:,:,i);
             PEA_level(:,:,lev)=(((rho(:,:,lev)-rhoa_mean).*9.81).*H(:,:,lev)).*thickness;
             end
             PEA(:,:,Day)=nansum(PEA_level,3)./bathymod;
             PEA(PEA<0)=0;
             VAR=PEA;
           
           
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
PEA_mean(:,:,:,year)=ncread(['/scratch/micdom/CMEMS_data_releaseDec2020/postproc/PEA_' y1 '.nc'],'PEA_mean');
PEA_max(:,:,:,year)=ncread(['/scratch/micdom/CMEMS_data_releaseDec2020/postproc/PEA_' y1 '.nc'],'PEA_max');
end

figure
for year=1:YY
         year,
         pippo=year;
    for Month=1:12
         Month,

h=subplot_tight(12,YY,pippo,[0.0001])
%pcolor(log10(squeeze(PEA_mean(:,:,Month,year))')); shading flat; colormap(jet); %caxis([0 500]);%cmocean('thermal',20);
pcolor(log10(squeeze(PEA_max(:,:,Month,year))')); shading flat; colormap(jet); %caxis([0 500]);%cmocean('thermal',20);
%Contours=[1 10 100 1000 3000];colorbar('YTick',log10(Contours),'YTickLabel',Contours);
caxis(log10([Contours(1) Contours(length(Contours))]));
axis off
pippo=pippo+YY
   end
 end

