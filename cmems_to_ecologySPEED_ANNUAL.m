% M. De Dominicis (2020)
% monthly averages from CMEMS data of horizontal speed for ORE Supergen

clear all, close all
path(path,'/login/micdom/matlab/m_map')
basedir=('/scratch/micdom/CMEMS_data_releaseDec2020/Velocity/')

root_name=('metoffice_foam1_amm7_NWS_CUR_dm')
Year=['1993';'1994';'1995';'1996';'1997';'1998';'1999';...
    '2000';'2001';'2002';'2003';'2004';'2005';'2006';'2007';'2008';'2009';'2010';'2011';'2012';'2013';'2014';'2015';'2016';'2017';'2018';'2019'];

YY=length(Year)%-1

%fmask='/scratch/micdom/CMEMS_data/NWS-MFC_004_001_b_mask_bathy.nc'
%bathy=ncread(fmask,'deptho');bathy(bathy==-999)=nan;
%tmask=double(ncread(fmask,'mask'));%this mask has 33 layers instead then 24....
%tmask(tmask==0)=nan;

myvar={'SPEED';'W'};
var=2
%for var=2%:length(myvar)
dummyfile='/scratch/micdom/CMEMS_data_releaseDec2020/Velocity/1998/01/metoffice_foam1_amm7_NWS_CUR_dm19980101.nc'
longitude=ncread(dummyfile,'longitude');
latitude=ncread(dummyfile,'latitude');
depth=ncread(dummyfile,'depth');
I=find(longitude>180);longitude(I)=longitude(I)-360;
[latitude,longitude]=meshgrid(latitude,longitude);
uo=ncread(dummyfile,'uo');

%adjust the last thickness
for lev=1:23
thickness(lev)=depth(lev+1)-depth(lev)
thick_mtx(:,:,lev)=ones(297,375)*thickness(lev);                 
end

bathymod=nan(297,375);
for i=1:297
 for j=1:375
 last=find(isnan(uo(i,j,:)),1,'first');                     
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
          
          
           fname = [basedir y1 '/' mm '/' root_name y1 mm dd '.nc']  
           uo=ncread(fname,'uo');
           vo=ncread(fname,'vo');
           
           if var==1
           speed=sqrt(uo.^2+vo.^2);
           
             for lev=1:23
             thickness=depth(lev+1)-depth(lev)
             speeda(:,:,lev)=speed(:,:,lev)*thickness;
             end
           speeda_mean(:,:,Day)=(nansum(speeda,3)./bathymod);
           VAR=speeda_mean;
           
           elseif var==2
             w=zeros(297,375,24);
            
             for i=1:297-1
              for j=1:375-1
                  for z=23:-1:1
              
                     dudx(i,j,z)=(uo(i+1,j,z)-uo(i,j,z))/(m_lldist([longitude(i+1,j) longitude(i,j)],[latitude(i+1,j) latitude(i,j)])*1000);
                     dvdy(i,j,z)=(vo(i,j+1,z)-vo(i,j,z))/(m_lldist([longitude(i,j+1) longitude(i,j)],[latitude(i,j+1) latitude(i,j)])*1000);
                     dwdz(i,j,z)=(dudx(i,j,z)+dvdy(i,j,z));
                     thickness=depth(z+1)-depth(z);
                     if isnan(dwdz(i,j,z))==0
                     w(i,j,z)=-dwdz(i,j,z)*thickness+w(i,j,z+1);
                     end
                         
                  end
              end
             end
             
             for lev=1:23
             thickness=depth(lev+1)-depth(lev)
             wa(:,:,lev)=w(:,:,lev)*thickness;
             end
             wa_mean(:,:,Day)=(nansum(wa,3)./bathymod);
             VAR=wa_mean;%w.*mask;
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
W_mean(:,:,:,year)=ncread(['/scratch/micdom/CMEMS_data_releaseDec2020/postproc/W_' y1 '.nc'],'W_mean');
W_max(:,:,:,year)=ncread(['/scratch/micdom/CMEMS_data_releaseDec2020/postproc/W_' y1 '.nc'],'W_max');
end

figure
for year=1:YY
         year,
         pippo=year;
    for Month=1:12
         Month,

h=subplot_tight(12,YY,pippo,[0.0001])
%pcolor(squeeze(W_mean(:,:,Month,year))'); shading flat; caxis([-0.0001 0.0001]);colormap(jet);%cmocean('thermal',20);
pcolor(squeeze(W_max(:,:,Month,year))'); shading flat; caxis([-0.0001 0.0001]);colormap(jet);%cmocean('thermal',20);

axis off
pippo=pippo+YY
   end
 end

