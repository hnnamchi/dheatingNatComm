;*************************************************
;*************************************************
begin
;==================================================
latS   = -15.
latN   =  15. 
lonL   = -60.
lonR   =  20.

t1 = 1
t2 = 11	;SST from month 1 to 11
;
p1 = 0
p2 = 10	;Q/SSH from month 0 to 10



nmos = 12
xlabel = ispan(0,nmos-1,1)
print(xlabel)
;
xlabel2 = ispan(0,nmos-2,1)
print(xlabel2)
  ;
km = 1000.

mm = 5		;index counts from zero, June is 5th month

p95 = 0.05

mxlag = 7
lag = 1	;first-order autocorrelation
;
days = (/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./)
days = days*24*60*60	;day,hours,minutes in each month
print(days)
;
;
const = 86400	;1 kg/m2/s = 86400 mm/day	;heat constants


xy_colors_dheating        = (/"red","greenyellow","deepskyblue","black","darksalmon","black"/)
xy_colors_thermocline        = (/"blue","yellow3","red","springgreen4","springgreen","black"/)
;=========================================================================================
;=========================================================================================
;=========================================================================================
;=========================================================================================
;=========================================================================================
;=========================================================================================
yrStrt = 1982
yrLast = 2015
;========================================================================================
;========================================================================================
;
f0      = addfile("/home/hyacinth/work/datasets/cfsr/itcz_wind_stress_cfsr_6hrs_4per_day_1979_2019.nc", "r")
  TIME   = f0->TIME   
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  i_cfsr    = dble2flt(f0->ITCZ(iYYYY,:))

  i_cfsr = lonFlip(i_cfsr)
  printVarSummary(i_cfsr) 
  printMinMax(i_cfsr,True)

      
  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set cfsre time dimension
;-------------------------------
f0      = addfile("/home/hyacinth/work/datasets/era5/itcz_era5_wind_stress_1979_2018.nc", "r")
  TIME   = f0->TIME    
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  i_era5    = dble2flt(f0->ITCZ(iYYYY,:))

  i_era5!1 = "LON"
  i_era5 = lonFlip(i_era5)
  printVarSummary(i_era5) 
  printMinMax(i_era5,True)

       
  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set cfsre time dimension
;-------------------------------
f0      = addfile("/home/hyacinth/work/datasets/erai/itcz_erai_turbulent_fluxes_wind_stress.nc", "r")
  TIME   = f0->INITIAL_TIME0_HOURS   
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  i_erai    = dble2flt(f0->ITCZ(iYYYY,:))
  
  i_erai!0 = "TIME"
  i_erai!1 = "LON"
  i_erai = lonFlip(i_erai)
  printVarSummary(i_erai) 
  printMinMax(i_erai,True)

       
  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set cfsre time dimension
  
;-------------------------------
f0      = addfile("/home/hyacinth/work/datasets/jra55/itcz_jra55_fcst_phy2m125.1958_2019.nc", "r")
  TIME   = f0->TIME    
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  i_jra55    = dble2flt(f0->ITCZ(iYYYY,:))

  i_jra55 = lonFlip(i_jra55)
  printVarSummary(i_jra55) 
  printMinMax(i_jra55,True)

       
  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set cfsre time dimension
;-------------------------------
f0      = addfile("/home/hyacinth/work/datasets/merra/itcz_merra2_MERRA2_100_400.tavgM_2d_ocn_Nx.198001_201906.nc4.nc", "r")
  TIME   = f0->TIME    
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  i_merra    = dble2flt(f0->ITCZ(iYYYY,:))

 ; i_merra = lonFlip(i_merra)
  printVarSummary(i_merra) 
  printMinMax(i_merra,True)

       
  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set cfsre time dimension
 
;---------------------------
itcz = new((/5,dimsizes(i_era5&TIME)/),float,i_era5@_FillValue)
printVarSummary(itcz) 
itcz(4,:) = i_cfsr(TIME|:,{LON|-28})
itcz(3,:) = i_era5(TIME|:,{LON|-28})
itcz(2,:) = i_erai(TIME|:,{LON|-28})
itcz(1,:) = i_jra55(TIME|:,{LON|-28})
itcz(0,:) = i_merra(TIME|:,{LON|-28})

itcz!0 = "lev"
itcz!1 = "time"
itcz&time = i_era5&TIME
itcz&lev = fspan(0,4,5)
print(itcz&time)
print(itcz&lev)
printVarSummary(itcz)
;print(itcz)
i_itcz = itcz(time|:,lev|:)
printVarSummary(i_itcz)
;---------------------
;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
 i_itcz_d = new(dimsizes(i_itcz),typeof(i_itcz),i_itcz@_FillValue)
 x   = i_itcz&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work = i_itcz(tt::nmos,:)
    work = dtrend_msg_n(x1,work,False,False,0)
    i_itcz_d(tt::nmos,:) = work
 end do

 printVarSummary(i_itcz_d)
 i_itcz = i_itcz_d

;===
delete([/tt,x,x1,i_itcz_d,work/])
;=========================================================================================
;=========================================================================================
;--------------------------
;
nyears = dimsizes(i_itcz&time)/nmos
print(nyears)
i_itcz_s     = onedtond(ndtooned(i_itcz),(/nyears,nmos,5/))
i_itcz_s!0 = "time"
i_itcz_s!1 = "month"
i_itcz_s!2 = "lev"
printVarSummary(i_itcz_s)
delete(i_itcz)
;---------------------
printMinMax(i_itcz_s,True)
;i_itcz_s = i_itcz_s/conform(i_itcz,days,1)
printMinMax(i_itcz_s,True)
;
;i_itcz_s = dtrend_n(i_itcz_s,False,0)
;
std_i_itcz_s = dim_stddev_n_Wrap(i_itcz_s,0)
printVarSummary(std_i_itcz_s)
;
mean_i_itcz_s = dim_avg_n_Wrap(i_itcz_s,0)
printVarSummary(mean_i_itcz_s)
;print(xxxxxxxxxxxxx)

;---------------------------
itcz_2sigma = new((/3,nmos/),float,i_itcz_s@_FillValue)
printVarSummary(i_itcz_s) 
itcz_2sigma(0,:) = dim_avg_n_Wrap(mean_i_itcz_s,1)-(2.*dim_avg_n_Wrap(std_i_itcz_s,1))
itcz_2sigma(1,:) = dim_avg_n_Wrap(mean_i_itcz_s,1)
itcz_2sigma(2,:) = dim_avg_n_Wrap(mean_i_itcz_s,1)+(2.*dim_avg_n_Wrap(std_i_itcz_s,1))
itcz_2sigma!0 = "lev"
itcz_2sigma!1 = "time"
printVarSummary(itcz_2sigma)
;print(xxxxxxxxxxxxxxxxxxxxxx)





;---------------------------for the basin
itczb = new((/5,dimsizes(i_era5&TIME)/),float,i_era5@_FillValue)
printVarSummary(itczb) 
itczb(4,:) = dim_avg_n_Wrap(i_cfsr(TIME|:,{LON|-35:-18}),1) 
itczb(3,:) = dim_avg_n_Wrap(i_era5(TIME|:,{LON|-35:-18}),1)
itczb(2,:) = dim_avg_n_Wrap(i_erai(TIME|:,{LON|-35:-18}),1)
itczb(1,:) = dim_avg_n_Wrap(i_jra55(TIME|:,{LON|-35:-18}),1)
itczb(0,:) = dim_avg_n_Wrap(i_merra(TIME|:,{LON|-35:-18}),1)

itczb!0 = "lev"
itczb!1 = "time"
itczb&time = i_era5&TIME
itczb&lev = fspan(0,4,5)
print(itczb&time)
print(itczb&lev)
printVarSummary(itczb)
;print(itczb)
i_itczb = itczb(time|:,lev|:)
printVarSummary(i_itczb)
;---------------------
delete(i_cfsr)
delete(i_era5)
delete(i_erai)
delete(i_jra55)
delete(i_merra)
delete(itczb)
;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
 i_itczb_d = new(dimsizes(i_itczb),typeof(i_itczb),i_itczb@_FillValue)
 x   = i_itczb&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work = i_itczb(tt::nmos,:)
    work = dtrend_msg_n(x1,work,False,False,0)
    i_itczb_d(tt::nmos,:) = work
 end do

 printVarSummary(i_itczb_d)
 i_itczb = i_itczb_d

;===
delete([/tt,x,x1,i_itczb_d,work/])
;=========================================================================================
;=========================================================================================
;--------------------------
;
nyears = dimsizes(i_itczb&time)/nmos
print(nyears)
i_itczb_s     = onedtond(ndtooned(i_itczb),(/nyears,nmos,5/))
i_itczb_s!0 = "time"
i_itczb_s!1 = "month"
i_itczb_s!2 = "lev"
printVarSummary(i_itczb_s)
delete(i_itczb)
;---------------------
printMinMax(i_itczb_s,True)
;i_itczb_s = i_itczb_s/conform(i_itczb,days,1)
printMinMax(i_itczb_s,True)
;
;i_itczb_s = dtrend_n(i_itczb_s,False,0)
;
std_i_itczb_s = dim_stddev_n_Wrap(i_itczb_s,0)
printVarSummary(std_i_itczb_s)
;
mean_i_itczb_s = dim_avg_n_Wrap(i_itczb_s,0)
printVarSummary(mean_i_itczb_s)
;print(xxxxxxxxxxxxx)

;---------------------------
itczb_2sigma = new((/3,nmos/),float,i_itczb_s@_FillValue)
printVarSummary(i_itczb_s) 
itczb_2sigma(0,:) = dim_avg_n_Wrap(mean_i_itczb_s,1)-(2.*dim_avg_n_Wrap(std_i_itczb_s,1))
itczb_2sigma(1,:) = dim_avg_n_Wrap(mean_i_itczb_s,1)
itczb_2sigma(2,:) = dim_avg_n_Wrap(mean_i_itczb_s,1)+(2.*dim_avg_n_Wrap(std_i_itczb_s,1))
itczb_2sigma!0 = "lev"
itczb_2sigma!1 = "time"
printVarSummary(itczb_2sigma)
;print(xxxxxxxxxxxxxxxxxxxxxx)

;=================================================
;=================================================
;======================================================================================
;======================================================================================


;=================================================
;=================================================
;======================================================================================
;======================================================================================
f0      = addfile("/home/hyacinth/work/geomar/sstP/R2/a_sst__1982_2015.nc", "r")
  TIME   = f0->time    
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  a_sst    = f0->a_sst(iYYYY)

  printVarSummary(a_sst) 
  printMinMax(a_sst,True)

       
  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set cfsre time dimension
  


;
nyears = dimsizes(a_sst&time)/nmos
print(nyears)
;
a_sst_s     = onedtond(ndtooned(a_sst),(/nyears,nmos/))
a_sst_s!0 = "time"
a_sst_s!1 = "month"
printVarSummary(a_sst_s)
printMinMax(a_sst_s,True)
;
a_sst_s = dtrend_n(a_sst_s,False,0)
;
std_sst_s = dim_stddev_n_Wrap(a_sst_s,0)
print(std_sst_s)

;
temp = dim_standardize_n(a_sst_s,0,0)      ; normalize
copy_VarMeta(a_sst_s(:,:),temp) 
printVarSummary(temp)

printVarSummary(temp)
atemp = temp(:,5)
atemp!0 = "time"
print(atemp)
delete(temp)
;=========================================================================================
;=========================================================================================
;=========================================================================================
;=========================================================================================
;===============================================================================================================
;================================================================
; ===============================================================
; Open the file: Read only the user specified period
; ==============================================================
f      = addfile ("/home/hyacinth/work/datasets/rain/GPCP_precip.mon.mean.nc", "r")
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

rain    = short2flt(f->precip(iYYYY,:,:))
rain = lonFlip(rain)
printVarSummary(rain)       ; variable overview

printMinMax(rain,True)
;print(xxxxxxxx)

rain!1  =  "lat"
rain!2  =  "lon"
delete([/YYYY,iYYYY,f,TIME/])

;===============Mask ou land points==============
   a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
   lsdata = a->LSMASK
   lsm  = landsea_mask(lsdata,rain&lat,rain&lon)

    rain = mask(rain,lsm.eq.1,False)


  delete(a)
  delete(lsm)
  delete(lsdata)
;==============================================================
;===========================
arain = wgt_areaave_Wrap(rain(time |:, {lon|-35:-18}, {lat | -3:3}),1.0, 1.0, 0)
printVarSummary(arain)
;
nyears = dimsizes(arain&time)/nmos
print(nyears)
;
arain_s     = onedtond(ndtooned(arain),(/nyears,nmos/))
arain_s!0 = "time"
arain_s!1 = "month"
printVarSummary(arain_s)
printMinMax(arain_s,True)
;
arain_s = dtrend_n(arain_s,False,0)
printMinMax(arain_s,True)
;print(xxxxxxxxxx)
;
std_sst_s = dim_stddev_n_Wrap(arain_s,0)
print(std_sst_s)

;
;=========================================================================================
;=========================================================================================
;
;arain_s = dim_standardize_n(arain_s,0,0)      ; normalize 
;
d_arain = regCoef(atemp(time|:),arain_s(month|:,time|:))
;
d_arain!0 = "month"
print(d_arain)
;print(xxxxxxxxxxxxxxxxxx)
;
d_arain = d_arain
print(d_arain)
;==============================================================
;============================================================================================
;==============================================================
;==============================================================

;--------------test for significance

 lag_ts = esacr(atemp(time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(arain_s(month|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

  tval = onedtond(d_arain@tval,dimsizes(d_arain))

       df = floattointeger(dimsizes(arain_s(:,0))*(1. -(lag_ts(lag)*lag_ts_map(:,lag)))\
        /(1. + (lag_ts_map(:,lag)*lag_ts(lag)))) 

printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)

 prob@long_name = "probability"
   printVarSummary(prob)
;

d_araind = d_arain
d_araind = where(prob.le.p95,d_araind,d_araind@_FillValue) ; Significant TEMP
printVarSummary(d_araind)
printMinMax(d_araind,True)
print(d_araind)

;delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(tval)
delete(df)
delete(prob)
delete(b)
delete(arain)
delete(arain_s)
;==========================================================================================
;==========================================================================================
;=============================================================================================================
 rain_d = new(dimsizes(rain),typeof(rain),rain@_FillValue)
 x   = rain&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =rain(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,False,False,0)
    rain_d(tt::nmos,:,:) = work
 end do

 printVarSummary(rain_d)
 rain = rain_d

;===
delete([/tt,x,x1,rain_d,work/])
printMinMax(rain,True)
;print(xxxxxxxxx)
;---------------------------------------------------------------------------
;
;---------------------------------------------------------------------------
;---------------------------------------------------------------------------
nyears = dimsizes(rain&time)/nmos
print(nyears)
rain_s     = onedtond(ndtooned(rain),(/nyears,nmos,dimsizes(rain&lat),dimsizes(rain&lon)/))
copy_VarCoords(rain(0,:,:),rain_s(0,0,:,:))                             ; copy coordinate arrays from slp to newslp
rain_s!0 = "time"
rain_s!1 = "month"
;
printVarSummary(rain_s)
;
;
;
xx = rain_s(time|:,month|:,lat|:,{lon|-28})
yy = rain_s		;keep for the basin average
printVarSummary(xx)
delete(rain_s)
rain_s = xx
delete(xx)
printVarSummary(rain_s)
;
;
d_GPCP = regCoef(atemp,rain_s(month|:,lat|:,time|:))
copy_VarCoords(rain_s(0,:,:),d_GPCP)
printVarSummary(d_GPCP)
;
a_rain_s = dim_avg_n_Wrap(rain_s,0)
printVarSummary(a_rain_s)
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;--------------test for significance
 d_GPCPd = d_GPCP
 
 df = d_GPCP

 lag_ts = esacr(atemp,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(rain_s(month|:,lat|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(d_GPCP@tval,dimsizes(d_GPCP))
     
    ;do ne = 0,neof-1,1
           df(:,:) = floattointeger(dimsizes(rain_s(:,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   printMinMax(df,True)
  
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)


   
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printVarSummary(prob)

  printMinMax(prob,True)
   
;d_GPCPd = where(abs(prob).lt.p95,d_GPCPd,d_GPCPd@_FillValue) 
d_GPCPd = where(abs(prob).lt.p95,prob,d_GPCPd@_FillValue) ; mask out stat sig.
printVarSummary(d_GPCPd)
printMinMax(d_GPCPd,True)

delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(rain_s)
delete(rain)
delete(nyears)
;print(xxxxxxxxxxxxxx)
;===============================================================================
;===============================================================================
rainb_s = dim_avg_n_Wrap(yy(time|:,month|:,lat|:,{lon|-35:-18}),3)	;now based on the basin average
printVarSummary(rainb_s)
;
;
d_GPCPb = regCoef(atemp,rainb_s(month|:,lat|:,time|:))
copy_VarCoords(rainb_s(0,:,:),d_GPCPb)
printVarSummary(d_GPCPb)
;
a_rainb_s = dim_avg_n_Wrap(rainb_s,0)
printVarSummary(a_rainb_s)
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;--------------test for significance
 d_GPCPbd = d_GPCPb
 
 df = d_GPCPb

 lag_ts = esacr(atemp,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(rainb_s(month|:,lat|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(d_GPCPb@tval,dimsizes(d_GPCPb))
     
    ;do ne = 0,neof-1,1
           df(:,:) = floattointeger(dimsizes(rainb_s(:,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   printMinMax(df,True)
  
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)


   
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printVarSummary(prob)

  printMinMax(prob,True)
   
;d_GPCPbd = where(abs(prob).lt.p95,d_GPCPbd,d_GPCPbd@_FillValue) ; 
d_GPCPbd = where(abs(prob).lt.p95,prob,d_GPCPbd@_FillValue) ; mask out stat. sig.
printVarSummary(d_GPCPbd)
printMinMax(d_GPCPbd,True)

delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(rainb_s)

print("==============================================================")
print("==============================================================")
;===============================================================================================================
;======================================
  wks = gsn_open_wks("eps","FIG_S5")         ; send graphics to PNG file
  plot = new(3,graphic)                ; create graphic array ; only needed if paneling
  
  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res@cnMinLevelValF       = -2.0            ; set min contour level
  res@cnMaxLevelValF       =  2.0             ; set max contour level
  res@cnLevelSpacingF      =  0.2            ; set contour spacing   
  res@cnFillOn             = True         ; turn on color fill
  res@cnFillPalette        = "cmp_b2r"  
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@cnInfoLabelOn        = False
  res@lbLabelBarOn         = True       ; turn off individual lb's 
  res@lbOrientation        = "vertical"
  res@lbTitleFontHeightF  = 0.017 
  res@lbLabelFontHeightF  = 0.017
  res@lbLeftMarginF       = 0.10
  res@lbRightMarginF      = 0.55
  res@lbLabelFont   	  = "Helvetica-Bold"     ; label font
  res@lbTitleString       = "Precipitation (mm day~S~-1~N~)"
  res@lbTitleAngleF       = 270.
  res@lbTitlePosition      = "Right"
  res@lbTitleDirection     = "Across"
  
  
  res@lbLabelStride        = 10 
  res@tmYLMode   = "Explicit"
  res@tmYLValues = ispan(-5,20,5)
  res@tmYLLabels = (/"-5", "0","5","10","15","20"/)
  res@tmYLMinorValues  = ispan(-5,20,1) 
  res@tmYLLabelFontHeightF = 0.017

  res@tmXBMode = "Explicit"
  res@tmXBValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
  res@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
  res@tiXAxisString = "Time lag (month)"                  ; yaxis
  res@tmXBLabelFontHeightF = 0.017


  res@tiXAxisFontHeightF = 0.020
  res@tiYAxisFontHeightF = 0.020
  res@vpHeightF 	  = 0.40        ; Changes the aspect ratio
  res@vpWidthF 		  = 0.80
  res@gsnCenterStringFontHeightF   = 0.030
  

  

  res@gsnCenterString  = ""   
  res@gsnLeftString  = ""  
  res@gsnRightString  = ""  
    
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
  res2                      =   True
  res2@gsnDraw              = False        ; don't draw yet
  res2@gsnFrame             = False        ; don't advance frame yet
  res2@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  res2@cnLevels		   = (/0./)   ; set levels
  res2@cnMonoLineThickness  = True
  res2@cnLineThicknessF     = 4.0
  res2@cnMonoLineColor	    = True
  res2@cnLineColor	    = "grey25"
  res2@gsnAddCyclic         = False   
  res2@cnLinesOn            = True       ; True is default
  res2@cnLineLabelsOn       = False       ; True is default
  res2@cnInfoLabelOn        = False
  res2@gsnRightString   = ""
  res2@gsnLeftString   = ""
  res2@cnLineLabelDensityF =  2.0
  res2@cnLabelMasking	   = True
  res2@gsnContourZeroLineThicknessF = 2.
  res2@cnLineDashPattern	=	16
  
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
  res3                      =   True
  res3@gsnDraw              = False        ; don't draw yet
  res3@gsnFrame             = False        ; don't advance frame yet
  res3@cnMinLevelValF       = -20.            ; set min contour level
  res3@cnMaxLevelValF       =  20.             ; set max contour level
  res3@cnLevelSpacingF      =  2.            ; set contour spacing   
  res3@cnMonoLineThickness  = True
  res3@cnLineThicknessF     = 1.5
  res3@cnMonoLineColor	    = True
  res3@cnLineColor	    = "white"
  res3@gsnAddCyclic         = False   
  res3@cnLinesOn            = True       ; True is default
  res3@cnLineLabelsOn       = True       ; True is default
  res3@cnInfoLabelOn        = False
  res3@gsnRightString   = ""
  res3@gsnLeftString   = ""
  res3@cnLineLabelDensityF =  1.0
  res3@cnLabelMasking	   = True
  res3@gsnContourZeroLineThicknessF = 0.
  res3@gsnContourNegLineDashPattern = 16
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
  ;-----------------------------------------------------------------------------
  res1                      =   True
  res1@gsnDraw              = False        ; don't draw yet
  res1@gsnFrame             = False        ; don't advance frame yet
  res1@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res1@cnMinLevelValF       = res@cnMinLevelValF            ; set min contour level
  res1@cnMaxLevelValF       = res@cnMaxLevelValF            ; set max contour level
  res1@cnLevelSpacingF      = res@cnLevelSpacingF           ; set contour spacing 
  res1@cnLinesOn            = False       ; True is default
  res1@cnLineLabelsOn       = False        ; True is default
  res1@cnInfoLabelOn        = False
  res1@gsnRightString       = ""
  res1@gsnLeftString        = ""
  res1@cnFillScaleF         =  1.5         ; add extra density
  res1@cnFillDotSizeF       = 0.001           ; increase dot size of contours
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
  opt = True
  opt@gsnShadeFillType = "pattern"        ; color is the default
  opt@gsnShadeLow = 17                     ; use fill pattern #4
;-----------------------------------------------------------------------------
;=============================================================================
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
   res@gsnCenterString  = ""  
   res@tiYAxisString = "Latitude (~S~o~N~)"                  ; yaxis
   plotA0 = gsn_csm_contour(wks,d_GPCP({lat|-7:17},month|:),res)            ; default hovmueller
   plotB0 = gsn_csm_contour(wks,d_GPCPd({lat|-7:17},month|:),res1)            ; default hovmueller
   ;plotB0 = gsn_contour_shade(plotB0,5.,p95,opt)
   plotB0 = gsn_contour_shade(plotB0,max(d_GPCP),min(d_GPCP),opt)
   plotC0 = gsn_csm_contour(wks,a_rain_s({lat|-7:17},month|:),res3)            ; default hovmueller

   
   overlay(plotA0,plotC0)
   overlay(plotA0,plotB0)
 
;
   plotA1 = gsn_csm_contour(wks,d_GPCPb({lat|-7:17},month|:),res)            ; default hovmueller
   plotB1 = gsn_csm_contour(wks,d_GPCPbd({lat|-7:17},month|:),res1)            ; default hovmueller
   ;plotB1 = gsn_contour_shade(plotB0,5.,p95,opt)
   plotB1 = gsn_contour_shade(plotB1,max(d_GPCPb),min(d_GPCPb),opt)
   plotC1 = gsn_csm_contour(wks,a_rainb_s({lat|-7:17},month|:),res3)            ; default hovmueller

   
   overlay(plotA1,plotC1)
   overlay(plotA1,plotB1)
 
;==================================================================================
   r0               = True                        ; plot mods desired
   r0@gsnDraw       = False                       ; do not draw
   r0@gsnFrame      = False                       ; do not advance frame
   r0@tiMainString  = ""                       ; title
   r0@tiYAxisString = "ITCZ feedback (~S~o~N~N)"                  ; yaxis
   r0@tiXAxisString = "Calendar month"                  ; yaxis
   r0@tiXAxisFontHeightF = 0.020
   r0@tiYAxisFontHeightF = 0.020
   r0@tmXBLabelFontHeightF = 0.020
   r0@tmYLLabelFontHeightF = 0.020
   r0@trYLog              = False                ; log scaling
   r0@trXLog              = False               ; log scaling
   r0@trYMinF             =   -1.5                 ; manually set lower limit
   r0@trYMaxF             =    1.5           

   r0@trXMinF             =  0.                 ; manually set lower limit
   r0@trXMaxF             =  11.                 ;   "          upper
  
   r0@xyLineColors        =  xy_colors_dheating
   r0@tmXTBorderOn            =  True
   r0@tmXTOn                  =  True
      
   r0@tmYLFormat = "0@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
 
   
   r0@vpWidthF                = 0.80
   r0@vpHeightF               = 0.45
     r0@xyLineThicknesses   = (/1.5,5.5,1.5,1.5,1.5,1.5/)       ; Define line thicknesses 
   r0@xyDashPatterns      = (/16.,16.,16.,0.,0.,0./)           ; Dash patterns 
   r0@tmYLLabelFontColor = "black"   ; and right axes are same.
   r0@tiYAxisFontColor  = "black"
   
 
   r0@tmXBMode = "Explicit"
   r0@tmXBValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
    r0@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)

     
  ;--------
  r0@tmXMajorGrid                = False          ; implement x grid 
  r0@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
  r0@tmXMajorGridLineDashPattern = 2             ; select short dash lines
  r0@tmXMajorGridLineColor	 = "black"
  
  r0@tmYMajorGrid                = False          ; implement y grid
  r0@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
  r0@tmYMajorGridLineDashPattern = 2             ; select short dash lines
  r0@tmYMajorGridLineColor	 = "black"

  r0@pmLegendDisplayMode    = "Always"            ; turn on legend
  r0@pmLegendSide           = "Top"               ; Change location of 
  r0@pmLegendParallelPosF   = 0.86                 ; move units right
  r0@pmLegendOrthogonalPosF = -1.05                ; move units down
  r0@pmLegendWidthF         = 0.12                ; Change width and
  r0@pmLegendHeightF        = 0.10                ; height of legend.
  r0@lgLabelFontHeightF     = 0.014                 ; change font height
  r0@lgTitleOn              = False               ; turn on legend title
  r0@lgTitleString          = "SST indices"           ; create legend title
  r0@lgTitleFontHeightF     = .003                ; font of legend title
  r0@lgPerimOn              = False               ; turn off box around
  r0@xyExplicitLegendLabels = (/"Mean ITCZ ~F18~S~F21~2~F33~s~N","Mean ITCZ","ERAI","ERA5","CFSR","SODA331"/)         ; explicit labels

 
  r0@xyLineColors        =  (/"black","black","black","deeppink4","deeppink4","deeppink4"/)

  plotD0  =  gsn_csm_xy(wks,xlabel,itcz_2sigma(lev|0:1,time|:),r0)
  plotD1  =  gsn_csm_xy(wks,xlabel,itczb_2sigma(lev|0:1,time|:),r0)
  r0@pmLegendDisplayMode    = "Never"            ; turn on legend
  plotE0  =  gsn_csm_xy(wks,xlabel,itcz_2sigma(lev|2:2,time|:),r0)
  plotE1  =  gsn_csm_xy(wks,xlabel,itczb_2sigma(lev|2:2,time|:),r0)


 r0@pmLegendDisplayMode    = "Always"            ; turn on legend
;-----------------------------------------------------------------------------
 
 
  overlay(plotA0,plotE0)
  overlay(plotA0,plotD0)
  plot(0)  =   plotA0 
;
  overlay(plotA1,plotE1)
  overlay(plotA1,plotD1)
  plot(1)  =   plotA1 

;=============================================================================
;-----------------------------------------------------------------------------
;=====================================================================
;=====================================================================
nplots = 2
;-----------------------------------------------------------------------------
  txres               = True
  txres@txAngleF      = 0.          ; Rotate string clockwise 
  txres@txFontHeightF = 0.025   ; Use same font height as left axis
  txres@txFont                = "Helvetica-Bold"
 ; Move text string to center/right edge of plot.;
  amres                  = True
  amres@amParallelPosF   =  -0.65 ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.
  amres@amOrthogonalPosF =  -0.52    ; This is the center of the plot.
  amres@amJust           = "TopLeft"
  

 tt = (/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o"/)
create_tt = gsn_create_text(wks,tt,txres);
tt_p0 = new(9,graphic)                ; create graphic array ; only needed if paneling
;---Attach the letters
     do i = 0,nplots-1,1
      tt_p0(i) = gsn_add_annotation(plot(i),create_tt(i),amres)
     end do 
 delete(i)


DATA_NAMES = (/"28~S~o~N~W","18 - 35~S~o~N~W","ORAS4","ORAS5","SODA3","Residual term"/)
txres@txFontHeightF = 0.020   ; Use same font height as left axis
txres@txFont                = "Helvetica"
 amres@amParallelPosF   =  -0.47 ; 0.5 is the right edge of the plot, so
 amres@amOrthogonalPosF =   -0.40    ; This is the center of the plot.
create_tt2 = gsn_create_text(wks,DATA_NAMES,txres);
tt_p1 = new(9,graphic) 
;---Attach the letters
     do i = 0,nplots-1,1
      tt_p1(i) = gsn_add_annotation(plot(i),create_tt2(i),amres)
     end do 
 delete(i)
;======================================================================== 
;========================================================================

;===============================================================
  latEQ = (/ 3.,  3.,  -3.,   -3., 3./)
  lonEQ = (/-20., 0., 0., -20.,-20./)
;---Resources for the polyline
  pres0                  = True           ; polyline mods desired
  pres0@gsLineThicknessF = 2.            ; line thickness   
  pres0@gsLineColor      = "black"         ; color of lines
  pres0@gsLineDashPattern = 0
  pres0@tfPolyDrawOrder   = "postDraw"
 
 ;dum0 = gsn_add_polyline(wks,plot(2),lonEQ,latEQ,pres0)
; dum1 = gsn_add_polyline(wks,plot(3),lonEQ,latEQ,pres0)
; dum2 = gsn_add_polyline(wks,plot(2),lonEQ,latEQ,pres0)


;===============================================================
;===============================================================
  
 
 
  x_0 = 5.
  x_1 = 5.
  y_0 = -40. 
  y_1 = 40.  

  ;y_0 = r2@trYMinF 
  ;y_1 = r2@trYMaxF  

;---Resources for the polyline
  pres1                  = True           ; polyline mods desired
  pres1@gsLineThicknessF = 0.5            ; line thickness   
  pres1@gsLineColor      = "black"         ; color of lines
  pres1@gsLineDashPattern = 2
  pres1@tfPolyDrawOrder   = "postDraw"
  pres1@gsLineColor      = "black"         ; color of lines
  dumL0 = gsn_add_polyline(wks,plot(0), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
  dumL1 = gsn_add_polyline(wks,plot(1), (/x_0,x_1/),(/y_0,y_1/) ,pres1)

;
  x0 = 0.
  x1 = 11.
  y0 = 0.
  dumL0x = gsn_add_polyline(wks,plot(0), (/x0,x1/),(/y0,y0/) ,pres1)
  dumL1x = gsn_add_polyline(wks,plot(1), (/x0,x1/),(/y0,y0/) ,pres1)

;===============================================================
;===============================================================
;===============================================================
;===============================================================
;=====================================================================
;===============================================================
   resP             			  = True                        ; panel mods desired
   resP@gsnDraw              		  = True       ; don't draw yet
   resP@gsnFrame         		  = False        ; don't advance frame yet
   resP@gsnPanelYWhiteSpacePercent       = 0.0    ; default is 1.0
   resP@gsnPanelXWhiteSpacePercent       = 0.0    ; default is 1.0
   resP@gsnPanelTop                      = 0.70
   resP@gsnPanelBottom                   = 0.41
   resP@gsnPanelLeft                     = 0.10
   resP@gsnPanelRight                    = 0.90
   resP@gsnPanelFigureStringsPerimOn 	  = False 
   resP@gsnPaperOrientation 		  = "portrait"
   resP@gsnPanelDebug      		  = True
   resP@amJust   			  = "BottomLeft"
   resP@gsnPanelFigureStringsPerimOn     = False
   resP@gsnPanelFigureStringsFontHeightF = 0.010
   resP@gsnPanelFigureStringsBackgroundFillColor = "transparent"
   resP@gsnPanelFigureStrings= (/"lag = -1","lag = +1","c","d","e","f","g","h","i","j","k","l"/) ; add strings to panel
   resP@gsnMaximize = True                        ; blow up plot
 ; gsn_panel(wks,plot(2:),(/1,1/),resP)               ; create panel plots
;;*******************************************
;===============================================================
 
;===============================================================
;===============================================================
;***********************************************
; panel plot only resources
  res_P                     		 = True            ; modify the panel plot
  res_P@gsnDraw            		 = True       ; don't draw yet
  res_P@gsnFrame            		 = False        ; don't advance frame yet
  res_P@gsnPanelYWhiteSpacePercent 	= 0.0    ; default is 1.0
  res_P@gsnPanelXWhiteSpacePercent 	= 1.0    ; default is 1.0
  res_P@gsnPanelTop                      = 0.85
  res_P@gsnPanelBottom                   = 0.05
  res_P@gsnPanelLeft                     = 0.01
  res_P@gsnPanelRight                    = 0.99
  res_P@gsnPaperOrientation = "portrait"   ; force portrait
 ; res_P@gsnPanelFigureStrings   = (/"(b)","~F0~Peak~F34~;~F0~12 yr~S~-1~N~"/)
  res_P@amJust   			  = "TopLeft"
  res_P@gsnPanelFigureStringsFontHeightF = 0.020
  res_P@gsnPanelFigureStringsPerimOn 	  = True 
  res_P@gsnPanelFigureStringsBackgroundFillColor = "grey"
  ;res_P@gsnPanelFigureStringsJust    =     "CenterRight"
  ;res_P@amJust               =   "CenterRight"
  res_P@gsnMaximize          = True            ; large format
  gsn_panel(wks,plot(0:1),(/2,1/),res_P)     ; now draw as one plot
;======================================================
;======================================================
;===============================================================
;===============================================================
;***********************************************
;==============================================================
frame(wks)
end
