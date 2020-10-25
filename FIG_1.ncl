;*************************************************
;*************************************************
begin
;==================================================
latS   = -15.
latN   =  15. 
lonL   = -60.
lonR   =  20.


nmos = 12
xlabel = ispan(0,nmos-1,1)
print(xlabel)
  
 

mm = 5		;index counts from zero, June is 5th month

p95 = 0.05

mxlag = 7
lag = 1	;first-order autocorrelation
;
days = (/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./)
days = days*24*60*60	;day,hours,minutes in each month
print(days)



xy_colors_dheating        = (/"red","blue","deepskyblue","black","deeppink4","black"/)
;xy_colors_thermocline        = (/"blue","yellow2","red","midnightblue","deepskyblue","black"/)
xy_colors_thermocline        = (/"red","yellow2","deepskyblue","blue","deeppink4","black"/)




yrStrt = 1982
yrLast = 2010
;=================================================
;=================================================

 f0      = addfile("atl3_sst_diab_cfsr_1982_2010_basin.nc", "r")
  TIME   = f0->time    
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  cfsr    = f0->diab(iYYYY,:)

  printVarSummary(cfsr) 
  printMinMax(cfsr,True)

  cfsr!0 = "time"
  cfsr!1 = "lev"

 
     
  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set cfsre time dimension
  
;=============================================================================================================
;=============================================================================================================
;=========================================================================================
;=========================================================================================
;=========================================================================================
nyears = dimsizes(cfsr&time)/nmos
print(nyears)
;
cfsr_s     = onedtond(ndtooned(cfsr),(/nyears,nmos,2/))
cfsr_s!0 = "time"
cfsr_s!1 = "month"
cfsr_s!2 = "lev"
printVarSummary(cfsr_s)
printMinMax(cfsr_s,True)
;cfsr_s = cfsr_s/conform(cfsr_s,days,1)
print(cfsr_s(:,:,0))
;
delete(cfsr)
;---------------------------
;
temp = cfsr_s(:,:,0)
printVarSummary(temp)
temp = dtrend_n(temp,False,0)
;print(temp)
;print(xxxxxxxxxxxxx)

;
xx = cfsr_s
delete(cfsr_s)
;
;
cfsr_s = xx(:,:,1)
delete(xx)
printVarSummary(cfsr_s)
cfsr_s = dtrend_n(cfsr_s,False,0)
cfsr_s = cfsr_s/conform(cfsr_s,days,1)
;print(cfsr_s)
;print(xxxxxxxxx)
;
;
std_cfsr_s = dim_stddev_n_Wrap(cfsr_s,0)
print(std_cfsr_s)
;print(xxxxxxxxx)
std_cfsr = std_cfsr_s
print(std_cfsr)

print(std_cfsr_s)
;print(xxxxxxxxx)

;
cfsr_s = dim_standardize_n(cfsr_s,0,0)      ; normalize
;print(temp)
;print(xxxxxxxxxxx)
;=========================================================================================
;=========================================================================================
;
;
;d_cfsr = regCoef_n(cfsr_s(:,4),temp(:,:),0,0)
d_cfsr = regCoef(cfsr_s(month|:,time|:),temp(month|:,time|:))
d_cfsr!0 = "month"
print(d_cfsr)
;print(xxxxxxxxxxxxxxxxxx)
;
;d_cfsr = d_cfsr*1000.
print(d_cfsr)
;============================================================================================
;============================================================================================
;==============================================================
;==============================================================

;--------------test for significance

 lag_ts = esacr(temp(month|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(cfsr_s(month|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

  tval = onedtond(d_cfsr@tval,dimsizes(d_cfsr))

       df = floattointeger(dimsizes(cfsr_s(:,0))*(1. -(lag_ts(:,lag)*lag_ts_map(:,lag)))\
        /(1. + (lag_ts_map(:,lag)*lag_ts(:,lag)))) 

printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)

 prob@long_name = "probability"
   printVarSummary(prob)
;

d_cfsrd = d_cfsr
d_cfsrd = where(prob.le.p95,d_cfsrd,d_cfsrd@_FillValue) ; Significant TEMP
printVarSummary(d_cfsrd)
printMinMax(d_cfsrd,True)
print(d_cfsrd)

;delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(tval)
delete(df)
delete(prob)
delete(b)
;
;
delete(temp)
delete(cfsr_s)
;

;print(xxxxxxxxxxxx)

;=========================================================================================
;=========================================================================================
yrStrt = 1982
yrLast = 2015
;=================================================
;=================================================
;======================================================================================
;======================================================================================
f0      = addfile("a_sst__1982_2015.nc", "r")
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
;print(xxxxxxxxxxxxx)
;
temp = dim_standardize_n(a_sst_s,0,0)      ; normalize
copy_VarMeta(a_sst_s(:,:),temp) 
printVarSummary(temp)
print(temp(:,5))
;

;d_sst = regCoef_n(temp(:,5),a_sst_s(:,:),0,0)
d_sst = regCoef(temp(month|:,time|:),a_sst_s(month|:,time|:))
d_sst!0 = "month"
print(d_sst)
;print(xxxxxxxxxxxxx)
;===============================================================================================================
;===============================================================================================================
;============================================================================================
;==============================================================
;==============================================================
;--------------test for significance

  lag_ts = esacr(temp(month|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 
 lag_ts_map   = esacr(a_sst_s(month|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

  tval = onedtond(d_sst@tval,dimsizes(d_sst))

       df = floattointeger(dimsizes(a_sst_s(:,0))*(1. -(lag_ts(:,lag)*lag_ts_map(:,lag)))\
        /(1. + (lag_ts_map(:,lag)*lag_ts(:,lag)))) 

print(df)
;print(xxxxxxxxxx)
printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)

prob@long_name = "probability"
printVarSummary(prob)
print(prob)

;

d_sstd = d_sst
d_sstd = where(prob.lt.p95,d_sstd,d_sstd@_FillValue) ; Significant TEMP
printVarSummary(d_sstd)
printMinMax(d_sstd,True)
print(d_sstd)
;print(xxxxxxxxxxxxxx)
;delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(tval)
delete(df)
delete(prob)
delete(b)
;
delete(a_sst)
delete(a_sst_s)
delete(temp)
;print(xxxxxxxxxxx)
;=============================================================================================================
;=============================================================================================================
;===============================================================================================================
;===============================================================================================================
;===============================================================================================================
;===============================================================================================================
;=========================================================================================
;=========================================================================================
yrStrt = 1993
yrLast = 2015
;=================================================
;=================================================
;=================================================
;=================================================
 f0      = addfile("atl3_sst_aviso_1993_2017_basin.nc", "r")
  TIME   = f0->time    
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  th    = f0->th(iYYYY,:)

  printVarSummary(th) 
  printMinMax(th,True)

  th!0 = "time"
  th!1 = "lev"

 
     
  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set the time dimension
  
;=============================================================================================================
;=============================================================================================================
;=========================================================================================
;=========================================================================================
;=========================================================================================
nyears = dimsizes(th&time)/nmos
print(nyears)
;
th_s     = onedtond(ndtooned(th),(/nyears,nmos,2/))
th_s!0 = "time"
th_s!1 = "month"
th_s!2 = "lev"
printVarSummary(th_s)
delete(th)
;
th_s = dtrend_n(th_s,False,0)
;
std_th_s = dim_stddev_n_Wrap(th_s,0)
print(std_th_s(:,1))
std_aviso = std_th_s(:,1)*1000.
print(std_aviso)
;print(xxxxxxxxx)
;
temp = th_s(:,:,0)      ; normalize
copy_VarMeta(th_s(:,:,0),temp) 
printVarSummary(temp)
;print(temp(:,5))
;
;
;
xx = th_s(:,:,1)
delete(th_s)
th_s = xx
delete(xx)
th_s = dim_standardize_n(th_s,0,0)      ; normalize
;copy_VarMeta(th_s(:,:),temp) 
printVarSummary(th_s)
;print(xxxxxxxxxx)
;=========================================================================================
;=========================================================================================
;d_aviso = regCoef_n(temp(:,5),th_s(:,:,1),0,0)
d_aviso = regCoef(th_s(month|:,time|:),temp(month|:,time|:))
d_aviso!0 = "month"
print(d_aviso)
;
;d_aviso = d_aviso*1000.
print(d_aviso)
;print(xxxxxxxxxxxxx)
;==============================================================

;============================================================================================
;==============================================================
;==============================================================

;--------------test for significance

 lag_ts = esacr(temp(month|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(th_s(month|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

  tval = onedtond(d_aviso@tval,dimsizes(d_aviso))

       df = floattointeger(dimsizes(th_s(:,0))*(1. -(lag_ts(:,lag)*lag_ts_map(:,lag)))\
        /(1. + (lag_ts_map(:,lag)*lag_ts(:,lag)))) 

printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)

 prob@long_name = "probability"
   printVarSummary(prob)
;

d_avisod = d_aviso
d_avisod = where(prob.le.p95,d_avisod,d_avisod@_FillValue) ; Significant TEMP
printVarSummary(d_avisod)
printMinMax(d_avisod,True)
print(d_avisod)

;delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(tval)
delete(df)
delete(prob)
delete(b)
;=========================================================================================
;==============================================================
delete(nyears)
delete(th_s)
delete(std_th_s)
delete(temp)
delete(yrStrt)
delete(yrLast)
;print(xxxxxxxxxxxxxxxxxx)




yrStrt = 1982
yrLast = 2015


;===============================================================================================================
;===============================================================================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================

 f0      = addfile("/home/hyacinth/work/datasets/sst/OISST_sst.mnmean.nc", "r")
  
  TIME   = f0->time    
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  oisst    = short2flt(f0->sst(iYYYY,:,:))
 
  printVarSummary(oisst) 
  printMinMax(oisst,True)

  ;gradLatLon  = grad_latlon_cfd(oisst,oisst&lat,oisst&lon,True,False)
  ;delete(oisst)                           
 ; T_grad_lat  = gradLatLon[0]    ; for clarity; explicitly extract variables from returned 'list' 
  ;oisst  = gradLatLon[1]
  ;delete(gradLatLon)


;==============================================================
; compute climatology and Anomalies
;==============================================================
  oisstClm   = clmMonTLL(oisst)                             ; (12,lat,lon)
  printVarSummary(oisstClm)

  oisst  = calcMonAnomTLL(oisst,oisstClm)                ; (time, lat,lon)
  printVarSummary(oisst)         
  printMinMax(oisst,True)
;=================================================================
 ;  oisst = oisst - runave_n_Wrap (oisst,nave,opt,dim)   

   oisst = lonFlip(oisst)

   ;avg_OISST = wgt_areaave_Wrap(oisst(time |:, {lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0)

   ;avg_OISST = rmAnnCycle1D(avg_OISST)

   ;avg_OISST  = avg_OISST - runave(avg_OISST,nave,1)
   
   delete(f0)
   delete(TIME) 
   delete(YYYY)
   delete(iYYYY)            ; Re-set the time dimension
;=============================================================================================================
;=============================================================================================================
 ;=========================================================================================
;=============================================================================================================
;=============================================================================================================
;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
 oisst_d = new(dimsizes(oisst),typeof(oisst),oisst@_FillValue)
 x   = oisst&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =oisst(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    oisst_d(tt::nmos,:,:) = work
 end do

 printVarSummary(oisst_d)
 oisst = oisst_d

;===
delete([/tt,x,x1,oisst_d,work/])

;---------------------------------------------------------------------------
;
;aOISST = wgt_areaave_Wrap(oisst(time |:, {lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0)
;system("/bin/rm -f aOISST.nc")   ; remove any pre-existing file
;ncdf = addfile("aOISST.nc" ,"c")  ; open output netCDF file
;ncdf->aOISST  = aOISST                         ; 4D              
;
;print(xxxxxxxxxxxxxxxxx)
;
;---------------------------------------------------------------------------
nyears = dimsizes(oisst&time)/nmos
print(nyears)
oisst_s     = onedtond(ndtooned(oisst),(/nyears,nmos,dimsizes(oisst&lat),dimsizes(oisst&lon)/))
copy_VarCoords(oisst(0,:,:),oisst_s(0,0,:,:))                             ; copy coordinate arrays from slp to newslp
oisst_s!0 = "time"
oisst_s!1 = "month"
;
printVarSummary(oisst_s)

;
aOISST_s = wgt_areaave_Wrap(oisst_s(time|:,month|:, {lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0)	;June 
printVarSummary(aOISST_s)

std_aOISST_s = dim_stddev_n_Wrap(aOISST_s,0)
print(std_aOISST_s)


aOISST = wgt_areaave_Wrap(oisst_s(time|:,month|mm, {lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0)	;June 
printVarSummary(aOISST)
aOISST = dim_standardize_n(aOISST,0,0)      ; normalize 
;
aOISST3 = wgt_areaave_Wrap(oisst_s(time|:,month|3, {lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0)	;June 
printVarSummary(aOISST3)
aOISST3 = dim_standardize_n(aOISST3,0,0)      ; normalize 
;
;print(xxxxxxxxxxxxxxxxxxxxxxxxxx)
;==========================================================================================
;==========================================================================================
d_OISST = regCoef_n(aOISST,oisst_s(:,mm,:,:),0,0)
copy_VarCoords(oisst_s(0,0,:,:),d_OISST)
printVarSummary(d_OISST)
;
d_OISST3 = regCoef_n(aOISST3,oisst_s(:,3,:,:),0,0)
copy_VarCoords(oisst_s(0,0,:,:),d_OISST3)
printVarSummary(d_OISST3)
;==============================================================
;==============================================================
;==============================================================
;==============================================================

;--------------test for significance
 d_OISSTd = d_OISST
 
 df = d_OISST

 lag_ts = esacr(aOISST,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(oisst_s(lat|:,lon|:,month|mm,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(d_OISST@tval,dimsizes(d_OISST))
      
    ;do ne = 0,neof-1,1
           df(:,:) = floattointeger(dimsizes(oisst_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)


   
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printVarSummary(prob)


d_OISSTd = where(abs(prob).le.p95,d_OISSTd,d_OISSTd@_FillValue) ; Significant TEMP
printVarSummary(d_OISSTd)

delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(oisst_s)
delete(nyears)
;===============================================================================
;===============================================================================








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

rain    = short2flt(f->precip(iYYYY,{latS:latN},:))
rain = lonFlip(rain)
printVarSummary(rain)       ; variable overview

printMinMax(rain,True)

rain!1  =  "lat"
rain!2  =  "lon"
delete([/YYYY,iYYYY,f,TIME/])
;==============================================================
;===========================
;=============================================================================================================
 rain_d = new(dimsizes(rain),typeof(rain),rain@_FillValue)
 x   = rain&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =rain(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    rain_d(tt::nmos,:,:) = work
 end do

 printVarSummary(rain_d)
 rain = rain_d

;===
delete([/tt,x,x1,rain_d,work/])

;---------------------------------------------------------------------------

;
;=========================================================================================
;=========================================================================================
;==========================================================================================
;==========================================================================================
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
d_GPCP = regCoef_n(aOISST,rain_s(:,:,:,:),0,0)
copy_VarCoords(rain_s(0,:,:,:),d_GPCP)
printVarSummary(d_GPCP)
;
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;--------------test for significance
 d_GPCPd = d_GPCP
 
 df = d_GPCP

 lag_ts = esacr(aOISST,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(rain_s(month|:,lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(d_GPCP@tval,dimsizes(d_GPCP))
      
    ;do ne = 0,neof-1,1
           df(:,:,:) = floattointeger(dimsizes(rain_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,:,lag)))\
        /(1. + (lag_ts_map(:,:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)


   
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printVarSummary(prob)


d_GPCPd = where(abs(prob).le.p95,d_GPCPd,d_GPCPd@_FillValue) ; Significant TEMP
printVarSummary(d_GPCPd)

delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(rain_s)
delete(rain)
delete(nyears)
;===============================================================================
;===============================================================================














;===============================================================================================================
;===============================================================================================================
; Open the file: Read only the user specified period
; ==============================================================
f      = addfile("/home/hyacinth/work/datasets/ens_reanalysis/REMAPPED_ssh_CMEMS_altimetry_1993_2015.nc", "r")
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100  
printVarSummary(YYYY)               ; entire file
iYYYY  = ind(YYYY.ge.yrStrt.and.YYYY.le.yrLast)
printVarSummary(iYYYY)
aviso    = f->ssh(iYYYY,:,:)

printVarSummary(aviso)       ; variable overview

printMinMax(aviso,True)

aviso!1  =  "lat"
aviso!2  =  "lon"
;==============================================================
; compute climatology and Anomalies
;==============================================================
  avisoClm   = clmMonTLL(aviso)                             ; (12,lat,lon)
  printVarSummary(avisoClm)

  aviso  = calcMonAnomTLL(aviso,avisoClm)                ; (time, lat,lon)
  printVarSummary(aviso)         
  printMinMax(aviso,True)
;=================================================================
 ;  aviso = aviso - runave_n_Wrap (aviso,nave,opt,dim)   

   aviso = lonFlip(aviso)


   delete(f)
   delete(TIME) 
   delete(YYYY)
   delete(iYYYY)            ; Re-set the time dimension
;=============================================================================================================
;=============================================================================================================
 ;=========================================================================================
;=============================================================================================================
;=============================================================================================================
;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
nyears = dimsizes(aviso&time)/nmos
 aviso_d = new(dimsizes(aviso),typeof(aviso),aviso@_FillValue)
 x   = aviso&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =aviso(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    aviso_d(tt::nmos,:,:) = work
 end do

 printVarSummary(aviso_d)
 aviso = aviso_d

;===
delete([/tt,x,x1,aviso_d,work/])
;=========================================================================================
;=========================================================================================
;
nyears = dimsizes(aviso&time)/nmos
print(nyears)
aviso_s     = onedtond(ndtooned(aviso),(/nyears,nmos,dimsizes(aviso&lat),dimsizes(aviso&lon)/))
copy_VarCoords(aviso(0,:,:),aviso_s(0,0,:,:))                             ; copy coordinate arrays from slp to newslp
aviso_s!0 = "time"
aviso_s!1 = "month"
printVarSummary(aviso_s)
delete(nyears)
;
;==========================================================================================
;==========================================================================================
d_AVISO = regCoef_n(aOISST(11:),aviso_s(:,5,:,:),0,0)
copy_VarCoords(aviso_s(0,0,:,:),d_AVISO)
printVarSummary(d_AVISO)
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;--------------test for significance
 d_AVISOd = d_AVISO
 
 df = d_AVISO

 lag_ts = esacr(aOISST,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(aviso_s(lat|:,lon|:,month|mm,time|11:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(d_AVISO@tval,dimsizes(d_AVISO))
      
    ;do ne = 0,neof-1,1
           df(:,:) = floattointeger(dimsizes(aviso_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)


   
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printVarSummary(prob)


d_AVISOd = where(prob.le.p95,d_AVISOd,d_AVISOd@_FillValue) ; Significant TEMP
printVarSummary(d_AVISOd)

delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(aviso_s)
delete(aviso)
;===============================================================================
;===============================================================================
;=========================================================================================









;================================================================
; ===============================================================
; Open the file: Read only the user specified period
; ==============================================================
f      = addfile("/home/hyacinth/work/datasets/ens_reanalysis/ustress_merra_jra55_erai_era5_cfsr_1982_2018.nc", "r")
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100  
iYYYY  = ind(YYYY.ge.yrStrt.and.YYYY.le.yrLast)
ustr    = f->ustress(iYYYY,:,:,:)

printVarSummary(ustr)       ; variable overview


ustr = lonFlip(ustr)
printVarSummary(ustr)               ; entire file
printMinMax(ustr,True)
   
   delete(f)
   delete(TIME) 
   delete(YYYY)
   delete(iYYYY)            ; Re-set the time dimension
;=============================================================================================================
;=============================================================================================================
 ;=========================================================================================
;=============================================================================================================
;=============================================================================================================
;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
nyears = dimsizes(ustr&time)/nmos
print(nyears)
 ustr_d = new(dimsizes(ustr),typeof(ustr),ustr@_FillValue)
 x   = ustr&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =ustr(tt::nmos,:,:,:)
    work = dtrend_msg_n(x1,work,False,False,0)
    ustr_d(tt::nmos,:,:,:) = work
 end do

 printVarSummary(ustr_d)
 ustr = ustr_d

;===
delete([/tt,x,x1,ustr_d,work/])
;=========================================================================================
;=========================================================================================
;
nyears = dimsizes(ustr&time)/nmos
print(nyears)
;
ustr_s     = onedtond(ndtooned(ustr),(/nyears,nmos,dimsizes(ustr&lev),dimsizes(ustr&lat),dimsizes(ustr&lon)/))
printVarSummary(ustr_s)

copy_VarCoords(ustr(0,:,:,:),ustr_s(0,0,:,:,:))                             ; copy coordinate arrays from slp to newslp
ustr_s!0 = "time"
ustr_s!1 = "month"
printVarSummary(ustr_s)
;
;
meanU = dim_avg_n_Wrap(ustr_s,(/0/))
printVarSummary(meanU)
ens_meanU = dim_avg_n_Wrap(meanU,(/1/))
printVarSummary(ens_meanU)
delete(meanU)

;
;==========================================================================================
;==========================================================================================
d_USTR = regCoef_n(aOISST,ustr_s(:,:,:,:,:),0,0)
copy_VarCoords(ustr_s(0,:,:,:,:),d_USTR)
printVarSummary(d_USTR)
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;==============================================================

;--------------test for significance
 d_USTRd = d_USTR
 
 df = d_USTR

 lag_ts = esacr(aOISST,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(ustr_s(month|:,lev|:,lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(d_USTR@tval,dimsizes(d_USTR))
      
    ;do ne = 0,neof-1,1
           df(:,:,:,:) = floattointeger(dimsizes(ustr_s(:,0,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,:,:,lag)))\
        /(1. + (lag_ts_map(:,:,:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   probU = betainc(df/(df+tval^2),df/2.0,b)       ; probU(nlat,nlon)


   
   probU@long_name   = "regression coefficient"
   probU@long_name = "probability"
   printVarSummary(probU)



;d_USTRd = where(prob.le.p95,d_USTRd,d_USTRd@_FillValue) ; Significant
;printVarSummary(d_USTRd)
printMinMax(probU,True)

delete([/tval,df,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(ustr_s)
delete(ustr)
delete(nyears)



;================================================================
; ===============================================================
; Open the file: Read only the user specified period
; ==============================================================
f      = addfile("/home/hyacinth/work/datasets/ens_reanalysis/vstress_merra_jra55_erai_era5_cfsr_1982_2018.nc", "r")
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100  
iYYYY  = ind(YYYY.ge.yrStrt.and.YYYY.le.yrLast)
vstr    = f->vstress(iYYYY,:,:,:)

printVarSummary(vstr)       ; variable overview


vstr = lonFlip(vstr)
printVarSummary(vstr)               ; entire file
printMinMax(vstr,True)
   
   delete(f)
   delete(TIME) 
   delete(YYYY)
   delete(iYYYY)            ; Re-set the time dimension
;=============================================================================================================
;=============================================================================================================
 ;=========================================================================================
;=============================================================================================================
;=============================================================================================================
;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
nyears = dimsizes(vstr&time)/nmos
print(nyears)
 vstr_d = new(dimsizes(vstr),typeof(vstr),vstr@_FillValue)
 x   = vstr&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =vstr(tt::nmos,:,:,:)
    work = dtrend_msg_n(x1,work,False,False,0)
    vstr_d(tt::nmos,:,:,:) = work
 end do

 printVarSummary(vstr_d)
 vstr = vstr_d

;===
delete([/tt,x,x1,vstr_d,work/])
;=========================================================================================
;=========================================================================================
;
nyears = dimsizes(vstr&time)/nmos
print(nyears)
;
vstr_s     = onedtond(ndtooned(vstr),(/nyears,nmos,dimsizes(vstr&lev),dimsizes(vstr&lat),dimsizes(vstr&lon)/))
printVarSummary(vstr_s)

copy_VarCoords(vstr(0,:,:,:),vstr_s(0,0,:,:,:))                             ; copy coordinate arrays from slp to newslp
vstr_s!0 = "time"
vstr_s!1 = "month"
printVarSummary(vstr_s)
;
;
meanV = dim_avg_n_Wrap(vstr_s,(/0/))
printVarSummary(meanV)
ens_meanV = dim_avg_n_Wrap(meanV,(/1/))
printVarSummary(ens_meanV)
delete(meanV)

;
;==========================================================================================
;==========================================================================================
d_VSTR = regCoef_n(aOISST,vstr_s(:,:,:,:,:),0,0)
copy_VarCoords(vstr_s(0,:,:,:,:),d_VSTR)
printVarSummary(d_VSTR)
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;==============================================================

;--------------test for significance
 d_VSTRd = d_VSTR
 
 df = d_VSTR

 lag_ts = esacr(aOISST,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(vstr_s(month|:,lev|:,lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(d_VSTR@tval,dimsizes(d_VSTR))
      
    ;do ne = 0,neof-1,1
           df(:,:,:,:) = floattointeger(dimsizes(vstr_s(:,0,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,:,:,lag)))\
        /(1. + (lag_ts_map(:,:,:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   probV = betainc(df/(df+tval^2),df/2.0,b)       ; probV(nlat,nlon)


   
   probV@long_name   = "regression coefficient"
   probV@long_name = "probability"
   printVarSummary(probV)



;d_VSTRd = where(prob.le.p95,d_VSTRd,d_VSTRd@_FillValue) ; Significant
;printVarSummary(d_VSTRd)
printMinMax(probV,True)

delete([/tval,df,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(vstr_s)
delete(vstr)
delete(nyears)

;================================================================
ens_probV = dim_avg_n_Wrap(probV,1)
printVarSummary(ens_probV)
ens_d_VSTR = dim_avg_n_Wrap(d_VSTR,1)
printVarSummary(ens_d_VSTR)

;
ens_probU = dim_avg_n_Wrap(probU,1)
printVarSummary(ens_probU)
ens_d_USTR = dim_avg_n_Wrap(d_USTR,1)
printVarSummary(ens_d_USTR)
;

     
      ens_d_USTR  = where(ens_probU.le.p95.or.ens_probV.le.p95,ens_d_USTR,ens_d_USTR@_FillValue) 
      printMinMax(ens_d_USTR,True)
       printVarSummary(ens_d_USTR)
;++++++++++++++++++++++++++++++++++++++
      ens_d_VSTR  = where(ens_probV.le.p95.or.ens_probU.le.p95,ens_d_VSTR,ens_d_VSTR@_FillValue) 
      printMinMax(ens_d_VSTR,True)
       printVarSummary(ens_d_VSTR)
;------------------------------------------------------------------
;print(xxxxxxxxxxx)










;================================================================
; ===============================================================
; Open the file: Read only the user specified period
; ==============================================================
f      = addfile("/home/hyacinth/work/datasets/ens_reanalysis/ssh_soda331_oras5_oras4_godas_gecco2_1982_2015.nc", "r")
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100  
iYYYY  = ind(YYYY.ge.yrStrt.and.YYYY.le.yrLast)
ssh    = f->ssh(iYYYY,:,:,:)
ssh = ssh*1000.
printVarSummary(ssh)       ; variable overview


;ssh = lonFlip(ssh)
printVarSummary(ssh)               ; entire file
printMinMax(ssh,True)
;==============================================================
; compute climatology and Anomalies
;==============================================================
  sshClm   = clmMonTLLL(ssh)                             ; (12,lat,lon)
  printVarSummary(sshClm)

  ssh  = calcMonAnomTLLL(ssh,sshClm)                ; (time, lat,lon)
  printVarSummary(ssh)         
  printMinMax(ssh,True)
;=================================================================
 ;  ssh = ssh - runave_n_Wrap (ssh,nave,opt,dim)   

   ssh = lonFlip(ssh)
   
   delete(f)
   delete(TIME) 
   delete(YYYY)
   delete(iYYYY)            ; Re-set the time dimension
;=============================================================================================================
;=============================================================================================================
 ;=========================================================================================
;=============================================================================================================
;=============================================================================================================
;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
nyears = dimsizes(ssh&time)/nmos
print(nyears)
 ssh_d = new(dimsizes(ssh),typeof(ssh),ssh@_FillValue)
 x   = ssh&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =ssh(tt::nmos,:,:,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    ssh_d(tt::nmos,:,:,:) = work
 end do

 printVarSummary(ssh_d)
 ssh = ssh_d

;===
delete([/tt,x,x1,ssh_d,work/])
;=========================================================================================
;=========================================================================================
;
nyears = dimsizes(ssh&time)/nmos
print(nyears)
;
ssh_s     = onedtond(ndtooned(ssh),(/nyears,nmos,dimsizes(ssh&lev),dimsizes(ssh&lat),dimsizes(ssh&lon)/))
printVarSummary(ssh_s)

copy_VarCoords(ssh(0,:,:,:),ssh_s(0,0,:,:,:))                             ; copy coordinate arrays from slp to newslp
ssh_s!0 = "time"
ssh_s!1 = "month"
printVarSummary(ssh_s)
;
;==========================================================================================
;==========================================================================================
d_SSH = regCoef_n(aOISST,ssh_s(:,5,:,:,:),0,0)
copy_VarCoords(ssh_s(0,0,:,:,:),d_SSH)
printVarSummary(d_SSH)
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;==============================================================

;--------------test for significance
 d_SSHd = d_SSH
 
 df = d_SSH

 lag_ts = esacr(aOISST,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(ssh_s(lev|:,lat|:,lon|:,month|mm,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(d_SSH@tval,dimsizes(d_SSH))
      
    ;do ne = 0,neof-1,1
           df(:,:,:) = floattointeger(dimsizes(ssh_s(:,0,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,:,lag)))\
        /(1. + (lag_ts_map(:,:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)


   
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printVarSummary(prob)



d_SSHd = where(prob.le.p95,d_SSHd,d_SSHd@_FillValue) ; Significant
printVarSummary(d_SSHd)
printMinMax(prob,True)

delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(ssh_s)
delete(ssh)
delete(nyears)

;===============================================================================
dimZ = dimsizes(d_SSH)
x = new((/6,dimZ(1),dimZ(2)/),typeof(d_SSH),d_SSH@_FillValue)
x(0:4,:,:) = d_SSH	
x(5,:,:) = d_AVISO	;cfsr
printVarSummary(x)
delete(d_SSH)
d_SSH = dim_avg_n_Wrap(x,0)
delete(x)
;-------------------------
x = new((/6,dimZ(1),dimZ(2)/),typeof(d_SSHd),d_SSHd@_FillValue)
x(0:4,:,:) = d_SSHd	
x(5,:,:) = d_AVISOd	;cfsr
printVarSummary(x)
delete(d_SSHd)
d_SSHd = dim_avg_n_Wrap(x,0)
delete(x)
;
;===============================================================================
;===============================================================================================================
;===============================================================================================================
;================================================================
; ===============================================================
; Open the file: Read only the user specified period
; ==============================================================
f      = addfile("/home/hyacinth/work/datasets/ens_reanalysis/REMAPPED_diab_Ddiabf01.gdas.1979_2011.grb2.nc", "r")
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100  
printVarSummary(YYYY)               ; entire file
iYYYY  = ind(YYYY.ge.yrStrt.and.YYYY.le.yrLast)
printVarSummary(iYYYY)
cfsr    = f->diab(iYYYY,:,:)

printVarSummary(cfsr)       ; variable overview

printMinMax(cfsr,True)

cfsr!1  =  "lat"
cfsr!2  =  "lon"
;==============================================================
; compute climatology and Anomalies
;==============================================================
  cfsrClm   = clmMonTLL(cfsr)                             ; (12,lat,lon)
  printVarSummary(cfsrClm)

  cfsr  = calcMonAnomTLL(cfsr,cfsrClm)                ; (time, lat,lon)
  printVarSummary(cfsr)         
  printMinMax(cfsr,True)
;=================================================================
 ;  cfsr = cfsr - runave_n_Wrap (cfsr,nave,opt,dim)   

   cfsr = lonFlip(cfsr)


   delete(f)
   delete(TIME) 
   delete(YYYY)
   delete(iYYYY)            ; Re-set the time dimension
;=============================================================================================================
;=============================================================================================================
 ;=========================================================================================
;=============================================================================================================
;=============================================================================================================
;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
 cfsr_d = new(dimsizes(cfsr),typeof(cfsr),cfsr@_FillValue)
 x   = cfsr&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =cfsr(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    cfsr_d(tt::nmos,:,:) = work
 end do

 printVarSummary(cfsr_d)
 cfsr = cfsr_d

;===
delete([/tt,x,x1,cfsr_d,work/])
;=========================================================================================
;=========================================================================================
;
nyears = dimsizes(cfsr&time)/nmos
print(nyears)
cfsr_s     = onedtond(ndtooned(cfsr),(/nyears,nmos,dimsizes(cfsr&lat),dimsizes(cfsr&lon)/))
copy_VarCoords(cfsr(0,:,:),cfsr_s(0,0,:,:))                             ; copy coordinate arrays from slp to newslp
cfsr_s!0 = "time"
cfsr_s!1 = "month"
printVarSummary(cfsr_s)
delete(nyears)
;

;
;==========================================================================================
;==========================================================================================
printMinMax(cfsr_s,True)
cfsr_s = cfsr_s/conform(cfsr_s,days,1)			;Converts to K
printMinMax(cfsr_s,True)
;==========================================================================================
;==========================================================================================
d_CFSR = regCoef_n(aOISST(0:28),cfsr_s(:,4,:,:),0,0)
copy_VarCoords(cfsr_s(0,0,:,:),d_CFSR)
printVarSummary(d_CFSR)
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;==============================================================

;--------------test for significance
 d_CFSRd = d_CFSR
 
 df = d_CFSR

 lag_ts = esacr(aOISST,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(cfsr_s(lat|:,lon|:,month|4,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(d_CFSR@tval,dimsizes(d_CFSR))
      
    ;do ne = 0,neof-1,1
           df(:,:) = floattointeger(dimsizes(cfsr_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)


   
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printVarSummary(prob)


d_CFSRd = where(prob.le.p95,d_CFSRd,d_CFSRd@_FillValue) ; Significant TEMP
printVarSummary(d_CFSRd)

delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(cfsr_s)
delete(cfsr)
;===============================================================================
;===============================================================================

;================================================================
; ===============================================================
; Open the file: Read only the user specified period
; ==============================================================
f      = addfile("/home/hyacinth/work/datasets/ens_reanalysis/diab_merra_jra55_erai_era5_1982_2018.nc", "r")
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100  
iYYYY  = ind(YYYY.ge.yrStrt.and.YYYY.le.yrLast)
qt    = f->diab(iYYYY,:,:,:)
printVarSummary(qt)       ; variable overview


;qt = lonFlip(qt)
printVarSummary(qt)               ; entire file
printMinMax(qt,True)
;==============================================================
; compute climatology and Anomalies
;==============================================================
  qtClm   = clmMonTLLL(qt)                             ; (12,lat,lon)
  printVarSummary(qtClm)

  qt  = calcMonAnomTLLL(qt,qtClm)                ; (time, lat,lon)
  printVarSummary(qt)         
  printMinMax(qt,True)
;=================================================================
 ;  qt = qt - runave_n_Wrap (qt,nave,opt,dim)   

   qt = lonFlip(qt)
   
   delete(f)
   delete(TIME) 
   delete(YYYY)
   delete(iYYYY)            ; Re-set the time dimension
;=============================================================================================================
;=============================================================================================================
 ;=========================================================================================
;=============================================================================================================
;=============================================================================================================
;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
 qt_d = new(dimsizes(qt),typeof(qt),qt@_FillValue)
 x   = qt&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =qt(tt::nmos,:,:,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    qt_d(tt::nmos,:,:,:) = work
 end do

 printVarSummary(qt_d)
 qt = qt_d

;===
delete([/tt,x,x1,qt_d,work/])
;=========================================================================================
;=========================================================================================
;
nyears = dimsizes(qt&time)/nmos
print(nyears)
;
qt_s     = onedtond(ndtooned(qt),(/nyears,nmos,dimsizes(qt&lev),dimsizes(qt&lat),dimsizes(qt&lon)/))
printVarSummary(qt_s)

copy_VarCoords(qt(0,:,:,:),qt_s(0,0,:,:,:))                             ; copy coordinate arrays from slp to newslp
qt_s!0 = "time"
qt_s!1 = "month"
printVarSummary(qt_s)
;
;==========================================================================================
printMinMax(qt_s,True)				;converts to K/s
qt_s = qt_s/conform(qt_s,days,1)
printMinMax(qt_s,True)
;==========================================================================================
d_QT = regCoef_n(aOISST,qt_s(:,4,:,:,:),0,0)
copy_VarCoords(qt_s(0,0,:,:,:),d_QT)
printVarSummary(d_QT)
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;==============================================================
;==============================================================

;--------------test for significance
 d_QTd = d_QT
 
 df = d_QT

 lag_ts = esacr(aOISST,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(qt_s(lev|:,lat|:,lon|:,month|4,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(d_QT@tval,dimsizes(d_QT))
      
    ;do ne = 0,neof-1,1
           df(:,:,:) = floattointeger(dimsizes(qt_s(:,0,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,:,lag)))\
        /(1. + (lag_ts_map(:,:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)


   
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printVarSummary(prob)



d_QTd = where(prob.le.p95,d_QTd,d_QTd@_FillValue) ; Significant
printVarSummary(d_QTd)
printMinMax(prob,True)

delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(qt_s)
delete(qt)



;===============================================================================
dimZ = dimsizes(d_QT)
x = new((/5,dimZ(1),dimZ(2)/),typeof(d_QT),d_QT@_FillValue)
x(0:3,:,:) = d_QT	
x(4,:,:) = d_CFSR	;cfsr
printVarSummary(x)
delete(d_QT)
d_QT = dim_avg_n_Wrap(x,0)
delete(x)
;-------------------------
x = new((/5,dimZ(1),dimZ(2)/),typeof(d_QTd),d_QTd@_FillValue)
x(0:3,:,:) = d_QTd	
x(4,:,:) = d_CFSRd	;cfsr
printVarSummary(x)
delete(d_QTd)
d_QTd = dim_avg_n_Wrap(x,0)
delete(x)
;
;===============================================================================
d_QTd = d_QTd*1000000.
d_QT = d_QT*1000000.
;===============================================================================
;===============================================================================
;===============================================================================
;yrStrt = 1982
;yrLast = 2015
;
;f0      = addfile("atl3_ssh_soda342_oras5_oras4_godas_gecco2_1982_2015_basin.nc", "r")
f0      = addfile("/home/hyacinth/work/geomar/sstP/R0/a_ssh_soda342_oras5_oras4_godas_gecco2_1982_2015_basin.nc", "r")

TIME   = f0->time    
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
tcline    = f0->ssh(iYYYY,:)
printVarSummary(tcline) 
printMinMax(tcline,True)

 tcline  = tcline*1000.

  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set the time dimension

;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
 tcline_d = new(dimsizes(tcline),typeof(tcline),tcline@_FillValue)
 x   = tcline&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =tcline(tt::nmos,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    tcline_d(tt::nmos,:) = work
 end do

 printVarSummary(tcline_d)
 tcline = tcline_d

;===
delete([/tt,x,x1,tcline_d,work/])
;=========================================================================================
;=========================================================================================
;
tcline_s     = onedtond(ndtooned(tcline),(/nyears,nmos,5/))
tcline_s!0 = "time"
tcline_s!1 = "month"
tcline_s!2 = "lev"
printVarSummary(tcline_s)
delete(tcline)
;
tcline_s = dtrend_n(tcline_s,False,0)
;
std_tcline_s = dim_stddev_n_Wrap(tcline_s,0)
print(std_tcline_s)
;
x = new((/nmos,6/),float,std_tcline_s@_FillValue)
x(:,0:4) = std_tcline_s	
x(:,5) = std_aviso	;dynamic topography
printVarSummary(x)
print(x(:,5))
delete(std_tcline_s)
std_tcline_s = x
delete(x)
;
std_tcline_s!0 = "month"
std_tcline_s!1 = "lev"
;=========================================================================================
;=========================================================================================
;
;tcline_s = dim_standardize_n(tcline_s,0,0)      ; normalize 
;=========================================================================================
;=========================================================================================
printVarSummary(aOISST_s)
printVarSummary(tcline_s)
;
tConform = conform_dims(dimsizes(tcline_s),aOISST_s,(/0,1/)) 
printVarSummary(tConform)
delete(aOISST_s)
copy_VarMeta(tcline_s,tConform) 
aOISST_s = tConform
delete(tConform)
printVarSummary(aOISST_s)
;
;
;
tcline_s = dim_standardize_n(tcline_s,0,0)      ; normalize
;
;d_tcline = regCoef_n(aOISST,tcline_s,0,0)
d_tcline = regCoef(tcline_s(month|:,lev|:,time|:),aOISST_s(month|:,lev|:,time|:))
d_tcline!0 = "month"
d_tcline!1 = "lev"
print(d_tcline)
print(d_tcline(:,4))
;
;print(xxxxxxxxxxxx)
;r_tcline = regCoef_n(aOISST,tcline_s,0,0)
;r_tcline!0 = "month"
;r_tcline!1 = "lev"
;==============================================================
;==============================================================

;--------------test for significance

 lag_ts = esacr(aOISST_s(month|:,lev|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(tcline_s(month|:,lev|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

  tval = onedtond(d_tcline@tval,dimsizes(d_tcline))

       df = floattointeger(dimsizes(tcline_s(:,0,0))*(1. -(lag_ts(:,:,lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(:,:,lag)))) 

printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)

 prob@long_name = "probability"
   printVarSummary(prob)

   

;

d_tclined = d_tcline
d_tclined = where(prob.le.p95,d_tclined,d_tclined@_FillValue) ; Significant TEMP
printVarSummary(d_tclined)
printMinMax(d_tclined,True)

;delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(tcline_s)
delete(tval)
delete(df)
delete(prob)
delete(b)
;=========================================================================================
;=========================================================================================
;=========================================================================================
;=========================================================================================
x = new((/nmos,6/),float,d_tcline@_FillValue)
x(:,0:4) = d_tcline	
x(:,5) = d_aviso	;dynamic topography
printVarSummary(x)
print(x(:,5))
delete(d_tcline)
d_tcline = x
delete(x)
;
;
x = new((/nmos,6/),float,d_tcline@_FillValue)
x(:,0:4) = d_tclined	
x(:,5) = d_avisod	;dynamic topography
printVarSummary(x)
print(x(:,5))
delete(d_tclined)
d_tclined = x
delete(x)
;print(xxxxxxxxxxxxxxxxx)
;=========================================================================================
;=========================================================================================
;============================================================================================
;=========================================================================================
;=========================================================================================
;
;yrStrt = 1982
;yrLast = 2010
;

;f0      = addfile("a_diab_merra_jra55_erai_era5_cfsr_1982_2010.nc", "r")
f0      = addfile("atl3_diab_merra_jra55_erai_era5_1982_2018_basin.nc", "r")
TIME   = f0->time    
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
diab    = f0->diab(iYYYY,0:3)
printVarSummary(diab) 
printMinMax(diab,True)


  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set the time dimension

;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
 diab_d = new(dimsizes(diab),typeof(diab),diab@_FillValue)
 x   = diab&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =diab(tt::nmos,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    diab_d(tt::nmos,:) = work
 end do

 printVarSummary(diab_d)
 diab = diab_d

;===
delete([/tt,x,x1,diab_d,work/])
;=========================================================================================
;=========================================================================================
;
diab_s     = onedtond(ndtooned(diab),(/nyears,nmos,4/))
diab_s!0 = "time"
diab_s!1 = "month"
diab_s!2 = "lev"
printVarSummary(diab_s)
delete(diab)


;
printMinMax(diab_s,True)
diab_s = diab_s/conform(diab_s,days,1)
printMinMax(diab_s,True)
;print(xxxxxxxxxxxxx)
;
diab_s = dtrend_n(diab_s,False,0)
;
std_diab_s = dim_stddev_n_Wrap(diab_s,0)
print(std_diab_s(month|:,lev|1))
printVarSummary(std_diab_s)
;print(xxxxxxxxxxxxx)
;
;diab_s = dim_standardize_n(diab_s,0,0)      ; normalize 
;=========================================================================================
;
;
;
diab_s = dim_standardize_n(diab_s,0,0)      ; normalize
;
;=========================================================================================
r_diab = regCoef(diab_s(month|:,lev|:,time|:),aOISST_s(month|:,lev|0:3,time|:))
;r_diab = regCoef_n(aOISST,diab_s,0,0)
r_diab!0 = "month"
r_diab!1 = "lev"
printVarSummary(r_diab)
;print(xxxxxxxxxxxxxx)
;============================================================================================
;==============================================================
;==============================================================

;--------------test for significance

 lag_ts = esacr(aOISST_s(month|:,lev|0:3,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(diab_s(month|:,lev|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

  tval = onedtond(r_diab@tval,dimsizes(r_diab))

       df = floattointeger(dimsizes(diab_s(:,0,0))*(1. -(lag_ts(:,:,lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(:,:,lag)))) 

printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)

 prob@long_name = "probability"
   printVarSummary(prob)

   

;

r_diabd = r_diab
r_diabd = where(prob.le.p95,r_diabd,r_diabd@_FillValue) ; Significant TEMP
printVarSummary(r_diabd)
printMinMax(r_diabd,True)



mean_sig = dim_avg_n_Wrap(abs(prob),1)
print(mean_sig)
mean_diab = dim_avg_n_Wrap(r_diab,1)
mean_diab = mean_diab*1000000.
print(mean_diab)
c_diab = dim_stddev_n_Wrap(r_diab,1)
c_diab = c_diab*1000000.
print(c_diab)
;
delete(mean_sig)
delete(mean_diab)
delete(c_diab)
;print(xxxxxxxxxx)

;delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(diab_s)
delete(tval)
delete(df)
delete(prob)
delete(b)

;=========================================================================================
;=========================================================================================
x = new((/nmos,5/),float,d_tcline@_FillValue)
x(:,0:3) = (/r_diab/)	
x(:,4) = (/d_cfsr/)		;ncep cfsr
printVarSummary(x)
print(x(:,4))
delete(r_diab)
r_diab = x
delete(x)
r_diab!0 = "month"
r_diab!1 = "lev"
printVarSummary(r_diab)
;
;
x = new((/nmos,5/),float,r_diab@_FillValue)
x(:,0:3) = (/r_diabd/)	
x(:,4) = (/d_cfsrd/)	;ncep cfsr
printVarSummary(x)
print(x(:,4))
delete(r_diabd)
r_diabd = x
delete(x)

r_diabd!0 = "month"
r_diabd!1 = "lev"
printVarSummary(r_diabd)
;
;
;===============================================================
;
;
x = new((/nmos,5/),float,std_diab_s@_FillValue)
x(:,0:3) = (/std_diab_s/)	
x(:,4) = (/std_cfsr_s/)	;ncep cfsr
printVarSummary(x)
print(x(:,4))
delete(std_diab_s)
std_diab_s = x
delete(x)

std_diab_s!0 = "month"
std_diab_s!1 = "lev"
printVarSummary(std_diab_s)
;
;=========================================================================================
;=========================================================================================
;===============================================================================================================
;===============================================================================================================
;==========================================================================================
;==========================================================================================
print("==============================================================")
print("==============================================================")
;===============================================================================================================
;======================================
  wks = gsn_open_wks("eps","FIG_1")         ; send graphics to PNG file
  plot = new(6,graphic)                ; create graphic array ; only needed if paneling
  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res@cnMinLevelValF       = -30.            ; set min contour level
  res@cnMaxLevelValF       =  30.             ; set max contour level
  res@cnLevelSpacingF      =  3.            ; set contour spacing   
  res@cnFillOn             = True         ; turn on color fill
  res@cnFillPalette        = "cmp_b2r"  
  ;res@cnFillPalette        = "BlueWhiteOrangeRed"  
  ;res@cnFillPalette        = "ncl_default"  
  ;res@cnFillPalette        = "gui_default" ; set color map
 ; res@cnFillPalette        = "so4_23"
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@cnInfoLabelOn        = False
  res@lbLabelBarOn         = True       ; turn off individual lb's 
  ;res@pmLabelBarParallelPosF = 0.50
  ;res@pmLabelBarOrthogonalPosF = -0.11
  res@lbOrientation        = "vertical"
  res@lbTitleFontHeightF  = 0.020 
  res@lbLabelFontHeightF  = 0.020
  res@lbLeftMarginF       = 0.10
  res@lbRightMarginF      = 0.55
  ;res@lbLabelOffsetF	  = 0.20
  res@lbLabelFont   	  = "Helvetica-Bold"     ; label font
  res@lbTitleString       = "SSH response (~F10~e~S~~F21~-3~N~m)"
  ;res@lbLabelPosition	  = "Top"
  res@lbTitleAngleF       = 270.
  res@lbTitlePosition      = "Right"
  res@lbTitleDirection     = "Across"
  
  
  res@lbLabelStride        = 10 
  
  res@mpFillOn             = True        ; turn off map fill
  res@mpMinLatF            = latS         ; zoom in on map
  res@mpMaxLatF            = latN
  res@mpMinLonF            = lonL
  res@mpMaxLonF            = lonR
  res@pmTickMarkDisplayMode = "Never"
  ;res@tmXBLabelsOn         = False      ;no lon ticks
  ;res@tmYLLabelsOn         = False      ;no lon ticks
  res@tmYLMode   = "Explicit"
  res@tmYLValues = ispan(-15,15,15)
  res@tmYLLabels = (/"15~S~o~N~S", "0~S~o~N","15~S~o~N~N"/)
  res@tmYLMinorValues  = ispan(-15,15,5) 
  res@tmYLLabelFontHeightF = 0.017

  res@tmXBMode   = "Explicit"
  res@tmXBValues = ispan(-60,20,20)
  res@tmXBLabels = (/"60~S~o~N~W","40~S~o~N~W","20~S~o~N~W","0~S~o~N~","20~S~o~N~E"/)
  ;res@tmXBMinorValues  = ispan(-50,20,5)
  res@tmXBLabelFontHeightF = 0.017
  ;res@mpLandFillColor	   =  "grey"
  ;res@tmXTOn = True ; Turn off tickmarks and labels
  ;res@tmYROn = False 	
  ;res@tmXBOn = True ; Turn off tickmarks and labels
  ;res@tmYLOn = False 	

      
 
  res@mpGeophysicalLineColor = "grey" 
  res@mpFillDrawOrder      = "postDraw"

  res@vpHeightF 	  = 0.80        ; Changes the aspect ratio
  res@vpWidthF 		  = 0.80
  res@gsnCenterStringFontHeightF   = 0.030
  

  

  res@gsnCenterString  = ""   
  res@gsnLeftString  = ""  
  res@gsnRightString  = ""  
    


  ;res@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ;res@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  ;res@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ;res@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
 

  ;res@tmXTBorderOn    = False 
  ;res@tmXBBorderOn    = False 
  ;res@tmXTOn	      = False 

  ;res@tmXTMinorLineColor = "transparent"    
  ;res@tmXBMinorLineColor = "transparent" 
  ;res@tmXTMajorLineColor = "transparent"    
  ;res@tmXBMajorLineColor = "transparent"    
  

  ;res@tmYLBorderOn    = False 
  ;res@tmYLOn	      = False    
  ;res@tmYRBorderOn    = False 
  ;res@tmYROn	      = False 

  res@vcRefAnnoOrthogonalPosF = -1.0             ; move ref vector up
  res@vcRefMagnitudeF         = 0.01             ; define vector ref mag
  res@vcRefAnnoString1          = "0.01 N m~S~-2~N~ "
  res@vcRefLengthF            = 0.045            ; define length of vec ref
 ; vcres@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
  res@vcMinDistanceF          = 0.045
  res@vcLineArrowThicknessF   = 1.
  res@vcLineArrowColor        = "blue"           ; change vector color
  res@gsnScalarContour	      = True
  res@vcLineArrowHeadMaxSizeF = 0.05
  res@vcLineArrowHeadMinSizeF = 0.01
 
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
  res2                      =   True
  res2@gsnDraw              = False        ; don't draw yet
  res2@gsnFrame             = False        ; don't advance frame yet
  res2@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  res2@cnLevels		   = (/0./)   ; set levels
  ;res2@cnLevelSelectionMode = "ManualLevels"   ; set explicit contour levels
  ;res2@cnMinLevelValF       =  -0.12            ; set min contour level
  ;res2@cnMaxLevelValF       =  0.12             ; set max contour level
  ;res2@cnLevelSpacingF      =  0.01            ; set contour spacing   
  res2@cnMonoLineThickness  = True
  res2@cnLineThicknessF     = 4.0
  res2@cnMonoLineColor	    = True
  ;res2@cnLineColor	    = "purple"
  ;res2@cnLineColor	    = "dodgerblue2"
  res2@cnLineColor	    = "grey25"
  ;res2@cnLineLabelFontColor = res2@cnLineColor
  res2@gsnAddCyclic         = False   
  res2@cnLinesOn            = True       ; True is default
  res2@cnLineLabelsOn       = False       ; True is default
  res2@cnInfoLabelOn        = False
  res2@gsnRightString   = ""
  res2@gsnLeftString   = ""
 ; res2@cnFillScaleF        = 2.         ; add extra density
  res2@cnLineLabelDensityF =  2.0
  res2@cnLabelMasking	   = True
  res2@gsnContourZeroLineThicknessF = 2.
 ; res2@gsnContourLineThicknessesScale = 1.
  ;res2@gsnContourNegLineDashPattern = 16
  res2@cnLineDashPattern	=	16
  
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
 ;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
  res3                      =   True
  res3@gsnDraw              = False        ; don't draw yet
  res3@gsnFrame             = False        ; don't advance frame yet
  ;res3@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  ;res3@cnLevels		   = (/10./)   ; set levels
  res3@cnMinLevelValF       = -5.            ; set min contour level
  res3@cnMaxLevelValF       =  5.             ; set max contour level
  res3@cnLevelSpacingF      =  0.5            ; set contour spacing   
  res3@cnMonoLineThickness  = True
  res3@cnLineThicknessF     = 1.5
  res3@cnMonoLineColor	    = True
  ;res3@cnLineColor	    = "purple"
  ;res3@cnLineColor	    = "dodgerblue2"
  res3@cnLineColor	    = "white"
 ; res3@cnLineLabelFontColor = res2@cnLineColor
  res3@gsnAddCyclic         = False   
  res3@cnLinesOn            = True       ; True is default
  res3@cnLineLabelsOn       = True       ; True is default
  res3@cnInfoLabelOn        = False
  res3@gsnRightString   = ""
  res3@gsnLeftString   = ""
 ;res3@cnFillScaleF        = 2.         ; add extra density
  res3@cnLineLabelDensityF =  1.0
  res3@cnLabelMasking	   = True
  res3@gsnContourZeroLineThicknessF = 0.
 ;res3@gsnContourLineThicknessesScale = 1.
  res3@gsnContourNegLineDashPattern = 16
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
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
  res1@cnMaxLevelValF       = res@cnMaxLevelValF             ; set max contour level
  res1@cnLevelSpacingF      = res@cnLevelSpacingF           ; set contour spacing           ; set contour spacing  
  res1@gsnAddCyclic         = False   
  res1@cnLinesOn            = False       ; True is default
  res1@cnLineLabelsOn       = False        ; True is default
  res1@cnInfoLabelOn        = False
  res1@gsnRightString       = ""
  res1@gsnLeftString        = ""
  res1@cnFillScaleF         =  1.5         ; add extra density
  ;res1@cnLineLabelDensityF =  -10.
  res1@cnFillDotSizeF       = 0.0001            ; increase dot size of contours
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
  opt = True
  opt@gsnShadeFillType = "pattern"        ; color is the default
 ; opt@gsnShadeLow = 17                     ; use fill pattern #4
  opt@gsnShadeHigh = 17                  ; use pattern #2   
  ;opt@gsnShadeMid = 2   
;-----------------------------------------------------------------------------
;=============================================================================
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;res@gsnCenterString  = "Thermocline feedback"  
;plotA0 = gsn_csm_contour_map_ce(wks,d_SSH,res)
plotA0 = gsn_csm_vector_scalar_map_ce(wks,ens_d_USTR(month|5,lat|:,lon|:),ens_d_VSTR(month|5,lat|:,lon|:),d_SSH,res)  ; create plot
plotB0 = gsn_csm_contour(wks,d_SSHd,res1)
plotB0 = gsn_contour_shade(plotB0,p95,d_SSHd@_FillValue,opt)
;plotC0 = gsn_csm_contour(wks,d_OISST,res2)
plotC0 = gsn_csm_contour(wks,ens_meanV(month|5,lat|:,lon|:),res2)
plotD0 = gsn_csm_contour(wks,d_GPCP(month|5,lat|:,lon|:),res3)
overlay(plotA0,plotB0)
;overlay(plotA0,plotC0)
overlay(plotA0,plotD0)
plot(4)  =   plotA0 
;

; res3@cnLevels		   = (/8./)   ; set levels
res@cnMinLevelValF       = -0.5            ; set min contour level
res@cnMaxLevelValF       =  0.5             ; set max contour level
res@cnLevelSpacingF      =  0.05            ; set contour spacing  
res1@cnMinLevelValF       = res@cnMinLevelValF            ; set min contour level
res1@cnMaxLevelValF       = res@cnMaxLevelValF             ; set max contour level
res1@cnLevelSpacingF      = res@cnLevelSpacingF           ; set contour spacing           ; set contour spacing  
res@lbTitleString       = "~F33~D~F10~Q~F21~ response (~F10~e~S~~F21~-6~N~ K s~S~-1~N~)"
;res@gsnCenterString  = "Diabatic heating"  
;plot(4) = gsn_csm_contour_map_ce(wks,d_QTd,res)
;plotA1 = gsn_csm_contour_map_ce(wks,d_QT,res)
plotA1 = gsn_csm_vector_scalar_map_ce(wks,ens_d_USTR(month|4,lat|:,lon|:),ens_d_VSTR(month|4,lat|:,lon|:),d_QT,res)  ; create plot
plotB1 = gsn_csm_contour(wks,d_QTd,res1)
plotB1 = gsn_contour_shade(plotB1,d_QTd@_FillValue,p95,opt)
;plotC1 = gsn_csm_contour(wks,d_OISST3,res2) 
plotC1 = gsn_csm_contour(wks,ens_meanV(month|4,lat|:,lon|:),res2)
plotD1 = gsn_csm_contour(wks,d_GPCP(month|4,lat|:,lon|:),res3)
overlay(plotA1,plotB1)
;overlay(plotA1,plotC1)
overlay(plotA1,plotD1)
plot(5)  =   plotA1
;

delete(res2@cnLevels)
delete(res2@gsnContourNegLineDashPattern)
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;===============================================================
;=============================================================== 
;***********************************************
   r2               = True                        ; plot mods desired
   r2@gsnDraw       = False                       ; do not draw
   r2@gsnFrame      = False                       ; do not advance frame
   r2@tiMainString  = ""                       ; title
   r2@tiXAxisFontHeightF = 0.017
   r2@tiYAxisFontHeightF = 0.017
   r2@tiYAxisString = "Correlation, r"                  ; yaxis
   r2@tiXAxisString = "Calendar month"                  ; yaxis
  ; r2@tiYAxisString = "SST (~S~o~N~C)"                  ; yaxis
   ;r2@tiYAxisAngleF  =  270.
   r2@tmXBLabelFontHeightF = 0.015
   r2@tmYLLabelFontHeightF = 0.015
   r2@trYLog              = False                ; log scaling
   r2@trXLog              = False               ; log scaling
   ;r2@trXMinF             = 1940.                 ; manually set lower limit
   ;r2@trXMaxF             = 2020.                ;   "          upper
  r2@trYMinF             =    0.                 ; manually set lower limit
  r2@trYMaxF             =   30.           

   r2@trXMinF             =  0.                 ; manually set lower limit
   r2@trXMaxF             =  11.                 ;   "          upper
  ; r2@gsnYRefLineThicknesses = (/1/)
  ; r2@gsnYRefLineDashPatterns = (/0./)
   ;r2@xyLineColors        = (/"springgreen4","springgreen","transparent"/)
   r2@xyLineColors        = xy_colors_thermocline
   r2@tmXTBorderOn            =  True
   r2@tmXTOn                  =  True
      
   r2@tmYLFormat = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
  ; r2@gsnXRefLine           = (/5./) 
  ; r2@gsnYRefLineDashPattern = 16

   r2@tmYLMode        = "Explicit"	
   r2@tmYLValues      = (/"-0.5","0.0","0.5","1.0"/)
   r2@tmYLLabels      =  r2@tmYLValues
   ;r2@tmYLMinorValues = fspan(-1.0,1.0,11)
   ;print(r2@tmYLMinorValues)
   
   r2@vpWidthF                = 0.60
   r2@vpHeightF               = 0.30
   ;r2@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r2@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   ;r2@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r2@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  ; r2@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; r2@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r2@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r2@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
 
   r2@xyLineThicknesses   = (/2.,2.,2.,2.,2.,2./)       ; Define line thicknesses 
   r2@xyDashPatterns      = (/0.,0.,0.,0.,0.,0./)           ; Dash patterns 
   r2@tmYLLabelFontColor = "foreground"   ; and right axes are same.
   r2@tiYAxisFontColor  = "foreground"
   ;r2@tmBorderLineColor = "deepskyblue"
   ;r2@tmYLMajorLineColor = "red" ; YR axis uses YL resources, can
   ;r2@tmYLMinorLineColor = "red" ; change by setting tmYUseLeft = False 
   ;r2@tmBorderLineColor =  "red"
 
   r2@tmXBMode = "Explicit"
   r2@tmXBValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
   ;r2@tmXBLabels = (/"J", "F","M","A","M","J","J","A","S","O","N","D"/)
   r2@tmXBLabels = (/"Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
   ;r2@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
   
   ;r2@tmXBMinorOn          = False                ; Turn off x-axis (bottom) minor ticks
   ;r2@tmXTMinorOn          = False                ; Turn off x-axis (bottom)minor ticks
   ;r2@tmXTOn               = False                ; Turn off x-axis (top) major ticks
   ;r2@tmYRMinorOn          = False                ; Turn off y-axis (right) minor ticks
   ;r2@tmYLMinorOn          = False                ; Turn off y-axis (left) minor ticks
   ;r2@tmYROn               = False                ; Turn off y-axis (right) major ticks

  ;===legends
   ;===legends
   r2@pmLegendDisplayMode    = "Always"            ; turn on legend
   r2@pmLegendSide           = "Top"               ; Change location of 
   r2@pmLegendParallelPosF   = 0.15                 ; move units right
   r2@pmLegendOrthogonalPosF = -0.70                ; move units down
   r2@pmLegendWidthF         = 0.14                ; Change width and
   r2@pmLegendHeightF        = 0.12                ; height of legend.
   r2@lgLabelFontHeightF     = 0.012                 ; change font height
   r2@lgTitleOn              = False               ; turn on legend title
   r2@lgTitleString          = "SST indices"           ; create legend title
   r2@lgTitleFontHeightF     = .003                ; font of legend title
   r2@lgPerimOn              = False               ; turn off box around
   r2@xyExplicitLegendLabels = (/"SODA3","ORAS5","ORAS4","GODAS","GECCO2","AVISO"/)         ; explicit labels 


;==================================================================================
   r0               = True                        ; plot mods desired
   r0@gsnDraw       = False                       ; do not draw
   r0@gsnFrame      = False                       ; do not advance frame
   r0@tiMainString  = ""                       ; title
   r0@tiXAxisString = "Calendar month"                  ; yaxis
   r0@tiXAxisFontHeightF = 0.017
   r0@tiYAxisFontHeightF = 0.017
   r0@tiYAxisString = "SST (K)"                  ; yaxis
  ; r0@tiYAxisString = "SST (~S~o~N~C)"                  ; yaxis
   ;r0@tiYAxisAngleF  =  270.
   r0@tmXBLabelFontHeightF = 0.017
   r0@tmYLLabelFontHeightF = 0.017
   r0@trYLog              = False                ; log scaling
   r0@trXLog              = False               ; log scaling
   ;r0@trXMinF             = 1940.                 ; manually set lower limit
   ;r0@trXMaxF             = 2020.                ;   "          upper
   r0@trYMinF             =  0.0                 ; manually set lower limit
   r0@trYMaxF             =  1.0           

   r0@trXMinF             =  0.0                 ; manually set lower limit
   r0@trXMaxF             =  11.0                 ;   "          upper
  ; r0@gsnYRefLineThicknesses = (/1/)
  ; r0@gsnYRefLineDashPatterns = (/0./)
   ;r0@xyLineColors        = (/"black","springgreen","transparent"/)
   r0@xyLineColors        =  xy_colors_dheating	;diabatic heating
   r0@tmXTBorderOn            =  True
   r0@tmXTOn                  =  True
      
   r0@tmYLFormat = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
  ; r0@gsnYRefLine           = (/0./) 
   
   r0@vpWidthF                = 0.60
   r0@vpHeightF               = 0.30
   ;r0@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r0@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   ;r0@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r0@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  ; r0@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; r0@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  ; r0@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; r0@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
 
   r0@xyLineThicknesses   = (/1.5,1.5,1.5,1.5,1.5,1.5/)       ; Define line thicknesses 
   r0@xyDashPatterns      = (/0.,0.,0.,0.,0.,0./)           ; Dash patterns 
   r0@tmYLLabelFontColor = "black"   ; and right axes are same.
   r0@tiYAxisFontColor  = "black"
   ;r0@tmBorderLineColor = "deepskyblue"
   ;r0@tmYLMajorLineColor = "red" ; YR axis uses YL resources, can
   ;r0@tmYLMinorLineColor = "red" ; change by setting tmYUseLeft = False 
   ;r0@tmBorderLineColor =  "red"
 
   r0@tmXBMode = "Explicit"
   r0@tmXBValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
   ;r0@tmXBLabels = (/"J", "F","M","A","M","J","J","A","S","O","N","D"/)
   ;r0@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
   r0@tmXBLabels = (/"Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
   
   
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
  r0@pmLegendParallelPosF   = 0.85                 ; move units right
  r0@pmLegendOrthogonalPosF = -0.67                ; move units down
  r0@pmLegendWidthF         = 0.14                ; Change width and
  r0@pmLegendHeightF        = 0.12                ; height of legend.
  r0@lgLabelFontHeightF     = 0.012                 ; change font height
  r0@lgTitleOn              = False               ; turn on legend title
  r0@lgTitleString          = "SST indices"           ; create legend title
  r0@lgTitleFontHeightF     = .003                ; font of legend title
  r0@lgPerimOn              = False               ; turn off box around
  r0@xyExplicitLegendLabels = (/"MERRA2","JRA55","ERAI","ERA5","CFSR","SODA3"/)         ; explicit labels
;-----------------------------------------------------------------------------
;===============================================================
;=============================================================== 
;***********************************************
   r1               = True                        ; plot mods desired
   r1@gsnDraw       = False                       ; do not draw
   r1@gsnFrame      = False                       ; do not advance frame
   r1@tiYAxisString = "SST (K)"      
   r1@tiYAxisAngleF  =  270. 
   r1@xyLineColors          = (/"brown","yellow3","deepskyblue","black","darksalmon","black"/)
   r1@xyLineThicknesses   = (/1.5,1.5,1.5,1.5,1.5,1.5/)       ; Define line thicknesses 
   r1@xyDashPatterns      = (/16.,0.,0.,0.,0.,0./)           ; Dash patterns 
   r1@tiXAxisFontHeightF = 0.017
   r1@tiYAxisFontHeightF = 0.017
   r1@tmXBLabelFontHeightF = 0.017
   r1@tmYRLabelFontHeightF = 0.017
     
   r1@tmYUseLeft          =  False
   r1@trYMinF             =  0.0                 ; manually set lower limit
   r1@trYMaxF             =  1.0                 ;   "          upper
   
   r1@trXMinF             =  0.0                 ; manually set lower limit
   r1@trXMaxF             =  11.0                 ;   "          upper
   
     
   r1@vpWidthF                = 0.60
   r1@vpHeightF               = 0.30
   r1@tmXTBorderOn            =  True
   r1@tmXTOn                  =  False
   r1@tmXTMode = "Explicit"
   r1@tmXTValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
   ;r1@tmXTLabels = (/"1", "2","3","4","5","6","7","8","10","10","11","12"/)
   ;r1@tmXBLabels = (/"J", "F","M","A","M","J","J","A","S","O","N","D"/)
   r1@tmXBLabels = (/"Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

  r1@tmYRLabelFontColor = "brown"   ; and right axes are same.
  r1@tiYAxisFontColor  = "brown"
  r1@tmXTLabelsOn        = False    ; have tick mark labels
  r1@tmXTOn              = True      ; have tick marks
  r1@tmYRLabelFontColor = "brown"   ; and right axes are same.
  r1@tiYAxisFontColor  = "brown"
 ; r1@tmBorderLineColor =  "deepskyblue"
  r1@tmYLBorderOn 	 = False
  r1@tmYRBorderOn	 = True
  r1@tmXTBorderOn	 = False
  r1@tmXBBorderOn	 = False
  ;r1@tmBorderLineColor = "deepskyblue"
  ;r1@tmYRMajorLineColor = "deepskyblue" ; YR axis uses YL r1ources, can
  ;r1@tmYRMinorLineColor = "deepskyblue" ; change by setting tmYUseLeft = False 
  ;r1@gsnYRefLine           = avg(RAPID)
  ;r1@gsnYRefLineColor     = "deepskyblue"
  r1@tmYRFormat   	=  "@*+^sg"          ; trim the uneccecary "0" in labels. 
    
     
   ;r1@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r1@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  ; r1@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; r1@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   ;r1@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r1@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   ;r1@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; r1@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
;========================================================
;========================================================================
;========================================================================
;========================================================================
;----------------------------------------------------------------------------------------------------------------
   delete(r2@tmYLValues)
   delete(r2@tmYLLabels)     
   ;r2@tmYLValues      = (/"2","4","6","8","10"/)
   ;r2@tmYLLabels      =  r2@tmYLValues   
   r2@tiYAxisString = "s.d. [SSH (~F10~e~S~~F21~-3~N~m)]"
  ; r2@tiYAxisString = "Explained variance, ~F10~Z~B~~F21~20~N~ (%)"
   ;r2@tiYAxisString = "s.d. (m)"
   r2@tiMainString  = "Thermocline feedback"                       ; title
   plot(0)  =  gsn_csm_xy(wks,xlabel,std_tcline_s(lev|:,month|:),r2)
   r2@tiMainString  = ""                       ; title
   r2@trYMinF             =  0.0                 ; manually set lower limit
   r2@trYMaxF             =  1.0       


   ;delete(r2@tmYLValues)
   ;delete(r2@tmYLLabels)     
   
   r2@xyExplicitLegendLabels = (/"MERRA2","JRA55","ERAI","ERA5","CFSR","SODA3"/)         ; explicit labels
   r2@pmLegendParallelPosF   = 0.88                 ; move units right
   ;r2@tmYLValues      = (/"0.0","0.5","1.0"/)
   r2@tmYLLabels      =  r2@tmYLValues
   ;r2@tiYAxisString = "Q~B~t~N~ (K s~S~-1~N~)"
   r2@tiYAxisString = "s.d. [~F33~D~F10~Q~F21~(~F10~e~S~~F21~-6~N~ K s~S~-1~N~)]"
 
   delete(r2@xyLineColors)
   r2@xyLineColors        = xy_colors_dheating	;diabatic heating
   std_diab_s = std_diab_s*1000000.
   r2@tmYLFormat = "0@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
   ;plot(5)  =  gsn_csm_xy(wks,xlabel,std_diab_s(lev|:,month|:),r2)

    
    r0@tmYLFormat	   = "0@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
    r1@tmYRFormat	   = "0@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
   ; r0@tiYAxisString = "s.d. [~F33~D~F21~Q (e~S~-6~N~ K s~S~-1~N~)]"
    ;r0@tiYAxisString = "s.d. [~F33~D~F10~Q (e~S~~F21~-6~N~ K s~S~-1~N~)]"
    r0@tiYAxisString = "s.d. [~F33~D~F10~Q~F21~(~F10~e~S~~F21~-6~N~ K s~S~-1~N~)]"
    r1@tiYAxisString = "s.d. [Atl3 SST (K)]"
    r0@tiMainString  = "Diabatic heating"                       ; title
    ;plot(1)  =  gsn_csm_x2y2(wks,xlabel,xlabel,std_diab_s(lev|:,month|:),std_sst_s,r0,r1)
    plot(1)  =  gsn_csm_xy(wks,xlabel,std_diab_s(lev|:,month|:),r0)
    r0@tiMainString  = ""                       ; title
    r2@tmYLFormat = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
    r2@pmLegendDisplayMode    = "Never"            ; turn on legend 
    r0@tmYLFormat	   = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
    r1@tmYRFormat	   = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
    r0@pmLegendDisplayMode    = "Never"            ; turn on legend

    ;r0@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
    ;r0@tiXAxisString = "Time lag (month)"                  ; yaxis
    r0@tiXAxisString = "Calendar month"                  ; yaxis
   ; r2@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
    ;r2@tiXAxisString = "Time lag (month)"                  ; yaxis
    r2@tiXAxisString = "Calendar month"                  ; yaxis
;----------------------------------------------------------------------------------------------------------------

   r2@xyLineThicknesses   = (/1.5,1.5,1.5,1.5,1.5,1.5/)       ; Define line thicknesses 
   ;r2@tiYAxisString = "SSH feedback (e~S~-3~N~m)"
   r2@tiYAxisString = "Thermocline feedback (K)"

   ;r2@tiYAxisString = "r (SST~B~lag=0~N~ | Z~B~20~N~)"                  ; yaxis
   r2@trYMinF             =  -0.4                 ; manually set lower limit
   r2@trYMaxF             =   0.9      
   
    r1@trYMinF             =   -0.4                 ; manually set lower limit
    r1@trYMaxF             =   0.9           
   ;delete(r2@tmYLValues)
   ;delete(r2@tmYLLabels)     
   ;r2@tmYLValues      = (/"-3","0","6"/)
   ;r2@tmYLLabels      =  r2@tmYLValues
   ;delete(r2@xyLineColors)
   r2@xyLineColors        = xy_colors_thermocline
   ;plotA2  =  gsn_csm_xy(wks,xlabel,d_tcline(lev|:,month|:),r2)
   plotA2  =  gsn_csm_x2y2(wks,xlabel,xlabel,d_tcline(lev|:,month|:),std_sst_s,r2,r1)


   ;delete(r2@tmYLValues)
   ;delete(r2@tmYLLabels)    
   ;r2@tiYAxisString = "r (SST~B~lag=0~N~ | ~F33~D~F21~Q)"                  ; yaxis
   r2@tiYAxisString = "~F33~D~F10~Q~F21~(e~S~-6~N~ K s~S~-1~N~)"
   r2@xyLineColors        = xy_colors_dheating	;diabatic heating
  ; r_diab = r_diab*1000000.
   ;r2@trYMinF             = -0.1                 ; manually set lower limit
   ;r2@trYMaxF             =  0.5      
   
;   plotA3  =  gsn_csm_xy(wks,xlabel,r_diab(lev|:,month|:),r2)
   
    r0@trYMinF             =  -0.4                ; manually set lower limit
    r0@trYMaxF             =   0.9      
   
    
    ;r0@tiYAxisString = "~F33~D~F21~Q forcing (e~S~-6~N~ K s~S~-1~N~)"
    r0@tiYAxisString = "~F33~D~F10~Q~F21~ feedback (K)"
    ;r1@tiYAxisString = "SST anomaly (K)"   
    plotA3  =  gsn_csm_x2y2(wks,xlabel,xlabel,r_diab(lev|:,month|:),std_sst_s,r0,r1)
    ; plotA3  =  gsn_csm_xy(wks,xlabel,r_diab(lev|:,month|:),r0)

    ;delete(r0@tmXBLabels)
    ;delete(r2@tiXAxisString)
;================================================================================================
   ;delete(r2@tmYLValues)
   ;delete(r2@tmYLLabels)    
   delete(r2@xyLineColors)
    r2@xyLineColors        = xy_colors_thermocline	
   r2@xyMarkLineMode	 =   "Markers"
   r2@xyMarkers     =  (/16.,16.,16.,16.,16.,16./)
   r2@xyMarkerColors = r2@xyLineColors
   r2@xyMarkerSizeF   = 0.010
   ;r2@trYMinF             = -3.0                 ; manually set lower limit
  ; r2@trYMaxF             =  6.0      
      plotB2  =  gsn_csm_xy(wks,xlabel,d_tclined(lev|:,month|:),r2)


   delete(r2@xyLineColors)
   r2@xyLineColors        = xy_colors_dheating	;diabatic heating
   r2@xyMarkerColors = r2@xyLineColors
   ;r_diabd = r_diabd*1000000.
   r2@trYMinF             = -0.1                 ; manually set lower limit
   r2@trYMaxF             =  0.6      
  ; plotB3  =  gsn_csm_xy(wks,xlabel,r_diabd(lev|:,month|:),r2)
   
   r0@xyMarkLineMode	 =   "Markers"
   r0@xyMarkers     =  (/16.,16.,16.,16.,16.,16./)
   r0@xyMarkerColors = r0@xyLineColors
   r0@xyMarkerSizeF   = 0.010

  ;
   r1@xyMarkLineMode	 =   "Markers"
   r1@xyMarkers     =  (/16.,16.,16.,16.,16.,16./)
   r1@xyMarkerColors = r1@xyLineColors
   r1@xyMarkerSizeF   = 0.010
  
    
  ;plotB3  =  gsn_csm_x2y2(wks,xlabel,xlabel,r_diabd(lev|:,month|:),d_sstd,r0,r1)
   plotB3  =  gsn_csm_xy(wks,xlabel,r_diabd(lev|:,month|:),r0)

  
  overlay(plotA2,plotB2)
  plot(2) = plotA2

;
  overlay(plotA3,plotB3)
  plot(3) = plotA3
;========================================================================

;=============================================================================
;-----------------------------------------------------------------------------
;---------------------------------------------------------------------
;var_names = (/"Developing Phase","Mature Phase","Dissipating Phase"/)
var_names = (/"Thermocline depth","Diabatic heating"/)
printVarSummary(var_names)
print(var_names)

  txres               = True         ; Text resources
  txres@txFontHeightF = 0.012
  txres@txJust        = "CenterLeft"
  txres@txAngleF      =  0.
  txres@txFont                = "Helvetica"

  labels0  = var_names
  xpos0    = fspan(0.25, 0.55,dimsizes(var_names))
; xpos0    = fspan(0.05, 0.70,dimsizes(var_names))

 do i = 0,dimsizes(labels0)-1,1
 ;  gsn_text_ndc(wks,labels0(i),xpos0(i),0.92,txres)
  end do
;delete(i)
;delete(txres)
;=====================================================================
;=====================================================================
;-----------------------------------------------------------------------------
  txres5               = True
  txres5@txAngleF      = 0.          ; Rotate string clockwise 
  txres5@txFontHeightF = 0.035   ; Use same font height as left axis
  txres5@txFont                = "Helvetica-Bold"
 ; txid = gsn_create_text(wks, "SST anomaly (~S~o~N~C)", txres)
   txid5 = gsn_create_text(wks, "e", txres5);
; Move text string to center/right edge of plot.;
  amres5                  = True
  amres5@amParallelPosF   =  -0.65 ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.
  amres5@amOrthogonalPosF =  -0.52    ; This is the center of the plot.
  amres5@amJust           = "TopLeft"
  annoid5i = gsn_add_annotation(plot(4),txid5,amres5)  ; Attach string to plot
 txid5b = gsn_create_text(wks, "f", txres5);
 annoid5ib = gsn_add_annotation(plot(5),txid5b,amres5)  ; Attach string to plot

 txres5@txFontHeightF = 0.025   ; Use same font height as left axis
 amres5@amParallelPosF   =  -0.67 ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.

 txid5c = gsn_create_text(wks, "a", txres5);

 txid5d = gsn_create_text(wks, "b", txres5);
 txid5e = gsn_create_text(wks, "c", txres5);
 txid5f = gsn_create_text(wks, "d", txres5);

 
 annoid5ic = gsn_add_annotation(plot(0),txid5c,amres5)  ; Attach string to plot
 annoid5id = gsn_add_annotation(plot(1),txid5d,amres5)  ; Attach string to plot
 annoid5ie = gsn_add_annotation(plot(2),txid5e,amres5)  ; Attach string to plot
 annoid5if = gsn_add_annotation(plot(3),txid5f,amres5)  ; Attach string to plot
;======================================================================== 
;========================================================================
  latEQ = (/ 3.,  3.,  -3.,   -3., 3./)
  lonEQ = (/-20., 0., 0., -20.,-20./)
;---Resources for the polyline
  pres0                  = True           ; polyline mods desired
  pres0@gsLineThicknessF = 2.            ; line thickness   
  pres0@gsLineColor      = "black"         ; color of lines
  pres0@gsLineDashPattern = 0
  pres0@tfPolyDrawOrder   = "postDraw"
 
 dum0 = gsn_add_polyline(wks,plot(4),lonEQ,latEQ,pres0)
 dum1 = gsn_add_polyline(wks,plot(5),lonEQ,latEQ,pres0)

  latEQ = (/ 3.,  3.,  -3.,   -3., 3./)
  lonEQ = (/-40., 5., 5., -40.,-40./)

  pres0@gsLineThicknessF = 2.5            ; line thickness   
 ; pres0@gsLineColor      = "white"         ; color of lines
  pres0@gsLineDashPattern = 2
 dum0x = gsn_add_polyline(wks,plot(4),lonEQ,latEQ,pres0)
 dum1x = gsn_add_polyline(wks,plot(5),lonEQ,latEQ,pres0)


;---------------------------------------------------------------

 latEQ = (/ 15.,  13.,  -15.,   -15., 15./)
 lonEQ = (/-28., -28., -28., -28.,-28./)
;---Resources for the polyline
 pres0@gsLineThicknessF = 3.            ; line thickness   
 pres0@gsLineDashPattern = 16
 pres0@gsLineColor      = "grey"         ; color of lines
 ;dum2 = gsn_add_polyline(wks,plot(4),lonEQ,latEQ,pres0)
; dum3 = gsn_add_polyline(wks,plot(5),lonEQ,latEQ,pres0)

 
;===============================================================
;===============================================================
  
 
 
  x_0 = 5.
  x_1 = 5.
  y_0 = -10. 
  y_1 = 40.  

  ;y_0 = r2@trYMinF 
  ;y_1 = r2@trYMaxF  

;---Resources for the polyline
  pres1                  = True           ; polyline mods desired
  pres1@gsLineThicknessF = 0.5            ; line thickness   
  pres1@gsLineColor      = "black"         ; color of lines
  pres1@gsLineDashPattern = 0
  pres1@tfPolyDrawOrder   = "postDraw"
  pres1@gsLineColor      = "black"         ; color of lines
  dumL0 = gsn_add_polyline(wks,plot(0), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
  dumL1 = gsn_add_polyline(wks,plot(1), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
  dumL2 = gsn_add_polyline(wks,plot(2), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
  dumL3 = gsn_add_polyline(wks,plot(3), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
;
  x0 = 0.
  x1 = 11.
  y0 = 0.
  dumL0x = gsn_add_polyline(wks,plot(2), (/x0,x1/),(/y0,y0/) ,pres1)
  dumL1x = gsn_add_polyline(wks,plot(3), (/x0,x1/),(/y0,y0/) ,pres1)
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
   resP@gsnPanelTop                      = 0.49
   resP@gsnPanelBottom                   = 0.20
   resP@gsnPanelLeft                     = 0.10
   resP@gsnPanelRight                    = 0.90
   ;resP@gsnPanelTop                      = 0.99
   ;resP@gsnPanelBottom                   = 0.60
   ;resP@gsnPanelLeft                     = 0.10
   ;resP@gsnPanelRight                    = 0.90
   resP@gsnPanelFigureStringsPerimOn 	  = False 
   resP@gsnPaperOrientation 		  = "portrait"
   resP@gsnPanelDebug      		  = True
   resP@amJust   			  = "BottomLeft"
   resP@gsnPanelFigureStringsPerimOn     = False
   resP@gsnPanelFigureStringsFontHeightF = 0.010
   resP@gsnPanelFigureStringsBackgroundFillColor = "transparent"
   resP@gsnPanelFigureStrings= (/"lag = 0","lag = -1","c","d","e","f","g","h","i","j","k","l"/) ; add strings to panel
   resP@gsnMaximize = True                        ; blow up plot
  gsn_panel(wks,plot(4:),(/1,2/),resP)               ; create panel plots
;*******************************************
;===============================================================
 
;===============================================================
;===============================================================
;***********************************************
; panel plot only resources
  res_P                     		 = True            ; modify the panel plot
  res_P@gsnDraw            		 = True       ; don't draw yet
  res_P@gsnFrame            		 = False        ; don't advance frame yet
  res_P@gsnPanelYWhiteSpacePercent 	= 0.0    ; default is 1.0
  res_P@gsnPanelXWhiteSpacePercent 	= 9.0    ; default is 1.0
  res_P@gsnPanelTop                      = 0.96
  res_P@gsnPanelBottom                   = 0.35
  res_P@gsnPanelLeft                     = 0.08
  res_P@gsnPanelRight                    = 0.88
  ;res_P@gsnPanelTop                      = 0.76
  ;res_P@gsnPanelBottom                   = 0.25
  ;res_P@gsnPanelLeft                     = 0.08
  ;res_P@gsnPanelRight                    = 0.88
  ;res_P@gsnPaperOrientation = "portrait"   ; force portrait
 ; res_P@gsnPanelFigureStrings   = (/"(b)","~F0~Peak~F34~;~F0~12 yr~S~-1~N~"/)
  res_P@amJust   			  = "TopLeft"
  res_P@gsnPanelFigureStringsFontHeightF = 0.020
  res_P@gsnPanelFigureStringsPerimOn 	  = True 
  res_P@gsnPanelFigureStringsBackgroundFillColor = "grey"
  ;res_P@gsnPanelFigureStringsJust    =     "CenterRight"
  ;res_P@amJust               =   "CenterRight"
  res_P@gsnMaximize          = True            ; large format
  gsn_panel(wks,plot(0:3),(/2,2/),res_P)     ; now draw as one plot
;======================================================
;======================================================
;===============================================================
;===============================================================
;***********************************************
;==============================================================
frame(wks)
end
