;*************************************************
;*************************************************
begin
;==================================================
latS   = -15.
latN   =  15. 
lonL   = 135.
lonR   = 300.


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



xy_colors_dheating        = (/"red","greenyellow","deepskyblue","black","darksalmon","black"/)
xy_colors_thermocline        = (/"blue","yellow3","red","springgreen4","springgreen","black"/)

;=================================================
;=================================================



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
 
 ;  oisst = lonFlip(oisst)

     
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
;aOISST_s = wgt_areaave_Wrap(oisst_s(time|:,month|:, {lon|210:270}, {lat | -5:5}),1.0, 1.0, 0)	;June 
;printVarSummary(aOISST_s)

;std_aOISST_s = dim_stddev_n_Wrap(aOISST_s,0)
;print(std_aOISST_s)


aOISST = wgt_areaave_Wrap(oisst_s(time|:,month|11, {lon|210:270}, {lat | -5:5}),1.0, 1.0, 0)	;June 
printVarSummary(aOISST)
aOISST = dim_standardize_n(aOISST,0,0)      ; normalize 
;
;aOISST3 = wgt_areaave_Wrap(oisst_s(time|:,month|3, {lon|210:270}, {lat | -5:5}),1.0, 1.0, 0)	;June 
;printVarSummary(aOISST3)
;aOISST3 = dim_standardize_n(aOISST3,0,0)      ; normalize 
;
;print(xxxxxxxxxxxxxxxxxxxxxxxxxx)
;==========================================================================================
;==========================================================================================
d_OISST = regCoef_n(aOISST,oisst_s(:,mm,:,:),0,0)
copy_VarCoords(oisst_s(0,0,:,:),d_OISST)
printVarSummary(d_OISST)
;
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
;rain = lonFlip(rain)
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

  ; aviso = lonFlip(aviso)


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
d_AVISO = regCoef_n(aOISST(11:),aviso_s(:,11,:,:),0,0)
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


;ustr = lonFlip(ustr)
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


;vstr = lonFlip(vstr)
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

  ; ssh = lonFlip(ssh)
   
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
d_SSH = regCoef_n(aOISST,ssh_s(:,11,:,:,:),0,0)
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

   ;cfsr = lonFlip(cfsr)


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
d_CFSR = regCoef_n(aOISST(0:28),cfsr_s(:,10,:,:),0,0)
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

  ; qt = lonFlip(qt)
   
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
d_QT = regCoef_n(aOISST,qt_s(:,10,:,:,:),0,0)
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
 
 lag_ts_map   = esacr(qt_s(lev|:,lat|:,lon|:,month|10,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
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
  wks = gsn_open_wks("eps","FIG_S3")         ; send graphics to PNG file
  plot = new(2,graphic)                ; create graphic array ; only needed if paneling
  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res@cnMinLevelValF       = -100.            ; set min contour level
  res@cnMaxLevelValF       =  100.             ; set max contour level
  res@cnLevelSpacingF      =  10.            ; set contour spacing   
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
  res@lbTitleFontHeightF  = 0.010 
  res@lbLabelFontHeightF  = 0.010
  res@lbLeftMarginF       = 0.10
  res@lbRightMarginF      = 0.55
  ;res@lbLabelOffsetF	  = 0.20
  res@lbLabelFont   	  = "Helvetica-Bold"     ; label font
 ; res@lbTitleString       = "SSH response (e~S~-3~N~m)"
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
  res@pmTickMarkDisplayMode = "Always"
  ;res@tmXBLabelsOn         = False      ;no lon ticks
  ;res@tmYLLabelsOn         = False      ;no lon ticks
  ;res@tmYLMode   = "Explicit"
  ;res@tmYLValues = ispan(-15,15,15)
  ;res@tmYLLabels = (/"15~S~o~N~S", "0~S~o~N","15~S~o~N~N"/)
  ;res@tmYLMinorValues  = ispan(-15,15,5) 
  res@tmYLLabelFontHeightF = 0.012

  ;res@tmXBMode   = "Explicit"
  ;res@tmXBValues = ispan(-60,20,20)
  ;res@tmXBLabels = (/"60~S~o~N~W","40~S~o~N~W","20~S~o~N~W","0~S~o~N~","20~S~o~N~E"/)
  ;res@tmXBMinorValues  = ispan(-50,20,5)
  res@tmXBLabelFontHeightF = 0.012
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
  res@vcRefMagnitudeF         = 0.02             ; define vector ref mag
  res@vcRefAnnoString1          = "0.02 Nm~S~-2~N~ "
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
plotA0 = gsn_csm_vector_scalar_map_ce(wks,ens_d_USTR(month|11,lat|:,lon|:),ens_d_VSTR(month|11,lat|:,lon|:),d_SSH,res)  ; create plot
plotB0 = gsn_csm_contour(wks,d_SSHd,res1)
plotB0 = gsn_contour_shade(plotB0,p95,d_SSHd@_FillValue,opt)
;plotC0 = gsn_csm_contour(wks,d_OISST,res2)
plotC0 = gsn_csm_contour(wks,ens_meanV(month|11,lat|:,lon|:),res2)
plotD0 = gsn_csm_contour(wks,d_GPCP(month|11,lat|:,lon|:),res3)
overlay(plotA0,plotB0)
;overlay(plotA0,plotC0)
overlay(plotA0,plotD0)
plot(0)  =   plotA0 
;

;delete(res1@cnMinLevelValF)       
;delete(res1@cnMaxLevelValF)       
;delete(res1@cnLevelSpacingF)     

res@cnMinLevelValF       = -0.5            ; set min contour level
res@cnMaxLevelValF       =  0.5             ; set max contour level
res@cnLevelSpacingF      =  0.05            ; set contour spacing  
res1@cnMinLevelValF       = res@cnMinLevelValF            ; set min contour level
res1@cnMaxLevelValF       = res@cnMaxLevelValF             ; set max contour level
res1@cnLevelSpacingF      = res@cnLevelSpacingF           ; set contour spacing           ; set contour spacing  
;res@lbTitleString       = "~F33~D~F21~Q response (e~S~-6~N~ K s~S~-1~N~)"
res@lbTitleString       = "~F33~D~F10~Q~F21~ response (~F10~e~S~~F21~-6~N~ K s~S~-1~N~)"
;res@gsnCenterString  = "Diabatic heating"  
;plot(4) = gsn_csm_contour_map_ce(wks,d_QTd,res)
;plotA1 = gsn_csm_contour_map_ce(wks,d_QT,res)
plotA1 = gsn_csm_vector_scalar_map_ce(wks,ens_d_USTR(month|10,lat|:,lon|:),ens_d_VSTR(month|10,lat|:,lon|:),d_QT,res)  ; create plot
plotB1 = gsn_csm_contour(wks,d_QTd,res1)
plotB1 = gsn_contour_shade(plotB1,d_QTd@_FillValue,p95,opt)
;plotC1 = gsn_csm_contour(wks,d_OISST3,res2) 
plotC1 = gsn_csm_contour(wks,ens_meanV(month|10,lat|:,lon|:),res2)
plotD1 = gsn_csm_contour(wks,d_GPCP(month|10,lat|:,lon|:),res3)
overlay(plotA1,plotB1)
;overlay(plotA1,plotC1)
overlay(plotA1,plotD1)
plot(1)  =   plotA1
;

delete(res2@cnLevels)
delete(res2@gsnContourNegLineDashPattern)
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;===============================================================
;=============================================================== 
;***********************************************
 ;=====================================================================
;=====================================================================
;-----------------------------------------------------------------------------
  txres5               = True
  txres5@txAngleF      = 0.          ; Rotate string clockwise 
  txres5@txFontHeightF = 0.020   ; Use same font height as left axis
  txres5@txFont                = "Helvetica-Bold"
 ; txid = gsn_create_text(wks, "SST anomaly (~S~o~N~C)", txres)
   txid5 = gsn_create_text(wks, "a", txres5);
; Move text string to center/right edge of plot.;
  amres5                  = True
  amres5@amParallelPosF   =  -0.60 ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.
  amres5@amOrthogonalPosF =  -0.52    ; This is the center of the plot.
  amres5@amJust           = "TopLeft"
  annoid5i = gsn_add_annotation(plot(0),txid5,amres5)  ; Attach string to plot
 txid5b = gsn_create_text(wks, "b", txres5);
 annoid5ib = gsn_add_annotation(plot(1),txid5b,amres5)  ; Attach string to plot

 
;======================================================================== 
;========================================================================
;======================================================== 
;========================================================
txr =  True
txr@txAngleF       =  0.
txr@txFontHeightF  =  0.015
txr@txFontColor  = "black"
txr@txFont    =  "Helvetica"

amres        =   True
amres@amParallelPosF  =   -0.18    ;left/right pos
amres@amOrthogonalPosF  =  0.45    ;up/down pos
amres@amJust            =  "CenterRight"

txid0 = gsn_create_text(wks,"Atmospheric forcing",txr)
txid0x = gsn_create_text(wks,"Atmospheric forcing",txr)
;txid1 = gsn_create_text(wks,"~F0~Peak~F0~=~F0~12-13 yr~S~-1~N~",txr)
;anid0 =  gsn_add_annotation(plot(2),txid0,amres)
;anid0x =  gsn_add_annotation(plot(3),txid0x,amres)
;
;txr@txFontColor  = "deepskyblue"
txid1 = gsn_create_text(wks,"Delayed oceanic response",txr)
txid1x = gsn_create_text(wks,"Delayed oceanic response",txr)
amres@amParallelPosF  =   0.48    ;left/right pos
;anid1 =  gsn_add_annotation(plot(2),txid1,amres)
;anid1x =  gsn_add_annotation(plot(3),txid1x,amres)
;
txr@txFontColor  = "seagreen3"
txid2 = gsn_create_text(wks,"1979/2016 - 1950/2016",txr)
amres@amOrthogonalPosF  =  0.43    ;up/down pos
amres@amParallelPosF  =   -0.09    ;left/right pos
;anid2 =  gsn_add_annotation(plot(1),txid2,amres)
;===============================================================
;==============================================================

;===============================================================
  latEQ = (/ 5.,  5.,  -5.,   -5., 5./)
  lonEQ = (/210., 270., 270., 210.,210./)
;---Resources for the polyline
  pres0                  = True           ; polyline mods desired
  pres0@gsLineThicknessF = 2.            ; line thickness   
  pres0@gsLineColor      = "black"         ; color of lines
  pres0@gsLineDashPattern = 0
  pres0@tfPolyDrawOrder   = "postDraw"
 
 dum0 = gsn_add_polyline(wks,plot(0),lonEQ,latEQ,pres0)
 dum1 = gsn_add_polyline(wks,plot(1),lonEQ,latEQ,pres0)

  latEQ = (/ 3.,  3.,  -3.,   -3., 3./)
  lonEQ = (/-40., 5., 5., -40.,-40./)

  pres0@gsLineThicknessF = 2.5            ; line thickness   
 ; pres0@gsLineColor      = "white"         ; color of lines
  pres0@gsLineDashPattern = 2
 ;dum0x = gsn_add_polyline(wks,plot(4),lonEQ,latEQ,pres0)
 ;dum1x = gsn_add_polyline(wks,plot(5),lonEQ,latEQ,pres0)


;---------------------------------------------------------------

 latEQ = (/ 15.,  13.,  -15.,   -15., 15./)
 lonEQ = (/-28., -28., -28., -28.,-28./)
;---Resources for the polyline
 pres0@gsLineThicknessF = 3.            ; line thickness   
 pres0@gsLineDashPattern = 16
 pres0@gsLineColor      = "grey"         ; color of lines
 dum2 = gsn_add_polyline(wks,plot(0),lonEQ,latEQ,pres0)
 dum3 = gsn_add_polyline(wks,plot(1),lonEQ,latEQ,pres0)
;===============================================================
;===============================================================
;===============================================================
;===============================================================

   resP             			  = True                        ; panel mods desired
   resP@gsnDraw              		  = True       ; don't draw yet
   resP@gsnFrame         		  = False        ; don't advance frame yet
   resP@gsnPanelYWhiteSpacePercent       = 0.0    ; default is 1.0
   resP@gsnPanelXWhiteSpacePercent       = 0.0    ; default is 1.0
  ; resP@gsnPanelTop                      = 0.49
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
   resP@gsnPanelFigureStringsFontHeightF = 0.012
   resP@gsnPanelFigureStringsBackgroundFillColor = "white"
   resP@gsnPanelFigureStrings= (/"lag = 0","lag = -1","c","d","e","f","g","h","i","j","k","l"/) ; add strings to panel
   resP@gsnMaximize = True                        ; blow up plot
  gsn_panel(wks,plot(0:1),(/2,1/),resP)               ; create panel plots
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
 ; gsn_panel(wks,plot(0:3),(/2,2/),res_P)     ; now draw as one plot
;======================================================
;======================================================
;===============================================================
;===============================================================
;***********************************************
;==============================================================
frame(wks)
end
