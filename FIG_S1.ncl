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
xy_colors_thermocline        = (/"red","yellow2","deepskyblue","blue","deeppink4","black"/)

yrStrt = 1982
yrLast = 2010
;=================================================
;=================================================

 f0      = addfile("atl3_sst_diab_cfsr_1982_2010_atl3.nc", "r")
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
 f0      = addfile("atl3_sst_aviso_1993_2017_atl3.nc", "r")
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

   oisst = lonFlip(oisst)


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
aOISST_s = wgt_areaave_Wrap(oisst_s(time|:,month|:, {lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0)	;June 
printVarSummary(aOISST_s)

std_aOISST_s = dim_stddev_n_Wrap(aOISST_s,0)
print(std_aOISST_s)


;===============================================================================
;===============================================================================
;===============================================================================
;yrStrt = 1982
;yrLast = 2015
;
f0      = addfile("atl3_ssh_soda342_oras5_oras4_godas_gecco2_1982_2015_atl3.nc", "r")
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
nyears = dimsizes(tcline&time)/nmos
print(nyears)

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
f0      = addfile("atl3_diab_merra_jra55_erai_era5_1982_2018_atl3.nc", "r")
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
  wks = gsn_open_wks("eps","FIG_S1")         ; send graphics to PNG file
  plot = new(4,graphic)                ; create graphic array ; only needed if paneling
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
   r2@pmLegendOrthogonalPosF = -0.55                ; move units down
   r2@pmLegendWidthF         = 0.12                ; Change width and
   r2@pmLegendHeightF        = 0.10                ; height of legend.
   r2@lgLabelFontHeightF     = 0.010                 ; change font height
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
  r0@pmLegendParallelPosF   = 0.88                 ; move units right
  r0@pmLegendOrthogonalPosF = -0.60                ; move units down
  r0@pmLegendWidthF         = 0.12                ; Change width and
  r0@pmLegendHeightF        = 0.10                ; height of legend.
  r0@lgLabelFontHeightF     = 0.010                 ; change font height
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
   r1@xyLineColors          = (/"brown","greenyellow","deepskyblue","black","darksalmon","black"/)
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
   ;r2@tiYAxisString = "s.d. [SSH (e~S~-3~N~m)]"
   r2@tiYAxisString = "s.d. [SSH (~F10~e~S~~F21~-3~N~m)]"
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
   r2@tiYAxisString = "s.d. [~F33~D~F21~Q (e~S~-6~N~ K s~S~-1~N~)]"
 
   delete(r2@xyLineColors)
   r2@xyLineColors        = xy_colors_dheating	;diabatic heating
   std_diab_s = std_diab_s*1000000.
   r2@tmYLFormat = "0@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
   ;plot(5)  =  gsn_csm_xy(wks,xlabel,std_diab_s(lev|:,month|:),r2)

    
    r0@tmYLFormat	   = "0@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
    r1@tmYRFormat	   = "0@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
    ;r0@tiYAxisString = "s.d. [~F33~D~F21~Q (e~S~-6~N~ K s~S~-1~N~)]"
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
   r2@tiYAxisString = "~F33~D~F21~Q (e~S~-6~N~ K s~S~-1~N~)"
   r2@xyLineColors        = xy_colors_dheating	;diabatic heating
  ; r_diab = r_diab*1000000.
   ;r2@trYMinF             = -0.1                 ; manually set lower limit
   ;r2@trYMaxF             =  0.5      
   
;   plotA3  =  gsn_csm_xy(wks,xlabel,r_diab(lev|:,month|:),r2)
   
    r0@trYMinF             =  -0.4                ; manually set lower limit
    r0@trYMaxF             =   0.9      
   
    
    ;r0@tiYAxisString = "~F33~D~F21~Q forcing (e~S~-6~N~ K s~S~-1~N~)"
    ;r0@tiYAxisString = "~F33~D~F21~Q feedback (K)"
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
; Move text string to center/right edge of plot.;
  amres5                  = True
  amres5@amParallelPosF   =  -0.65 ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.
  amres5@amOrthogonalPosF =  -0.52    ; This is the center of the plot.
  amres5@amJust           = "TopLeft"
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
 ; gsn_panel(wks,plot(4:),(/1,2/),resP)               ; create panel plots
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
