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
;days = days*1000*4000	;heat constants
;print(days)

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

;f0      = addfile("/home/hyacinth/work/datasets/cfsr/itcz_wind_stress_1979_2019.nc", "r")
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
itcz(4,:) = i_cfsr(TIME|:,{LON|-20})
itcz(3,:) = i_era5(TIME|:,{LON|-20})
itcz(2,:) = i_erai(TIME|:,{LON|-20})
itcz(1,:) = i_jra55(TIME|:,{LON|-20})
itcz(0,:) = i_merra(TIME|:,{LON|-20})

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
delete(i_cfsr)
delete(i_era5)
delete(i_erai)
delete(i_jra55)
delete(i_merra)
delete(itcz)
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



;Calculate the 95% confidence intervals
;
mnmx      = new ((/2,nmos/), float)
mnmx(0,:) = dim_avg_n_Wrap(mean_i_itcz_s,1)+(2.*dim_avg_n_Wrap(std_i_itcz_s,1))
mnmx(1,:) = dim_avg_n_Wrap(mean_i_itcz_s,1)-(2.*dim_avg_n_Wrap(std_i_itcz_s,1))
print(mnmx)
;=========================================================================================
;=========================================================================================
;============================================================================================
;=========================================================================================
;=========================================================================================
f0      = addfile("/home/hyacinth/work/geomar/sstP/R2/a_diab_merra_jra55_erai_era5_1982_2018_basin.nc", "r")
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
    work = dtrend_msg_n(x1,work,False,False,0)
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
;diab_s = dtrend_n(diab_s,False,0)
;
std_diab_s = dim_stddev_n_Wrap(diab_s,0)
print(std_diab_s(month|:,lev|1))
printVarSummary(std_diab_s)
;print(xxxxxxxxxxxxx)
;
mean_diab_s = dim_avg_n_Wrap(diab_s,0)
printVarSummary(mean_diab_s)

;=========================================================================================
;=========================================================================================
;=========================================================================================
;=========================================================================================
;=============================================================================================================
diab_s = dim_standardize_n(diab_s,0,0)      ; normalize 
;i_itcz_s = dim_standardize_n(i_itcz_s,0,0)      ; normalize 
r_itcz = regCoef(diab_s(lev|:,month|:,time|:),i_itcz_s(lev|0:3,month|:,time|:))
r_itcz!0 = "lev"
r_itcz!1 = "month"
printVarSummary(r_itcz)
printMinMax(r_itcz,True)

do ii = 0,3,1
print(r_itcz(month|:,lev|ii))
end do
delete(ii)
;print(xxxxxxxxxxxxx)
;==============================================================

;--------------test for significance

 lag_ts = esacr(i_itcz_s(lev|0:3,month|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(diab_s(lev|:,month|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

  tval = onedtond(r_itcz@tval,dimsizes(r_itcz))

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

r_itczd = r_itcz
r_itczd = where(prob.le.p95,r_itczd,r_itczd@_FillValue) ; Significant TEMP
printVarSummary(r_itczd)
printMinMax(r_itczd,True)
;
do ii = 0,3,1
print(r_itczd(month|:,lev|ii))
end do

;
delete(lag_ts)
delete(lag_ts_map)
delete(tval)
delete(df)
delete(prob)
delete(b)
;=========================================================================================
;=========================================================================================
;=========================================================================================
f0      = addfile("/home/hyacinth/work/geomar/sstP/R2/a_sst_diab_cfsr_1982_2010_basin.nc", "r")
TIME   = f0->time    
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
cfsr    = f0->diab(iYYYY,1)
printVarSummary(cfsr) 
printMinMax(cfsr,True)

  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set the time dimension

;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
 cfsr_d = new(dimsizes(cfsr),typeof(cfsr),cfsr@_FillValue)
 x   = cfsr&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =cfsr(tt::nmos)
    work = dtrend_msg_n(x1,work,False,False,0)
    cfsr_d(tt::nmos) = work
 end do

 printVarSummary(cfsr_d)
 cfsr = cfsr_d
;===
delete([/tt,x,x1,cfsr_d,work/])
;=========================================================================================
;=========================================================================================
;
cfsr_s     = onedtond(ndtooned(cfsr),(/nyears,nmos/))
cfsr_s!0 = "time"
cfsr_s!1 = "month"
printVarSummary(cfsr_s)
delete(cfsr)


;
printMinMax(cfsr_s,True)
cfsr_s = cfsr_s/conform(cfsr_s,days,1)
printMinMax(cfsr_s,True)
;print(xxxxxxxxxxxxx)
;
;
mean_cfsr_s = dim_avg_n_Wrap(cfsr_s,0)
printVarSummary(mean_cfsr_s)
;
cfsr_s = dim_standardize_n_Wrap(cfsr_s,0,0)  ; Or, overwrite with deviations

;=========================================================================================
;=========================================================================================
;=========================================================================================
;=========================================================================================
;=============================================================================================================
r_cfsr = regCoef(cfsr_s(month|:,time|:),i_itcz_s(lev|4,month|:,time|0:33))
r_cfsr!0 = "month"
printVarSummary(r_cfsr)
printMinMax(r_cfsr,True)
print(r_cfsr)
;==============================================================

;--------------test for significance

 lag_ts = esacr(i_itcz_s(lev|4,month|:,time|0:33),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(cfsr_s(month|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

  tval = onedtond(r_cfsr@tval,dimsizes(r_cfsr))

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

r_cfsrd = r_cfsr
r_cfsrd = where(prob.le.p95,r_cfsrd,r_cfsrd@_FillValue) ; Significant TEMP
printVarSummary(r_cfsrd)
printMinMax(r_cfsrd,True)

;delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(cfsr_s)
delete(tval)
delete(df)
delete(prob)
delete(b)
;====================
printVarSummary(r_itcz)
x = new((/5,nmos/),float,r_itcz@_FillValue)
x(0:3,:) = r_itcz	
x(4,:) = r_cfsr		;cfsr
printVarSummary(x)
print(x(:,0))
delete(r_itcz)
r_itcz = x
delete(x)
;
;
x = new((/5,nmos/),float,r_itcz@_FillValue)
x(0:3,:) = r_itczd	
x(4,:) = r_cfsrd	;cfsr
printVarSummary(x)
print(x(:,4))
delete(r_itczd)
r_itczd = x
delete(x)
;=========================================================================================
delete(ii)
delete(mean_diab_s)
delete(mean_cfsr_s)
;delete(std_i_itcz_s)
;
















;=========================================================================================
;=========================================================================================
yrStrt = 1982
yrLast = 2015
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

;
;===============================================================================================================
;===============================================================================================================
;============================================================================================

yrStrt = 1982
yrLast = 2015

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
tcline_s = wgt_areaave_Wrap(ustr_s(time|:,month|:,lev|:, {lon|-40:-20}, {lat | -3:3}),1.0, 1.0, 0)	;June 
tcline_s  = tcline_s*1000.
printVarSummary(tcline_s)


;print(xxxxxxxxxxx)
;
;f0      = addfile("u_weq_stress_merra_jra55_erai_era5_1982_2018.nc", "r")
;TIME   = f0->time    
;YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
;iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
;tcline    = f0->ustress(iYYYY,:)
;printVarSummary(tcline) 
;printMinMax(tcline,True)

 ;tcline  = tcline*1000.

;  delete(f0)
;  delete(TIME) 
;  delete(YYYY)
;  delete(iYYYY)            ; Re-set the time dimension

;;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
; tcline_d = new(dimsizes(tcline),typeof(tcline),tcline@_FillValue)
; x   = tcline&time
; x1  = x(::12)
; do tt = 0,nmos-1
;    work =tcline(tt::nmos,:)
;    work = dtrend_msg_n(x1,work,False,False,0)
;    tcline_d(tt::nmos,:) = work
; end do

; printVarSummary(tcline_d)
; tcline = tcline_d

;===
;delete([/tt,x,x1,tcline_d,work/])
;=========================================================================================
;=========================================================================================
;nyears = dimsizes(tcline&time)/nmos
;
;tcline_s     = onedtond(ndtooned(tcline),(/nyears,nmos,5/))
;tcline_s!0 = "time"
;tcline_s!1 = "month"
;tcline_s!2 = "lev"
;printVarSummary(tcline_s)
;delete(tcline)
;
;tcline_s = dtrend_n(tcline_s,False,0)
;
;mean_tcline = dim_avg_n_Wrap(tcline_s,0)
;printVarSummary(mean_tcline)
;=========================================================================================
;=========================================================================================
;
;tcline_s = dim_standardize_n(tcline_s,0,0)      ; normalize 
;=========================================================================================
;=============================================================================================================
sst = new((/nyears,nmos,5/),float,temp@_FillValue)
sst(:,:,0) = temp
sst(:,:,1) = temp
sst(:,:,2) = temp
sst(:,:,3) = temp
sst(:,:,4) = temp
;sst(:,:,5) = temp
sst!2 = "lev"
printVarSummary(sst)
printVarSummary(temp)
atemp = temp(:,5)
atemp!0 = "time"
print(atemp)
delete(temp)
;print(xxxxxxxxxxx)
;=============================================================================================================
;=========================================================================================
;r_tcline = regCoef(sst(lev|:,month|:,time|:),tcline_s(lev|0:,month|:,time|:))
r_tcline = regCoef(atemp,tcline_s(lev|0:,month|:,time|:))
;
r_tcline!0 = "lev"
r_tcline!1 = "month"
printVarSummary(r_tcline)
printMinMax(r_tcline,True)

do ii=0,4,1
print(r_tcline(ii,:))
end do
delete(ii)
;print(xxxxxxx)
;==============================================================
;==============================================================

;--------------test for significance

 lag_ts = esacr(atemp(time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(tcline_s(lev|:,month|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

  tval = onedtond(r_tcline@tval,dimsizes(r_tcline))

       df = floattointeger(dimsizes(tcline_s(:,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(lag)))) 

printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)

 prob@long_name = "probability"
   printVarSummary(prob)

   

;

r_tclined = r_tcline
r_tclined = where(prob.le.p95,r_tclined,r_tclined@_FillValue) ; Significant TEMP
printVarSummary(r_tclined)
printMinMax(r_tclined,True)

;delete([/tval,df,prob,b/])
delete(lag_ts)
delete(lag_ts_map)
delete(tcline_s)
delete(tval)
delete(df)
delete(prob)
delete(b)
delete(sst)
;=========================================================================================
;=========================================================================================
;=========================================================================================
;=========================================================================================
;===============================================================================================================
;================================================================
; ===============================================================
; Open the file: Read only the user specified period
; ==============================================================
;f      = addfile ("/home/hyacinth/work/datasets/rain/GPCP_precip.mon.mean.nc", "r")
f      = addfile ("/home/hyacinth/work/datasets/olr/olr.mon.mean.nc", "r")
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

;rain    = short2flt(f->precip(iYYYY,:,:))
rain    = f->olr(iYYYY,:,:)
rain = lonFlip(rain)
rain = (rain*-1.)
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
arain = wgt_areaave_Wrap(rain(time |:, {lon|-40:5}, {lat | -3:3}),1.0, 1.0, 0)
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
xx = rain_s(time|:,month|:,lat|:,{lon|-20})
;xx = dim_avg_n_Wrap(rain_s(time|:,month|:,lat|:,{lon|-65:20}),3)
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
   
;d_GPCPd = where(abs(prob).lt.p95,d_GPCPd,d_GPCPd@_FillValue) ; Significant TEMP
d_GPCPd = where(abs(prob).lt.p95,prob,d_GPCPd@_FillValue) ; Significant TEMP
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
;===============================================================================
;===============================================================================
;=========================================================================================
;==============================================================
;=========================================================================================
;=========================================================================================
print("==============================================================")
print("==============================================================")
;===============================================================================================================
;======================================
  wks = gsn_open_wks("eps","FIG_S4")         ; send graphics to PNG file
  plot = new(3,graphic)                ; create graphic array ; only needed if paneling
  
  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res@cnMinLevelValF       = -10.0            ; set min contour level
  res@cnMaxLevelValF       =  10.0             ; set max contour level
  res@cnLevelSpacingF      =  1.0            ; set contour spacing   
  res@cnFillOn             = True         ; turn on color fill
  res@cnFillPalette        = "cmp_b2r"  
  ;res@cnFillPalette        = "BlueWhiteOrangeRed"  
  ;res@cnFillPalette        = "ncl_default"  
  ;res@cnFillPalette        = "gui_default" ; set color map
 ; res@cnFillPalette        = "so4_23"

 ; cmap = read_colormap_file("cmp_b2r")
  ;cmap_r = cmap(::-1,:)
  ;res@cnFillPalette 	   = cmap_r ; for the reverse color map 
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@cnInfoLabelOn        = False
  res@lbLabelBarOn         = True       ; turn off individual lb's 
  ;res@pmLabelBarParallelPosF = 0.50
  ;res@pmLabelBarOrthogonalPosF = -0.11
  res@lbOrientation        = "vertical"
  res@lbTitleFontHeightF  = 0.017 
  res@lbLabelFontHeightF  = 0.017
  res@lbLeftMarginF       = 0.10
  res@lbRightMarginF      = 0.55
  ;res@lbLabelOffsetF	  = 0.20
  res@lbLabelFont   	  = "Helvetica-Bold"     ; label font
  ;res@lbTitleString       = "Precipitation (mm day~S~-1~N~)"
  res@lbTitleString       = "OLR (W m~S~-2~N~)"
  ;res@lbLabelPosition	  = "Top"
  res@lbTitleAngleF       = 270.
  res@lbTitlePosition      = "Right"
  res@lbTitleDirection     = "Across"
  
  
  res@lbLabelStride        = 10 
  
  ;res@mpFillOn             = True        ; turn off map fill
  ;res@mpMinLatF            = latS         ; zoom in on map
  ;res@mpMaxLatF            = latN
  ;res@mpMinLonF            = lonL
  ;res@mpMaxLonF            = lonR
  ;res@pmTickMarkDisplayMode = "Never"
  ;res@tmXBLabelsOn         = False      ;no lon ticks
  ;res@tmYLLabelsOn         = False      ;no lon ticks
  res@tmYLMode   = "Explicit"
  res@tmYLValues = ispan(-5,20,5)
  ;res@tmYLLabels = (/"10~S~o~N~S", "0~S~o~N","10~S~o~N~N"/)
  res@tmYLLabels = (/"-5", "0","5","10","15","20"/)
  res@tmYLMinorValues  = ispan(-5,20,1) 
  res@tmYLLabelFontHeightF = 0.017

  res@tmXBMode = "Explicit"
  res@tmXBValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
  ;res@tmXBLabels = (/"J", "F","M","A","M","J","J","A","S","O","N","D"/)
  res@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
  res@tiXAxisString = "Time lag (month)"                  ; yaxis
  res@tmXBLabelFontHeightF = 0.017


  res@tiXAxisFontHeightF = 0.020
  res@tiYAxisFontHeightF = 0.020


   ;res@tmYLMode = "Explicit"
   ;res@tmYLValues = (/"-5", "0","5","10","15","20"/)
   ;res@tmYLLabels = (/"-5", "0","5","10","15","20"/)
   

  ;res@tmXBMode   = "Explicit"
  ;res@tmXBValues = ispan(-60,20,20)
  ;res@tmXBLabels = (/"60~S~o~N~W","40~S~o~N~W","20~S~o~N~W","0~S~o~N~","20~S~o~N~E"/)
  ;res@tmXBMinorValues  = ispan(-50,20,5)
  ;res@tmXBLabelFontHeightF = 0.017
  ;res@mpLandFillColor	   =  "grey"
  ;res@tmXTOn = True ; Turn off tickmarks and labels
  ;res@tmYROn = False 	
  ;res@tmXBOn = True ; Turn off tickmarks and labels
  ;res@tmYLOn = False 	

      
 
  ;res@mpGeophysicalLineColor = "grey" 
  ;res@mpFillDrawOrder      = "postDraw"

  res@vpHeightF 	  = 0.45        ; Changes the aspect ratio
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
  res3                      =   True
  res3@gsnDraw              = False        ; don't draw yet
  res3@gsnFrame             = False        ; don't advance frame yet
  ;res3@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  ;res3@cnLevels		   = (/10./)   ; set levels
  res3@cnMinLevelValF       = -300.            ; set min contour level
  res3@cnMaxLevelValF       =  300.             ; set max contour level
  res3@cnLevelSpacingF      =  10.            ; set contour spacing   
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
  res1                      =   True
  res1@gsnDraw              = False        ; don't draw yet
  res1@gsnFrame             = False        ; don't advance frame yet
  res1@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res1@cnMinLevelValF       = res@cnMinLevelValF            ; set min contour level
  res1@cnMaxLevelValF       = res@cnMaxLevelValF            ; set max contour level
  res1@cnLevelSpacingF      = res@cnLevelSpacingF           ; set contour spacing 
  ;res1@gsnAddCyclic         = False   
  res1@cnLinesOn            = False       ; True is default
  res1@cnLineLabelsOn       = False        ; True is default
  res1@cnInfoLabelOn        = False
  res1@gsnRightString       = ""
  res1@gsnLeftString        = ""
  res1@cnFillScaleF         =  1.5         ; add extra density
  ;res1@cnLineLabelDensityF =  -10.
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
   res@tiYAxisString = "Latitude (~S~o~N~N)"                  ; yaxis
   plotA0 = gsn_csm_contour(wks,d_GPCP({lat|-7:17},month|:),res)            ; default hovmueller
   plotB0 = gsn_csm_contour(wks,d_GPCPd({lat|-7:17},month|:),res1)            ; default hovmueller
   plotB0 = gsn_contour_shade(plotB0,5.,p95,opt)
   ;plotB0 = gsn_contour_shade(plotB0,max(d_GPCP),min(d_GPCP),opt)
   plotC0 = gsn_csm_contour(wks,a_rain_s({lat|-7:17},month|:),res3)            ; default hovmueller

   
   overlay(plotA0,plotC0)
   overlay(plotA0,plotB0)
  ; plot(0)  =   plotA0 
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;===============================================================
;===============================================================




   r2               = True                        ; plot mods desired
   r2@gsnDraw       = False                       ; do not draw
   r2@gsnFrame      = False                       ; do not advance frame
   r2@tiMainString  = ""                       ; title
   r2@tiXAxisFontHeightF = 0.015
   r2@tiYAxisFontHeightF = 0.015
   r2@tiYAxisString = "Mean ITCZ (~S~o~N~N)"                  ; yaxis
   r2@tiXAxisString = "Calendar month"                  ; yaxis
  ; r2@tiYAxisString = "SST (~S~o~N~C)"                  ; yaxis
   ;r2@tiYAxisAngleF  =  270.
   r2@tmXBLabelFontHeightF = 0.017
   r2@tmYLLabelFontHeightF = 0.017
   r2@trYLog              = False                ; log scaling
   r2@trXLog              = False               ; log scaling
   ;r2@trXMinF             = 1940.                 ; manually set lower limit
   ;r2@trXMaxF             = 2020.                ;   "          upper
  r2@trYMinF             =   -10.                 ; manually set lower limit
  r2@trYMaxF             =    10.           

   r2@trXMinF             =  0.                 ; manually set lower limit
   r2@trXMaxF             =  11.                 ;   "          upper
  ; r2@gsnYRefLineThicknesses = (/1/)
  ; r2@gsnYRefLineDashPatterns = (/0./)
   ;r2@xyLineColors        = (/"springgreen4","springgreen","transparent"/)
   r2@xyLineColors        = xy_colors_dheating
   r2@tmXTBorderOn            =  True
   r2@tmXTOn                  =  True
      
   r2@tmYLFormat = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
  ; r2@gsnXRefLine           = (/5./) 
  ; r2@gsnYRefLineDashPattern = 16

   ;r2@tmYLMode        = "Explicit"	
  ; r2@tmYLValues      = (/"-0.5","0.0","0.5","1.0"/)
   ;r2@tmYLLabels      =  r2@tmYLValues
   ;r2@tmYLMinorValues = fspan(-1.0,1.0,11)
   ;print(r2@tmYLMinorValues)
   
   r2@vpWidthF                = 0.80
   r2@vpHeightF               = 0.40
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
   r2@tmXBLabels = (/"J", "F","M","A","M","J","J","A","S","O","N","D"/)
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
   r2@pmLegendOrthogonalPosF = -0.45                ; move units down
   r2@pmLegendWidthF         = 0.12                ; Change width and
   r2@pmLegendHeightF        = 0.10                ; height of legend.
   r2@lgLabelFontHeightF     = 0.010                 ; change font height
   r2@lgTitleOn              = False               ; turn on legend title
   r2@lgTitleString          = "SST indices"           ; create legend title
   r2@lgTitleFontHeightF     = .003                ; font of legend title
   r2@lgPerimOn              = False               ; turn off box around
   r2@xyExplicitLegendLabels = (/"MERRA2","JRA55","ERAI","ERA5","CFSR","SODA331"/)         ; explicit labels


;==================================================================================
   r0               = True                        ; plot mods desired
   r0@gsnDraw       = False                       ; do not draw
   r0@gsnFrame      = False                       ; do not advance frame
   r0@tiMainString  = ""                       ; title
   r0@tiYAxisString = "ITCZ feedback (~S~o~N~N)"                  ; yaxis
   r0@tiXAxisString = "Calendar month"                  ; yaxis
   r0@tiXAxisFontHeightF = 0.020
   r0@tiYAxisFontHeightF = 0.020
   ; r0@tiYAxisString = "SST (~S~o~N~C)"                  ; yaxis
   ;r0@tiYAxisAngleF  =  270.
   r0@tmXBLabelFontHeightF = 0.020
   r0@tmYLLabelFontHeightF = 0.020
   r0@trYLog              = False                ; log scaling
   r0@trXLog              = False               ; log scaling
   ;r0@trXMinF             = 1940.                 ; manually set lower limit
   ;r0@trXMaxF             = 2020.                ;   "          upper
   r0@trYMinF             =   -1.5                 ; manually set lower limit
   r0@trYMaxF             =    1.5           

   r0@trXMinF             =  0.                 ; manually set lower limit
   r0@trXMaxF             =  11.                 ;   "          upper
  ; r0@gsnYRefLineThicknesses = (/1/)
  ; r0@gsnYRefLineDashPatterns = (/0./)
   ;r0@xyLineColors        = (/"black","springgreen","transparent"/)
   r0@xyLineColors        =  xy_colors_dheating
   r0@tmXTBorderOn            =  True
   r0@tmXTOn                  =  True
      
   r0@tmYLFormat = "0@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
  ; r0@gsnYRefLine           = (/0./) 
   
   r0@vpWidthF                = 0.80
   r0@vpHeightF               = 0.45
   ;r0@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r0@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   ;r0@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r0@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  ; r0@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; r0@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  ; r0@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; r0@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
 
   r0@xyLineThicknesses   = (/1.5,5.5,1.5,1.5,1.5,1.5/)       ; Define line thicknesses 
   r0@xyDashPatterns      = (/16.,16.,16.,0.,0.,0./)           ; Dash patterns 
   r0@tmYLLabelFontColor = "black"   ; and right axes are same.
   r0@tiYAxisFontColor  = "black"
   ;r0@tmBorderLineColor = "deepskyblue"
   ;r0@tmYLMajorLineColor = "red" ; YR axis uses YL resources, can
   ;r0@tmYLMinorLineColor = "red" ; change by setting tmYUseLeft = False 
   ;r0@tmBorderLineColor =  "red"
 
   r0@tmXBMode = "Explicit"
   r0@tmXBValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
   ;r0@tmXBLabels = (/"J", "F","M","A","M","J","J","A","S","O","N","D"/)
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
  ;r0@xyExplicitLegendLabels = (/"MERRA2","JRA55","ERAI","ERA5","CFSR","SODA331"/)         ; explicit labels
  ;r0@xyExplicitLegendLabels = (/"-2~F33~s~F21~ ITCZ","Mean ITCZ","ERAI","ERA5","CFSR","SODA331"/)         ; explicit labels
  r0@xyExplicitLegendLabels = (/"Mean ITCZ ~F18~S~F21~2~F33~s~N","Mean ITCZ","ERAI","ERA5","CFSR","SODA331"/)         ; explicit labels
;-----------------------------------------------------------------------------
;===============================================================
;=============================================================== 
;***********************************************
   r1               = True                        ; plot mods desired
   r1@gsnDraw       = False                       ; do not draw
   r1@gsnFrame      = False                       ; do not advance frame
   r1@tiYAxisString = "Mean ITCZ (~S~o~N~N)"                  ; yaxis
   r1@tiXAxisString = ""                  ; yaxis
  ; r1@tiYAxisAngleF  =  270. 
   r1@xyLineColors        =  (/"deeppink4","deeppink4","deeppink4","deeppink4","deeppink4","deeppink4"/)
   r1@xyLineThicknesses   = (/1.5,1.5,1.5,1.5,1.5,1.5/)       ; Define line thicknesses 
   r1@xyDashPatterns      = (/0.,0.,0.,0.,0.,0./)           ; Dash patterns 
   r1@tiXAxisFontHeightF = 0.020
   r1@tiYAxisFontHeightF = 0.020
   r1@tmXBLabelFontHeightF = 0.020
   r1@tmYRLabelFontHeightF = 0.020
     
   r1@tmYUseLeft          =  False
   r1@trYMinF             = -15.                 ; manually set lower limit
   r1@trYMaxF             =  15.                 ;   "          upper
   
   r1@trXMinF             =  0.0                 ; manually set lower limit
   r1@trXMaxF             =  11.0                 ;   "          upper
   
     
   r1@vpWidthF                = 0.80
   r1@vpHeightF               = 0.40
   r1@tmXTBorderOn            =  True
   ;r1@tmXTOn                  =  False
   r1@tmXBMode = "Explicit"
   r1@tmXBValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
   ;r0@tmXBLabels = (/"J", "F","M","A","M","J","J","A","S","O","N","D"/)
   r1@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)

     
 ; r1@tmYRLabelFontColor = "deeppink4"   ; and right axes are same.
  ;r1@tiYAxisFontColor  = "deeppink4"
  r1@tmXTLabelsOn        = False    ; have tick mark labels
  r1@tmXTOn              = True      ; have tick marks
 ; r1@tmYRLabelFontColor = "deeppink4"   ; and right axes are same.
  ;r1@tiYAxisFontColor  = "deeppink4"
 ; r1@tmBorderLineColor =  "deepskyblue"
  ;r1@tmYLBorderOn 	 = False
  r1@tmYRBorderOn	 = True
  ;r1@tmXTBorderOn	 = False
  ;r1@tmXBBorderOn	 = False
  ;r1@tmBorderLineColor = "deepskyblue"
  ;r1@tmYRMajorLineColor = "deepskyblue" ; YR axis uses YL r1ources, can
  ;r1@tmYRMinorLineColor = "deepskyblue" ; change by setting tmYUseLeft = False 
  ;r1@gsnYRefLine           = avg(RAPID)
  ;r1@gsnYRefLineColor     = "deepskyblue"
  r1@tmYRFormat   	=  "0@*+^sg"          ; trim the uneccecary "0" in labels. 
    
     
   ;r1@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r1@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
  ; r1@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; r1@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   ;r1@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   ;r1@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   ;r1@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ; r1@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward

;-----------------------------------------------------------------------------
;=============================================================================
;print(mean_tcline)
 ;plot(1)  =  gsn_csm_xy(wks,xlabel,r_tcline(lev|:,month|:),r2)
 ; plotA0  =  gsn_csm_x2y2(wks,xlabel,xlabel,r_itcz(lev|:,month|:),mean_i_itcz_s(lev|0:3,month|:),r0,r1)

 ; plotD0  =  gsn_csm_x2y2(wks,xlabel,xlabel,r_itcz(lev|:,month|:),mean_i_itcz_s(lev|0:3,month|:),r0,r1)


  ;r0@xyLineColors        =   r1@xyLineColors

   
 ;r1@xyLineColors        =  (/"deeppink4","deeppink4","deeppink4","deeppink4","deeppink4","deeppink4"/)

  r0@xyLineColors        =  (/"black","black","black","deeppink4","deeppink4","deeppink4"/)

  plotD0  =  gsn_csm_xy(wks,xlabel,itcz_2sigma(lev|0:1,time|:),r0)
   r0@pmLegendDisplayMode    = "Never"            ; turn on legend
  plotE0  =  gsn_csm_xy(wks,xlabel,itcz_2sigma(lev|2:2,time|:),r0)
  r0@pmLegendDisplayMode    = "Always"            ; turn on legend
  r0@xyLineColors        =	 (/"deeppink4","deeppink4","deeppink4","deeppink4","deeppink4","deeppink4"/)
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
  rts          			= True
  rts@gsnDraw  			= False       ; don't draw yet
  rts@gsnFrame 			= False       ; don't advance frame yet
  rts@gsnCenterString 		= ""  
  rts@gsnLeftString   		= ""  
  rts@gsnRightString   		= ""  
  rts@xyDashPatterns  		=   (/0,0,0,0,0,0,0/)   ; choose dash patterns
  rts@xyLineThicknesses         =  2.0
  ;rts@gsnYRefLine	  =      (/0.3/)
  rts@xyLineColors        = (/"red","deepskyblue4","cyan","blue","deeppink4","green","cyan","red","cyan"/)
  ;rts@trXMinF 	     		= -60.                   ; min value on y-axis
  ;rts@trXMaxF  	      		=  20.                   ; max value on y-axis
  ;rts@trYMinF 	     		= -.3                   ; min value on y-axis
  ;rts@trYMaxF  	      		=  1.5                   ; max value on y-axis
  rts@tmYLFormat       		= "0@*+^sg" 
  rts@xyLineThicknesses         =  5.0
  rts@xyLineColors        = (/"red","deepskyblue4","cyan","blue","deeppink4","green","cyan","red","cyan"/)
  rts@xyMarkLineMode = "MarkLines"               ; Markers *and* lines
  rts@xyMarkers      = (/4,12,6,7,8,9, 16, 16,16/)               ; 3 different markers
; rts@xyMarkers      = (/4,4,4,4/)               ; 3 different markers
  rts@xyMarkerColors := rts@xyLineColors ; 3 different colors
  rts@xyMarkerSizeF  = 0.008
  rts@tiYAxisString = ""  

  delete(rts@xyMarkerColors)
  delete(rts@xyMarkers)
  delete(rts@xyMarkLineMode)
 ; delete(rts@xyLineThicknesses)

  ; Create a plot with the area between both curves filled in pink.
  rts@gsnXYFillColors = "red"
  rts@xyLineColors     = (/-1,-1,-1,-1,-1,-1,-1,-1,-1/)                         ; We don't want the line, so make it transparent. 
  rts@gsnXYFillOpacities	= (/0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5/)  
 ; plotE0  = gsn_csm_xy (wks,xlabel,mnmx(0:1,:),rts)  ; Create another filled XY plot.


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





  overlay(plotA0,plotE0)
  overlay(plotA0,plotD0)
 
  plot(0)  =   plotA0 

  delete(r0@gsnXYFillColors)
  delete(r0@xyLineColors) 
  delete(r0@gsnXYFillOpacities) 

 ;r1@xyLineColors        =  (/"red","deeppink4","deeppink4","deeppink4","deeppink4","deeppink4"/)
 r0@xyDashPatterns      = (/0.,0.,0.,0.,0.,0./)           ; Dash patterns 
 r0@xyExplicitLegendLabels = (/"MERRA2","JRA55","ERAI","ERA5","CFSR","SODA3"/)         ; explicit labels
 r0@xyLineThicknesses   = (/1.5,1.5,1.5,1.5,1.5,1.5/)       ; Define line thicknesses 
 r0@xyLineColors        =  xy_colors_dheating
;  plotA1  =  gsn_csm_x2y2(wks,xlabel,xlabel,r_tcline(lev|:,month|:),mean_i_itcz_s(lev|:,month|:),r0,r1)

 r2@pmLegendDisplayMode    = "Never"            ; turn on legend
 r2@trYMinF             =   -1.5                 ; manually set lower limit
 r2@trYMaxF             =    1.5           
 r2@tmYLFormat = "0@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
; plotA1  =  gsn_csm_xy(wks,xlabel,r_itcz(lev|:,month|:),r2)


  r0@trYMinF             =    -2.                 ; manually set lower limit
  r0@trYMaxF             =     10.           
  r0@tiYAxisString = "~F33~t~B~c~N~~F21~ (e~S~-3~N~ N m~S~-2~N~)"                  ; yaxis
;
  r1@trYMinF             =   -2.                 ; manually set lower limit
  r1@trYMaxF             =    8.           
  r1@xyLineThicknesses   = (/1.5,1.5,1.5,1.5,1.5,1.5/)       ; Define line thicknesses 
  r1@tiYAxisString = "OLR (W m~S~-2~N~)"                  ; yaxis
  r1@tmYLFormat = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
  r1@tiXAxisString = "Time lag (month)"                  ; yaxis
  ;r1@tiYAxisString = "Precipitation (mm day~S~-1~N~)"                  ; yaxis

  r0@pmLegendOrthogonalPosF = -0.3                ; move units down
  r0@pmLegendParallelPosF   = 0.12                 ; move units right
  r0@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
  r0@tiXAxisString = "Time lag (month)"                  ; yaxis
  r0@tmYLFormat = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
 ;plotA1  =  gsn_csm_x2y2(wks,xlabel,xlabel,r_tcline(lev|:,month|:),d_arain(month|:),r0,r1)

  plotA1  =  gsn_csm_xy(wks,xlabel,d_arain(month|:),r1)


r0@pmLegendDisplayMode    = "Never"            ; turn on legend
;------------------------------------------------------------------------
;----------------------------------------------------------------------------------------------------------------
   delete(r0@xyLineColors)
   r0@xyLineColors        = xy_colors_dheating
   r0@xyMarkLineMode	 =   "Markers"
   r0@xyMarkers     =  (/16.,16.,16.,16.,16.,16./)
   r0@xyMarkerColors = r0@xyLineColors
   r0@xyMarkerSizeF   = 0.010
   plotB0  =  gsn_csm_xy(wks,xlabel,r_itczd(lev|:,month|:),r0)
  
 
   delete(r1@xyLineColors)
   r1@xyLineColors        =  (/"deeppink4","deeppink4","deeppink4","deeppink4","deeppink4","deeppink4"/)
   r1@xyMarkLineMode	 =   "Markers"
   r1@xyMarkers     =  (/16.,16.,16.,16.,16.,16./)
   r1@xyMarkerColors = r1@xyLineColors
   r1@xyMarkerSizeF   = 0.010
  ;plotB1  =  gsn_csm_x2y2(wks,xlabel,xlabel,r_tclined(lev|:,month|:),d_araind(month|:),r0,r1)

   plotB1  =  gsn_csm_xy(wks,xlabel,d_araind(month|:),r1)
;
 ;  overlay(plotA0,plotB0)
  ; plot(0) = plotA0

  overlay(plotA1,plotB1)
   plot(1) = plotA1
print(d_araind)
;========================================================================
;========================================================================
;========================================================================
;========================================================================
;=============================================================================
;-----------------------------------------------------------------------------
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

; do i = 0,dimsizes(labels0)-1,1
 ;  gsn_text_ndc(wks,labels0(i),xpos0(i),0.92,txres)
;  end do
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
   txid5 = gsn_create_text(wks, "a", txres5);
; Move text string to center/right edge of plot.;
  amres5                  = True
  amres5@amParallelPosF   =  -0.62 ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.
  amres5@amOrthogonalPosF =  -0.57    ; This is the center of the plot.
  amres5@amJust           = "TopLeft"
  annoid5i = gsn_add_annotation(plot(0),txid5,amres5)  ; Attach string to plot
 txid5b = gsn_create_text(wks, "b", txres5);
 annoid5ib = gsn_add_annotation(plot(1),txid5b,amres5)  ; Attach string to plot

 txres5@txFontHeightF = 0.035   ; Use same font height as left axis
 amres5@amParallelPosF   =  -0.67 ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.

 txid5c = gsn_create_text(wks, "c", txres5);

 txid5d = gsn_create_text(wks, "d", txres5);
 txid5e = gsn_create_text(wks, "e", txres5);
 txid5f = gsn_create_text(wks, "f", txres5);

 
; annoid5ic = gsn_add_annotation(plot(2),txid5c,amres5)  ; Attach string to plot
 ;annoid5id = gsn_add_annotation(plot(3),txid5d,amres5)  ; Attach string to plot
 ;annoid5ie = gsn_add_annotation(plot(4),txid5e,amres5)  ; Attach string to plot
 ;annoid5if = gsn_add_annotation(plot(5),txid5f,amres5)  ; Attach string to plot

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
;txid1 = gsn_create_text(wks,"~F0~Peak~F0~=~F0~12-13 yr~S~-1~N~",txr)
;anid0 =  gsn_add_annotation(plot(1),txid0,amres)
;
;txr@txFontColor  = "deepskyblue"
txid1 = gsn_create_text(wks,"Delayed oceanic response",txr)
amres@amParallelPosF  =   0.48    ;left/right pos
;anid1 =  gsn_add_annotation(plot(1),txid1,amres)
;
txr@txFontColor  = "seagreen3"
txid2 = gsn_create_text(wks,"1979/2016 - 1950/2016",txr)
amres@amOrthogonalPosF  =  0.43    ;up/down pos
amres@amParallelPosF  =   -0.09    ;left/right pos
;anid2 =  gsn_add_annotation(plot(1),txid2,amres)
;===============================================================
;==============================================================
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
 ; dumL2 = gsn_add_polyline(wks,plot(2), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
  ;dumL3 = gsn_add_polyline(wks,plot(3), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
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
  res_P@gsnPanelTop                      = 0.95
  res_P@gsnPanelBottom                   = 0.15
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
