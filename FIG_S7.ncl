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


;xy_colors_dheating        = (/"red","greenyellow","deepskyblue","black","darksalmon","black"/)
xy_colors_thermocline        = (/"deeppink4","greenyellow","deepskyblue","red","black","springgreen4"/)
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
f0      = addfile("/home/hyacinth/work/datasets/cfsr/itcz_cfsr_wind_stress_ocnf.gdas.1979_2019.grb2.nc", "r")
;f0      = addfile("/home/hyacinth/work/datasets/cfsr/itcz_wind_stress_1979_2019.nc", "r")
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
do ii =0,4,1
print(std_i_itcz_s(:,ii))
end do
delete(ii)
;
mean_i_itcz_s = dim_avg_n_Wrap(i_itcz_s,0)
printVarSummary(mean_i_itcz_s)
;print(xxxxxxxxxxxxx)
;=========================================================================================
;=========================================================================================
;============================================================================================
;=========================================================================================
;=========================================================================================
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
;
std_rain_s = dim_stddev_n_Wrap(arain_s,0)
print(std_rain_s)

;print(xxxxxxxxx)
;
;=========================================================================================
;=========================================================================================
;===============================================================================
;===============================================================================
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

wustr_s = wgt_areaave_Wrap(ustr_s(time |:,month|:,lev|:, {lon|-40:-20}, {lat | -3:3}),1.0, 1.0, 0)
printVarSummary(wustr_s)
;
;
std_wustr_s = dim_stddev_n_Wrap(wustr_s,0)
printVarSummary(std_wustr_s)
std_wustr_s = std_wustr_s*1000.
do ii =0,4,1
print(std_wustr_s(:,ii))
end do
delete(ii)

;print(xxxxxxxxxx)
;
























;=========================================================================================
;=========================================================================================
yrStrt = 1984
yrLast = 2009
;========================================================================================
;========================================================================================
;===============================================================================================================
;============================================================

;===============================================================================================================
;===============================================================================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
;f0      = addfile("/home/hyacinth/work/datasets/tropflux/NETFLX_tropflux_1m_1979_2017.nc", "r")
f0      = addfile("/home/hyacinth/work/datasets/oaflux/qnet_oaflux.nc", "r")
  TIME   = f0->TIME
 ; TIME   = f0->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt.and. YYYY.le.yrLast)
  ;oaflux    = f0->netflux(iYYYY,:,:)
  oaflux    = f0->QNET(iYYYY,:,:)
 
  printVarSummary(oaflux) 
  printMinMax(oaflux,True)

   oaflux = lonFlip(oaflux)

   oaflux!0  =  "time"
  oaflux!1  =  "lat"
  oaflux!2  =  "lon"
   
   delete(f0)
   delete(TIME) 
   delete(YYYY)
   delete(iYYYY)            ; Re-set the time dimension
;=============================================================================================================
;=============================================================================================================
;================================================================
; ===============================================================
; Open the file: Read only the user specified period
; ==============================================================
f      = addfile("/home/hyacinth/work/datasets/oaflux/sw_isccp.nc", "r")
  TIME   = f->TIME
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt.and. YYYY.le.yrLast)

  SW    = f->NSWRS(iYYYY,:,:)
  
   SW!0  =  "time"
   SW!1  =  "lat"
   SW!2  =  "lon"
  
 ;  SW = -SW
   SW = lonFlip(SW)
  
printVarSummary(SW)       ; variable overview 
; ==============================================================
mean_SW = dim_avg_n_Wrap(SW,0)
;===================================
;===============================================================================================================
;===============================================================================================================
;===============================================================================================================
delete([/YYYY,iYYYY,f,TIME/])


;print(xxxxxxxxxxxxxxxxxx)
;================================================================
; ===============================================================
; Open the file: Read only the user specified period
; ==============================================================
  f      = addfile("/home/hyacinth/work/datasets/oaflux/lw_isccp.nc", "r")
  TIME   = f->TIME
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt.and.YYYY.le.yrLast)

  LW    = f->NLWRS(iYYYY,:,:)
 
  
   LW!0  =  "time"
   LW!1  =  "lat"
   LW!2  =  "lon"
  

   LW = -LW
   LW = lonFlip(LW)
    printVarSummary(LW)       ; variable overview
; ==============================================================
mean_LW = dim_avg_n_Wrap(LW,0)
;===================================
;===============================================================================================================
;===============================================================================================================
;===============================================================================================================
delete([/YYYY,iYYYY,f,TIME/])
;print(xxxxxxxxxxxxxxxxxx)
;================================================================
; ===============================================================
; Open the file: Read only the user specified period
; ==============================================================
  f      = addfile("/home/hyacinth/work/datasets/oaflux/lh_oaflux.nc", "r")
  TIME   = f->TIME
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

  LH    = f->LHTFL(iYYYY,:,:)
 
   LH!0  =  "time"
   LH!1  =  "lat"
   LH!2  =  "lon"
  
   LH = lonFlip(LH)
   LH = -LH
    printVarSummary(LH)       ; variable overview
; ==============================================================
mean_LH = dim_avg_n_Wrap(LH,0)
;===================================

;=====================================
;===============================================================================================================
;===============================================================================================================
;===============================================================================================================
delete([/YYYY,iYYYY,f,TIME/])


;print(xxxxxxxxxxxxxxxxxx)
;================================================================
; ===============================================================
; Open the file: Read only the user specified period
; ==============================================================
  f      = addfile("/home/hyacinth/work/datasets/oaflux/sh_oaflux.nc", "r")
  TIME   = f->TIME
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

  SH    = f->SHTFL(iYYYY,:,:)

   SH!0  =  "time"
   SH!1  =  "lat"
   SH!2  =  "lon"

  SH = -SH
  SH = lonFlip(SH)
 printVarSummary(SH)       ; variable overview
; ==============================================================
; ==============================================================
mean_SH = dim_avg_n_Wrap(SH,0)
;--------------------------------------------------------------
delete([/YYYY,iYYYY,f,TIME/])
;===============================================================================================================
;===============================================================================================================
;===============================================================================================================
;==========================================================================================
;==========================================================================================
qnet = new((/dimsizes(SH(:,0,0)),5/),float,SH@_FillValue)
qnet(:,4) = wgt_areaave_Wrap(oaflux(time |:,{lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0) 
qnet(:,3) = wgt_areaave_Wrap(SW(time |:,{lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0) 
qnet(:,2) = wgt_areaave_Wrap(LW(time |:, {lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0)
qnet(:,1) = wgt_areaave_Wrap(LH(time |:,{lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0) 
qnet(:,0) = wgt_areaave_Wrap(SH(time |:,{lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0)
;
qnet!0 = "time"
qnet!1 = "lev"
qnet&time = SH&time
qnet&lev = fspan(0,4,5)
;print(qnet&time)
print(qnet&lev)
printVarSummary(qnet)
;=========================================================================================
;=========================================================================================
;
qnet_s     = onedtond(ndtooned(qnet),(/nyears,nmos,5/))
qnet_s!0 = "time"
qnet_s!1 = "month"
qnet_s!2 = "lev"
printVarSummary(qnet_s)
delete(qnet)


;
;printMinMax(qnet_s,True)
;qnet_s = qnet_s/conform(qnet_s,days,1)
;printMinMax(qnet_s,True)
;print(xxxxxxxxxxxxx)
;
qnet_s = dtrend_n(qnet_s,False,0)
;
std_qnet_s = dim_stddev_n_Wrap(qnet_s,0)
;print(std_qnet_s(month|:,lev|1))
printVarSummary(std_qnet_s)
;print(xxxxxxxxxxxxx)
;
;qnet_s = dim_standardize_n(qnet_s,0,0)      ; normalize 
do ii=0,4,1
print(std_qnet_s(:,ii))
end do
delete(ii)
;print(xxxxxxxxxxxxxxxxxx)
;=========================================================================================
;=========================================================================================

;=========================================================================================
;==============================================================
;=========================================================================================
;=========================================================================================
print("==============================================================")
print("==============================================================")
;===============================================================================================================
;======================================
  wks = gsn_open_wks("eps","FIG_S7")         ; send graphics to PNG file
  plot = new(4,graphic)                ; create graphic array ; only needed if paneling
  
;===============================================================

   r2               = True                        ; plot mods desired
   r2@gsnDraw       = False                       ; do not draw
   r2@gsnFrame      = False                       ; do not advance frame
   r2@tiMainString  = ""                       ; title
   r2@tiXAxisFontHeightF = 0.015
   r2@tiYAxisFontHeightF = 0.015
   r2@tiYAxisString = "s.d. [ITCZ (~S~o~N~N)]"                  ; yaxis
   r2@tiXAxisString = "Calendar month"                  ; yaxis
  ; r2@tiYAxisString = "SST (~S~o~N~C)"                  ; yaxis
   ;r2@tiYAxisAngleF  =  270.
   r2@tmXBLabelFontHeightF = 0.017
   r2@tmYLLabelFontHeightF = 0.017
   r2@trYLog              = False                ; log scaling
   r2@trXLog              = False               ; log scaling
   ;r2@trXMinF             = 1940.                 ; manually set lower limit
   ;r2@trXMaxF             = 2020.                ;   "          upper
  r2@trYMinF             =   0.                 ; manually set lower limit
  r2@trYMaxF             =    2.4          

   r2@trXMinF             =  0.                 ; manually set lower limit
   r2@trXMaxF             =  11.                 ;   "          upper
  ; r2@gsnYRefLineThicknesses = (/1/)
  ; r2@gsnYRefLineDashPatterns = (/0./)
   ;r2@xyLineColors        = (/"springgreen4","springgreen","transparent"/)
   r2@xyLineColors        = xy_colors_thermocline
   r2@tmXTBorderOn            =  True
   r2@tmXTOn                  =  True
      
   r2@tmYLFormat = "0@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
  ; r2@gsnXRefLine           = (/5./) 
  ; r2@gsnYRefLineDashPattern = 16

   ;r2@tmYLMode        = "Explicit"	
  ; r2@tmYLValues      = (/"-0.5","0.0","0.5","1.0"/)
   ;r2@tmYLLabels      =  r2@tmYLValues
   ;r2@tmYLMinorValues = fspan(-1.0,1.0,11)
   ;print(r2@tmYLMinorValues)
   
   r2@vpWidthF                = 0.60
   r2@vpHeightF               = 0.35
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
   r2@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
   
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
   r2@pmLegendParallelPosF   = 0.90                 ; move units right
   r2@pmLegendOrthogonalPosF = -0.45                ; move units down
   r2@pmLegendWidthF         = 0.12                ; Change width and
   r2@pmLegendHeightF        = 0.10                ; height of legend.
   r2@lgLabelFontHeightF     = 0.010                 ; change font height
   r2@lgTitleOn              = False               ; turn on legend title
   r2@lgTitleString          = "SST indices"           ; create legend title
   r2@lgTitleFontHeightF     = .003                ; font of legend title
   r2@lgPerimOn              = False               ; turn off box around
   r2@xyExplicitLegendLabels = (/"~F10~Q~B~~F21~SH~N","~F10~Q~B~~F21~LH~N","~F10~Q~B~~F21~LW~N","~F10~Q~B~~F21~SW~N","~F10~Q~B~~F21~net~N","SODA331"/)         ; explicit labels 
;==================================================================================
  

;-----------------------------------------------------------------------------
;=============================================================================

r2@trYMinF             =   0.                 ; manually set lower limit
  r2@trYMaxF             =   21.           
r2@tiYAxisString = "s. d. [~F10~Q~B~~F21~net~N~ (W m~S~-2~N~)]"                  ; yaxis
plot(1)  =  gsn_csm_xy(wks,xlabel,std_qnet_s(lev|:,month|:),r2)


 

 r2@xyLineColors        = xy_colors_dheating
     
r2@xyExplicitLegendLabels = (/"MERRA2","JRA55","ERAI","ERA5","CFSR","SODA331"/)         ; explicit labels

 r2@tiYAxisString = "s.d. [ITCZ (~S~o~N~N)]"        
r2@trYMinF             =   0.                 ; manually set lower limit
  r2@trYMaxF             =    2.4                    ; yaxis
 plot(2)  =  gsn_csm_xy(wks,xlabel,std_i_itcz_s(lev|:,month|:),r2)



r2@trYMinF             =   0.                 ; manually set lower limit
  r2@trYMaxF             =  15.           
r2@tiYAxisString = "s.d. [Zonal wind stress (~F10~e~S~~F21~-3~N~ N m~S~-2~N~)]"                  ; yaxis
 r2@tmYLFormat = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
plot(3)  =  gsn_csm_xy(wks,xlabel,std_wustr_s(lev|:,month|:),r2)

r2@tmYLFormat = "0@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
r2@tiYAxisString = "s.d. [Precipitation (mm day~S~-1~N~)]"                  ; yaxis
r2@xyLineColors        = (/"grey35","yellow3","red","springgreen4","springgreen","black"/)
r2@trYMinF             =   0.                 ; manually set lower limit
  r2@trYMaxF             =    2.4       
  r2@pmLegendDisplayMode    = "Never"            ; turn on legend   
 plot(0)  =  gsn_csm_xy(wks,xlabel,std_rain_s,r2)

;  plotA1  =  gsn_csm_x2y2(wks,xlabel,xlabel,r_tcline(lev|:,month|:),mean_i_itcz_s(lev|:,month|:),r0,r1)

 r2@pmLegendDisplayMode    = "Never"            ; turn on legend
 r2@trYMinF             =   -1.5                 ; manually set lower limit
 r2@trYMaxF             =    1.5           
 r2@tmYLFormat = "0@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
; plotA1  =  gsn_csm_xy(wks,xlabel,r_itcz(lev|:,month|:),r2)


  ;r0@trYMinF             =    -2.                 ; manually set lower limit
  ;r0@trYMaxF             =    10.           
  ;r0@tiYAxisString = "Wind stress feedback (e~S~-3~N~ N m~S~-2~N~)"                  ; yaxis
;
 ; r1@trYMinF             =   -.3                 ; manually set lower limit
  ;r1@trYMaxF             =    1.5           
  ;r1@xyLineThicknesses   = (/4.5,1.5,1.5,1.5,1.5,1.5/)       ; Define line thicknesses 
  ;r1@tiYAxisString = "Precipitation feedback (mm day~S~-1~N~)"                  ; yaxis
;  r0@pmLegendOrthogonalPosF = -0.3                ; move units down
 ; r0@pmLegendParallelPosF   = 0.12                 ; move units right
 ; r0@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
 ; r0@tiXAxisString = "Time lag (month)"                  ; yaxis
 ;plotA1  =  gsn_csm_x2y2(wks,xlabel,xlabel,r_tcline(lev|:,month|:),d_arain(month|:),r0,r1)
;r0@pmLegendDisplayMode    = "Never"            ; turn on legend
;------------------------------------------------------------------------
;----------------------------------------------------------------------------------------------------------------
 

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
  txres5@txFontHeightF = 0.020   ; Use same font height as left axis
  txres5@txFont                = "Helvetica-Bold"
 ; txid = gsn_create_text(wks, "SST anomaly (~S~o~N~C)", txres)
   txid5 = gsn_create_text(wks, "a", txres5);
; Move text string to center/right edge of plot.;
  amres5                  = True
  amres5@amParallelPosF   =  -0.64 ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.
  amres5@amOrthogonalPosF =  -0.52    ; This is the center of the plot.
  amres5@amJust           = "TopLeft"
  annoid5i = gsn_add_annotation(plot(0),txid5,amres5)  ; Attach string to plot
 txid5b = gsn_create_text(wks, "b", txres5);
 annoid5ib = gsn_add_annotation(plot(1),txid5b,amres5)  ; Attach string to plot

 ;txres5@txFontHeightF = 0.035   ; Use same font height as left axis
 ;amres5@amParallelPosF   =  -0.67 ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.

 txid5c = gsn_create_text(wks, "c", txres5);

 txid5d = gsn_create_text(wks, "d", txres5);
 txid5e = gsn_create_text(wks, "e", txres5);
 txid5f = gsn_create_text(wks, "f", txres5);

 
 annoid5ic = gsn_add_annotation(plot(2),txid5c,amres5)  ; Attach string to plot
 annoid5id = gsn_add_annotation(plot(3),txid5d,amres5)  ; Attach string to plot
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
 
; dum0 = gsn_add_polyline(wks,plot(2),lonEQ,latEQ,pres0)
 ;dum1 = gsn_add_polyline(wks,plot(3),lonEQ,latEQ,pres0)
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
  dumL2 = gsn_add_polyline(wks,plot(2), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
  dumL3 = gsn_add_polyline(wks,plot(3), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
;
  x0 = 0.
  x1 = 11.
  y0 = 0.
 ; dumL0x = gsn_add_polyline(wks,plot(0), (/x0,x1/),(/y0,y0/) ,pres1)
 ; dumL1x = gsn_add_polyline(wks,plot(1), (/x0,x1/),(/y0,y0/) ,pres1)
;===============================================================
;===============================================================
;===============================================================
;===============================================================
;=====================================================================
;===============================================================
   resP             			  = True                        ; panel mods desired
   resP@gsnDraw              		  = True       ; don't draw yet
   resP@gsnFrame         		  = False        ; don't advance frame yet
   resP@gsnPanelYWhiteSpacePercent       = 3.0    ; default is 1.0
   resP@gsnPanelXWhiteSpacePercent       = 3.0    ; default is 1.0
   resP@gsnPanelTop                      = 0.95
   resP@gsnPanelBottom                   = 0.20
   resP@gsnPanelLeft                     = 0.05
   resP@gsnPanelRight                    = 0.95
   resP@gsnPanelFigureStringsPerimOn 	  = False 
   resP@gsnPaperOrientation 		  = "portrait"
   resP@gsnPanelDebug      		  = True
   resP@amJust   			  = "BottomLeft"
   resP@gsnPanelFigureStringsPerimOn     = False
   resP@gsnPanelFigureStringsFontHeightF = 0.010
   resP@gsnPanelFigureStringsBackgroundFillColor = "transparent"
   ;resP@gsnPanelFigureStrings= (/"lag = -1","lag = +1","c","d","e","f","g","h","i","j","k","l"/) ; add strings to panel
   resP@gsnMaximize = True                        ; blow up plot
  gsn_panel(wks,plot(:),(/2,2/),resP)               ; create panel plots

  ;resP@gsnPanelCenter   = False                 ; =True .... must move labels
  ;resP@gsnPanelRowSpec  = True                  ; tell panel what order to plot
  ;gsn_panel(wks,plot,(/1,2/),resP)
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
  res_P@gsnPanelXWhiteSpacePercent 	= 1.0    ; default is 1.0
  res_P@gsnPanelTop                      = 0.96
  res_P@gsnPanelBottom                   = 0.65
  res_P@gsnPanelLeft                     = 0.08
  res_P@gsnPanelRight                    = 0.88
  res_P@gsnPaperOrientation = "portrait"   ; force portrait
 ; res_P@gsnPanelFigureStrings   = (/"(b)","~F0~Peak~F34~;~F0~12 yr~S~-1~N~"/)
  res_P@amJust   			  = "TopLeft"
  res_P@gsnPanelFigureStringsFontHeightF = 0.020
  res_P@gsnPanelFigureStringsPerimOn 	  = True 
  res_P@gsnPanelFigureStringsBackgroundFillColor = "grey"
  ;res_P@gsnPanelFigureStringsJust    =     "CenterRight"
  ;res_P@amJust               =   "CenterRight"
  res_P@gsnMaximize          = True            ; large format
 ; gsn_panel(wks,plot(0:1),(/1,2/),res_P)     ; now draw as one plot
;======================================================
;======================================================
;===============================================================
;===============================================================
;***********************************************
;==============================================================
frame(wks)
end
