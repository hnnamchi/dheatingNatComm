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



xy_colors_dheating        = (/"red","slateblue","deepskyblue","black","darksalmon","black"/)
;xy_colors_thermocline        = (/"blue","yellow2","red","slateblue","deepskyblue","black"/)
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

print(cfsr_s(:,:,0))
;
delete(cfsr)
;---------------------------
;
temp = cfsr_s(:,:,0)
printVarSummary(temp)
temp = dtrend_n(temp,False,0)

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


;=========================================================================================
;=========================================================================================
;
;
d_cfsr = escorc(cfsr_s(month|:,time|:),temp(month|:,time|:))
d_cfsr!0 = "month"
print(d_cfsr)
delete(temp)
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
d_sst = escorc(temp(month|:,time|:),a_sst_s(month|:,time|:))
d_sst!0 = "month"
print(d_sst)
delete(temp)
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
d_aviso = escorc(th_s(month|:,time|:),temp(month|:,time|:))
d_aviso!0 = "month"
print(d_aviso)
;==============================================================
;==============================================================


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

;==============================================================
; compute climatology and Anomalies
;==============================================================
  oisstClm   = clmMonTLL(oisst)                             ; (12,lat,lon)
  printVarSummary(oisstClm)

  oisst  = calcMonAnomTLL(oisst,oisstClm)                ; (time, lat,lon)
  printVarSummary(oisst)         
  printMinMax(oisst,True)
;=================================================================

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
aOISST_s = wgt_areaave_Wrap(oisst_s(time|:,month|:, {lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0)	
printVarSummary(aOISST_s)

;print(xxxxxxxxxxxxxxxxxxxxxxxxxx)
;==========================================================================================
;==========================================================================================
;===============================================================================
;===============================================================================
;===============================================================================
;===============================================================================
;=========================================================================================
;===============================================================================
;===============================================================================
;yrStrt = 1982
;yrLast = 2015
;


f0      = addfile("atl3_d23_soda342_oras5_oras4_godas_gecco2_en4_1982_2015_basin.nc", "r")
TIME   = f0->time    
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
d23    = f0->d23(iYYYY,:)
printVarSummary(d23) 
printMinMax(d23,True)

  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set the time dimension


;=========================================================================================
;=========================================================================================
;
d23_s     = onedtond(ndtooned(d23),(/nyears,nmos,6/))
d23_s!0 = "time"
d23_s!1 = "month"
d23_s!2 = "lev"
printVarSummary(d23_s)
delete(d23)
;
d23_s = dtrend_n(d23_s,False,0)
d23_s = dim_standardize_n(d23_s,0,0)  
;print(xxxxxxxxxxxxxx)
;
;=========================================================================================
;=========================================================================================
;
printVarSummary(aOISST_s)
printVarSummary(d23_s)
;
tConform = conform_dims(dimsizes(d23_s),aOISST_s,(/0,1/)) 
printVarSummary(tConform)
delete(aOISST_s)
copy_VarMeta(d23_s,tConform) 
aOISST_s = tConform
delete(tConform)
printVarSummary(aOISST_s)
;


r_d23 = escorc(aOISST_s(month|:,lev|:,time|:),d23_s(month|:,lev|:,time|:))
r_d23!0 = "month"
r_d23!1 = "lev"
;print(r_d23)
printVarSummary(r_d23)
;
r_d23 = (r_d23^2)*100.
do ii=0,4,1
print(r_d23(:,ii))
end do
delete(ii)
;==========================================================================================
;==========================================================================================
;===============================================================================
;===============================================================================
;yrStrt = 1982
;yrLast = 2015
;


f0      = addfile("atl3_d20_soda342_oras5_oras4_godas_gecco2_en4_1982_2015_basin.nc", "r")
TIME   = f0->time    
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
d20    = f0->d20(iYYYY,:)
printVarSummary(d20) 
printMinMax(d20,True)

  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set the time dimension


;=========================================================================================
;=========================================================================================
;
d20_s     = onedtond(ndtooned(d20),(/nyears,nmos,6/))
d20_s!0 = "time"
d20_s!1 = "month"
d20_s!2 = "lev"
printVarSummary(d20_s)
delete(d20)
;
d20_s = dtrend_n(d20_s,False,0)
d20_s = dim_standardize_n(d20_s,0,0)  
;print(xxxxxxxxxxxxxx)
;
;=========================================================================================
;=========================================================================================
;

r_d20 = escorc(aOISST_s(month|:,lev|:,time|:),d20_s(month|:,lev|:,time|:))
r_d20!0 = "month"
r_d20!1 = "lev"
;print(r_d20)
printVarSummary(r_d20)
;
r_d20 = (r_d20^2)*100.
do ii=0,4,1
print(r_d20(:,ii))
end do
delete(ii)
;==========================================================================================
;==========================================================================================





;f0      = addfile("atl3_ssh_soda342_oras5_oras4_godas_gecco2_1982_2015_basin.nc", "r")
;
f0      = addfile("/home/hyacinth/work/geomar/sstP/R0/a_ssh_soda342_oras5_oras4_godas_gecco2_1982_2015_basin.nc", "r")
;
TIME   = f0->time    
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
ssh    = f0->ssh(iYYYY,:)
printVarSummary(ssh) 
printMinMax(ssh,True)
;print(xxxxxxxxxxxx)
 
  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set the time dimension

;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
 ssh_d = new(dimsizes(ssh),typeof(ssh),ssh@_FillValue)
 x   = ssh&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =ssh(tt::nmos,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    ssh_d(tt::nmos,:) = work
 end do

 printVarSummary(ssh_d)
 ssh = ssh_d

;===
delete([/tt,x,x1,ssh_d,work/])
;=========================================================================================
;=========================================================================================
;
ssh_s     = onedtond(ndtooned(ssh),(/nyears,nmos,5/))
ssh_s!0 = "time"
ssh_s!1 = "month"
ssh_s!2 = "lev"
printVarSummary(ssh_s)
delete(ssh)
;
ssh_s = dtrend_n(ssh_s,False,0)
;
str_ssh_s = dim_stddev_n_Wrap(ssh_s,0)
print(str_ssh_s)
;
x = new((/nmos,6/),float,str_ssh_s@_FillValue)
x(:,0:4) = str_ssh_s	
x(:,5) = std_aviso	;dynamic topography
printVarSummary(x)
print(x(:,5))
delete(str_ssh_s)
str_ssh_s = x
delete(x)
;
str_ssh_s!0 = "month"
str_ssh_s!1 = "lev"
;=========================================================================================
;=========================================================================================
;
;ssh_s = dim_standardize_n(ssh_s,0,0)      ; normalize 
;=========================================================================================
;=========================================================================================

;
;
ssh_s = dim_standardize_n(ssh_s,0,0)      ; normalize
;
;r_ssh = regCoef_n(aOISST,ssh_s,0,0)
r_ssh = escorc(aOISST_s(month|:,lev|0:4,time|:),ssh_s(month|:,lev|:,time|:))
r_ssh!0 = "month"
r_ssh!1 = "lev"
printVarSummary(r_ssh)
print(r_ssh(:,4))
;=========================================================================================
;=========================================================================================
;=========================================================================================
;=========================================================================================
x = new((/nmos,6/),float,r_ssh@_FillValue)
x(:,0:4) = r_ssh	
x(:,5) = d_aviso	;dynamic topography
printVarSummary(x)
print(x(:,5))
delete(r_ssh)
r_ssh = x
delete(x)
;
;
;
r_ssh = (r_ssh^2)*100.

do ii=0,5,1
print(r_ssh(:,ii))
end do
delete(ii)

;=========================================================================================
;=========================================================================================
;============================================================================================
;=========================================================================================
;=========================================================================================
;
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
r_diab = escorc(aOISST_s(month|:,lev|0:3,time|:),diab_s(month|:,lev|:,time|:))
;r_diab = regCoef_n(aOISST,diab_s,0,0)
r_diab!0 = "month"
r_diab!1 = "lev"
printVarSummary(r_diab)
;=========================================================================================
;=========================================================================================
x = new((/nmos,5/),float,r_ssh@_FillValue)
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
;
r_diab = (r_diab^2)*100.
do ii=0,4,1
print(r_diab(:,ii))
end do
delete(ii)



print("==============================================================")
print("==============================================================")
;===============================================================================================================
;======================================
  wks = gsn_open_wks("eps","FIG_2")         ; send graphics to PNG file
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
  r2@trYMaxF             =   80.           

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

   r2@tmYLMode        = "Explicit"	
   r2@tmYLValues      = (/"-0.5","0.0","0.5","1.0"/)
   r2@tmYLLabels      =  r2@tmYLValues
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
 
   r2@xyLineThicknesses   = (/1.5,1.5,1.5,1.5,1.5,1.5/)       ; Define line thicknesses 
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
   r2@pmLegendParallelPosF   = 0.88                 ; move units right
   r2@pmLegendOrthogonalPosF = -0.50                ; move units down
   r2@pmLegendWidthF         = 0.10                ; Change width and
   r2@pmLegendHeightF        = 0.12                ; height of legend.
   r2@lgLabelFontHeightF     = 0.010                 ; change font height
   r2@lgTitleOn              = False               ; turn on legend title
   r2@lgTitleString          = "SST indices"           ; create legend title
   r2@lgTitleFontHeightF     = .003                ; font of legend title
   r2@lgPerimOn              = False               ; turn off box around
   
   r2@xyExplicitLegendLabels = (/"MERRA2","JRA55","ERAI","ERA5","CFSR","SODA3"/)         ; explicit labels

;==================================================================================
 
;========================================================================
;========================================================================
;========================================================================
;----------------------------------------------------------------------------------------------------------------
   delete(r2@tmYLValues)
   delete(r2@tmYLLabels)     
   ;r2@tmYLValues      = (/"2","4","6","8","10"/)
   ;r2@tmYLLabels      =  r2@tmYLValues   
   r2@tiYAxisString = "Explained variance, ~F33~D~F10~Q  (%)"
   r2@pmLegendParallelPosF   = 0.12                 ; move units right
   ;r2@tiYAxisString = "s.d. (m)"
   ;r2@tiMainString  = "Thermocline feedback"                       ; title
   plot(0)  =  gsn_csm_xy(wks,xlabel,r_diab(lev|:,month|:),r2)
   r2@xyLineColors  =  xy_colors_thermocline
   r2@tiYAxisString = "Explained variance, SSH (%)"
   r2@xyExplicitLegendLabels = (/"SODA3","ORAS5","ORAS4","GODAS","GECCO2","AVISO"/)         ; explicit labels 
   
   plot(1)  =  gsn_csm_xy(wks,xlabel,r_ssh(lev|:,month|:),r2)
   r2@tiYAxisString = "Explained variance, ~F10~Z~B~~F21~20~N~ (%)"
   r2@xyExplicitLegendLabels = (/"SODA3","ORAS5","ORAS4","GODAS","GECCO2","EN4"/)         ; explicit labels 
   plot(2)  =  gsn_csm_xy(wks,xlabel,r_d20(lev|:,month|:),r2)
   r2@tiYAxisString = "Explained variance, ~F10~Z~B~~F21~23~N~ (%)"
   plot(3)  =  gsn_csm_xy(wks,xlabel,r_d23(lev|:,month|:),r2)
;=====================================================================
;=====================================================================
;-----------------------------------------------------------------------------
;=====================================================================
;=====================================================================
nplots = 4
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



DATA_NAMES = (/"~F33~D~F21~Q explained variance" ,"SSH explained variance","Z~B~20~N~ explained variance","Z~B~23~N~ explained variance","SODA3","Residual term"/)
txres@txFontHeightF = 0.017   ; Use same font height as left axis
txres@txFont                = "Helvetica"
 amres@amParallelPosF   =   -0.60 ; 0.5 is the right edge of the plot, so
 amres@amOrthogonalPosF =  -0.50    ; This is the center of the plot.
create_tt2 = gsn_create_text(wks,DATA_NAMES,txres);
tt_p1 = new(9,graphic) 
;---Attach the letters
  ;   do i = 0,nplots-1,1
     ; tt_p1(i) = gsn_add_annotation(plot(i),create_tt2(i),amres)
      
 ;    end do 
 ;delete(i)
;======================================================================== 
;========================================================================
  x_0 = 5.
  x_1 = 5.
  y_0 = -20. 
  y_1 = 100.  

;---Resources for the polyline
  pres1                  = True           ; polyline mods desired
  pres1@gsLineThicknessF = 1.0            ; line thickness   
  pres1@gsLineColor      = "black"         ; color of lines
  pres1@gsLineDashPattern = 2
  pres1@tfPolyDrawOrder   = "postDraw"
  pres1@gsLineColor      = "black"         ; color of lines
;---------- 
 
  dum0 = new(9,graphic)  
  dum1 = new(9,graphic)    
  
;---Attach the VERTICAL polyline
     do i = 0,nplots-1,1
      dum0(i) = gsn_add_polyline(wks,plot(i),(/x_0,x_1/),(/y_0,y_1/),pres1)
     end do 
	delete(i)
;
  x0 = 0.
  x1 = 11.
  y0 = 0.
  ;---Attach the ZERO polyline
     do i = 0,nplots-1,1
      dum1(i) = gsn_add_polyline(wks,plot(i), (/x0,x1/),(/y0,y0/) ,pres1)
     end do 
	delete(i)
;
;===============================================================
;===============================================================
;=====================================================================
;===============================================================
   resP             			  = True                        ; panel mods desired
   resP@gsnDraw              		  = True       ; don't draw yet
   resP@gsnFrame         		  = False        ; don't advance frame yet
   resP@gsnPanelYWhiteSpacePercent       = 0.0    ; default is 1.0
   resP@gsnPanelXWhiteSpacePercent       = 0.0    ; default is 1.0
   resP@gsnPanelTop                      = 0.99
   resP@gsnPanelBottom                   = 0.20
   resP@gsnPanelLeft                     = 0.01
   resP@gsnPanelRight                    = 0.99
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
  ; resP@gsnPanelFigureStrings= (/"lag = 0","lag = -1","c","d","e","f","g","h","i","j","k","l"/) ; add strings to panel
   resP@gsnMaximize = True                        ; blow up plot
  gsn_panel(wks,plot(:),(/2,2/),resP)               ; create panel plots
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
  ;gsn_panel(wks,plot(0:3),(/2,2/),res_P)     ; now draw as one plot
;======================================================
;======================================================
;===============================================================
;===============================================================
;***********************************************
;==============================================================
frame(wks)
end
