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


xy_colors_dheating        = (/"red","slateblue","deepskyblue","black","darksalmon","black"/)
;xy_colors_dheating        = (/"red","greenyellow","deepskyblue","black","darksalmon","black"/)
xy_colors_thermocline        = (/"blue","yellow3","red","springgreen4","springgreen","black"/)
;=========================================================================================
;=========================================================================================
;=========================================================================================
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
 ;  oisst = oisst - runave_n_Wrap (oisst,nave,opt,dim)   

  ; oisst = lonFlip(oisst)

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
temp = wgt_areaave_Wrap(oisst_s(time|:,month|:, {lon|210:270}, {lat | -5:5}),1.0, 1.0, 0)	;June 
printVarSummary(temp)
std_temp = dim_stddev_n_Wrap(temp,0)
print(std_temp)
temp = dim_standardize_n(temp,0,0)      ; normalize 

;print(xxxxxxxxxxxxxxxxxxxxxxxxxx)
;==========================================================================================
;==========================================================================================
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
;
;
xy = ustr(time|6:401,lev|:,lat|:,lon|:)
printVarSummary(xy)
;
delete(ustr)
ustr = xy
delete(xy)
;
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
tcline_s = wgt_areaave_Wrap(ustr_s(time |:,month|:,lev|:,{lon|160:210}, {lat | -5:5}),1.0, 1.0, 0)	
tcline_s  = tcline_s*1000.
printVarSummary(tcline_s)
;
;=========================================================================================
;=============================================================================================================
sst = new((/nyears+1,nmos,5/),float,temp@_FillValue)
sst(:,:,0) = temp(time |:,month|:)
sst(:,:,1) = temp(time |:,month|:)
sst(:,:,2) = temp(time |:,month|:)
sst(:,:,3) = temp(time |:,month|:)
sst(:,:,4) = temp(time |:,month|:)
;sst(:,:,5) = temp
sst!2 = "lev"
printVarSummary(sst)
printVarSummary(temp)
atemp = temp(:,11)	;December
atemp!0 = "time"
print(atemp)
delete(temp)
;print(xxxxxxxxxxx)
;=============================================================================================================
;=========================================================================================
;r_tcline = regCoef(sst(lev|:,month|:,time|:),tcline_s(lev|0:,month|:,time|:))
r_tcline = regCoef(atemp(time|0:32),tcline_s(lev|0:,month|:,time|:))
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
f      = addfile ("/home/hyacinth/work/datasets/rain/GPCP_precip.mon.mean.nc", "r")
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

rain    = short2flt(f->precip(iYYYY,:,:))
;rain = lonFlip(rain)
printVarSummary(rain)       ; variable overview

printMinMax(rain,True)
;print(xxxxxxxx)

rain!1  =  "lat"
rain!2  =  "lon"
delete([/YYYY,iYYYY,f,TIME/])
;==============================================================
;===========================
arain = wgt_areaave_Wrap(rain(time|:,{lon|160:270}, {lat | -5:5}),1.0, 1.0, 0)
printVarSummary(arain)
;
xy = arain(time|6:401)
printVarSummary(xy)
;
delete(arain)
arain = xy
delete(xy)
;
nyears = dimsizes(arain&time)/nmos
print(nyears)
;print(xxxxxxxxxxxx)
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
d_arain = regCoef(atemp(time|0:32),arain_s(month|:,time|0:32))
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
;=========================================================================================
;=========================================================================================
print("==============================================================")
print("==============================================================")
;===============================================================================================================
;======================================
  wks = gsn_open_wks("eps","FIG_4")         ; send graphics to PNG file
  plot = new(1,graphic)                ; create graphic array ; only needed if paneling
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
   r0@vpHeightF               = 0.50
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
  r0@pmLegendParallelPosF   = 0.12                 ; move units right
  r0@pmLegendOrthogonalPosF = -0.25                ; move units down
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
   r1@tiYAxisAngleF  =  270. 
   r1@xyLineColors        =  (/"deeppink4","deeppink4","deeppink4","deeppink4","deeppink4","deeppink4"/)
   r1@xyLineThicknesses   = (/1.5,1.5,1.5,1.5,1.5,1.5/)       ; Define line thicknesses 
   r1@xyDashPatterns      = (/1.,0.,0.,0.,0.,0./)           ; Dash patterns 
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
   r1@tmXTOn                  =  False
   r1@tmXTMode = "Explicit"
   r1@tmXTValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
   ;r1@tmXTLabels = (/"1", "2","3","4","5","6","7","8","10","10","11","12"/)
   
  r1@tmYRLabelFontColor = "deeppink4"   ; and right axes are same.
  r1@tiYAxisFontColor  = "deeppink4"
  r1@tmXTLabelsOn        = False    ; have tick mark labels
  r1@tmXTOn              = True      ; have tick marks
  r1@tmYRLabelFontColor = "deeppink4"   ; and right axes are same.
  r1@tiYAxisFontColor  = "deeppink4"
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

;-----------------------------------------------------------------------------
;=============================================================================

  r0@xyLineColors        =  (/"black","black","black","deeppink4","deeppink4","deeppink4"/)

 

 r0@xyLineColors        =	 (/"deeppink4","deeppink4","deeppink4","deeppink4","deeppink4","deeppink4"/)


 r0@xyDashPatterns      = (/0.,0.,0.,0.,0.,0./)           ; Dash patterns 
 r0@xyExplicitLegendLabels = (/"MERRA2","JRA55","ERAI","ERA5","CFSR","SODA331"/)         ; explicit labels
 r0@xyLineThicknesses   = (/1.5,1.5,1.5,1.5,1.5,1.5/)       ; Define line thicknesses 
 r0@xyLineColors        =  xy_colors_dheating


  r0@trYMinF             =      0.                 ; manually set lower limit
  r0@trYMaxF             =      30.           
  ;r0@tiYAxisString = "~F33~t~B~c~N~~F21~ (e~S~-3~N~ N m~S~-2~N~)"                  ; yaxis
  r0@tiYAxisString = "Zonal wind stress (~F10~e~S~~F21~-3~N~ N m~S~-2~N~)"                  ; yaxis
;
  r1@trYMinF             =    0.0                 ; manually set lower limit
  r1@trYMaxF             =    3.0           
  r1@xyLineThicknesses   = (/1.5,1.5,1.5,1.5,1.5,1.5/)       ; Define line thicknesses 
  r1@tiYAxisString = "Precipitation (mm day~S~-1~N~)"                  ; yaxis
  ;r0@pmLegendParallelPosF   = 0.12                 ; move units right
  r0@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
  r0@tiXAxisString = "Time lag (month)"                  ; yaxis
  r0@tmYLFormat = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
  r1@tmYRFormat = "0@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
  plotA1  =  gsn_csm_x2y2(wks,xlabel,xlabel,r_tcline(lev|:,month|:),d_arain(month|:),r0,r1)
  r0@pmLegendDisplayMode    = "Never"            ; turn on legend
;------------------------------------------------------------------------
;----------------------------------------------------------------------------------------------------------------
   delete(r0@xyLineColors)
   r0@xyLineColors        = xy_colors_dheating
   r0@xyMarkLineMode	 =   "Markers"
   r0@xyMarkers     =  (/16.,16.,16.,16.,16.,16./)
   r0@xyMarkerColors = r0@xyLineColors
   r0@xyMarkerSizeF   = 0.010
 
  ; delete(r1@xyLineColors)
  
   r1@xyMarkLineMode	 =   "Markers"
   r1@xyMarkers     =  (/16.,16.,16.,16.,16.,16./)
   r1@xyMarkerColors = r1@xyLineColors
   r1@xyMarkerSizeF   = 0.010
  plotB1  =  gsn_csm_x2y2(wks,xlabel,xlabel,r_tclined(lev|:,month|:),d_araind(month|:),r0,r1)
;
  overlay(plotA1,plotB1)
   plot(0) = plotA1
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
;  annoid5i = gsn_add_annotation(plot(0),txid5,amres5)  ; Attach string to plot
 txid5b = gsn_create_text(wks, "b", txres5);
; annoid5ib = gsn_add_annotation(plot(1),txid5b,amres5)  ; Attach string to plot

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
  ;dumL1 = gsn_add_polyline(wks,plot(1), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
 ; dumL2 = gsn_add_polyline(wks,plot(2), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
  ;dumL3 = gsn_add_polyline(wks,plot(3), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
;
  x0 = 0.
  x1 = 11.
  y0 = 0.
  dumL0x = gsn_add_polyline(wks,plot(0), (/x0,x1/),(/y0,y0/) ,pres1)
  ;dumL1x = gsn_add_polyline(wks,plot(1), (/x0,x1/),(/y0,y0/) ,pres1)
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
  gsn_panel(wks,plot(0:0),(/1,1/),res_P)     ; now draw as one plot
;======================================================
;======================================================
;===============================================================
;===============================================================
;***********************************************
;==============================================================
frame(wks)
end
