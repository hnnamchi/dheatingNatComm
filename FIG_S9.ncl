;*************************************************
;*************************************************
begin
;==================================================


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

xy_colors_thermocline        = (/"red","yellow2","deepskyblue","blue","deeppink4","black"/)
xy_colors_dheating        = (/"red","blue","deepskyblue","black","deeppink4","black"/)
;xy_colors_dheating        = (/"red","greenyellow","deepskyblue","black","darksalmon","black"/)
;xy_colors_thermocline        = (/"blue","yellow3","red","springgreen4","springgreen","black"/)


yrStrt = 1982
yrLast = 2015

yyyymm = yyyymm_time(yrStrt,yrLast,"integer")
;print(yyyymm) 

ndata = 6

;=============================================================================================================
;===============================================================================================================
;===============================================================================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
f0      = addfile("/home/hyacinth/work/geomar/sstP/R0/atl3_soda342_oras5_oras4_godas_gecco2_mld_1984_2009.nc", "r")
  TIME   = f0->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt.and. YYYY.le.yrLast)
  mld    = f0->mld(iYYYY,0:ndata-1)
 
  printVarSummary(mld) 
  printMinMax(mld,True)

    
  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set the time dimension


;=============================================================================================================
;=============================================================================================================
;===============================================================================================================
;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
 mld_d = new(dimsizes(mld),typeof(mld),mld@_FillValue)
 x   = mld&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =mld(tt::nmos,:)
    work = dtrend_msg_n(x1,work,False,False,0)
    mld_d(tt::nmos,:) = work
 end do

 printVarSummary(mld_d)
 mld = mld_d

;===
delete([/tt,x,x1,mld_d,work/])
 printMinMax(mld,True)
;print(xxxxxxx)
;=========================================================================================
;=========================================================================================
;=========================================================================================
;=========================================================================================
;
nyears = dimsizes(mld&time)/nmos
print(nyears)

mld_s     = onedtond(ndtooned(mld),(/nyears,nmos,ndata/))
mld_s!0 = "time"
mld_s!1 = "month"
mld_s!2 = "lev"
printVarSummary(mld_s)
delete(mld)


;
;mld_s = dtrend_n(mld_s,False,0)
;
;mld_s = dim_standardize_n(mld_s,0,0)      ; normalize 
;
std_mld_s = dim_stddev_n_Wrap(mld_s,0)
printVarSummary(std_mld_s)
;
;
mean_mld_s = dim_avg_n_Wrap(mld_s,0)
printVarSummary(mean_mld_s)
;
print("========================================")
do ii=0,ndata-1,1
print(std_mld_s(:,ii))
end do
delete(ii)
print("========================================")


;===============================================================================

print("==============================================================")
print("==============================================================")
;===============================================================================================================
;======================================
  wks = gsn_open_wks("eps","FIG_S9")         ; send graphics to PNG file
  plot = new(2,graphic)                ; create graphic array ; only needed if paneling
 
;-----------------------------------------------------------------------------
;=============================================================================
;===============================================================
;=============================================================== 
;***********************************************
   r2               = True                        ; plot mods desired
   r2@gsnDraw       = False                       ; do not draw
   r2@gsnFrame      = False                       ; do not advance frame
   r2@tiMainString  = ""                       ; title
   r2@tiXAxisFontHeightF = 0.013
   r2@tiYAxisFontHeightF = 0.013
   r2@tiYAxisString = "Atl3 sea surface temperature (K)"                  ; yaxis
   ;r2@tiYAxisString = "Mean [h (m)]" 
   r2@tiYAxisString = "Mean [ocean mixed lay depth (m)]" 
   r2@tiXAxisString = "Calendar month"                  ; yaxis
   ;r2@tiYAxisString = "Q~B~net~N~ terms (W m~S~-2~N~)"                  ; yaxis
  ; r2@tiYAxisString = "SST (~S~o~N~C)"                  ; yaxis
   ;r2@tiYAxisAngleF  =  270.
   r2@tmXBLabelFontHeightF = 0.013
   r2@tmYLLabelFontHeightF = 0.013
   r2@trYLog              = False                ; log scaling
   r2@trXLog              = False               ; log scaling
   r2@trYReverse           = True               ; log scaling
   ;r2@trXMinF             = 1940.                 ; manually set lower limit
   ;r2@trXMaxF             = 2020.                ;   "          upper
  r2@trYMinF             =   0.                 ; manually set lower limit
  r2@trYMaxF             =   0.8           

   r2@trXMinF             =  0.                 ; manually set lower limit
   r2@trXMaxF             =  11.                 ;   "          upper
  ; r2@gsnYRefLineThicknesses = (/1/)
  ; r2@gsnYRefLineDashPatterns = (/0./)
   r2@xyLineColors        = xy_colors_thermocline
   

   r2@tmXTBorderOn            =  True
   r2@tmXTOn                  =  True
      
   r2@tmYLFormat = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
  ; r2@gsnXRefLine           = (/5./) 
  ; r2@gsnYRefLineDashPattern = 16

   ;r2@tmYLMode        = "Explicit"	
   ;r2@tmYLValues      = (/"-0.5","0.0","0.5","1.0"/)
   ;r2@tmYLLabels      =  r2@tmYLValues
   ;r2@tmYLMinorValues = fspan(-1.0,1.0,11)
   ;print(r2@tmYLMinorValues)
   
   r2@vpWidthF                = 0.50
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
   ;r2@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
   r2@tmXBLabels = (/"Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
   
   ;r2@tmXBMinorOn          = False                ; Turn off x-axis (bottom) minor ticks
   ;r2@tmXTMinorOn          = False                ; Turn off x-axis (bottom)minor ticks
   ;r2@tmXTOn               = False                ; Turn off x-axis (top) major ticks
   ;r2@tmYRMinorOn          = False                ; Turn off y-axis (right) minor ticks
   ;r2@tmYLMinorOn          = False                ; Turn off y-axis (left) minor ticks
   ;r2@tmYROn               = False                ; Turn off y-axis (right) major ticks

  ;===legends
   ;===legends
   r2@pmLegendDisplayMode    = "Never"            ; turn on legend
   r2@pmLegendSide           = "Bottom"               ; Change location of 
   r2@pmLegendParallelPosF   = 0.15                 ; move units right
   r2@pmLegendOrthogonalPosF = -0.55                ; move units down
   r2@pmLegendWidthF         = 0.12                ; Change width and
   r2@pmLegendHeightF        = 0.10                ; height of legend.
   r2@lgLabelFontHeightF     = 0.010                 ; change font height
   r2@lgTitleOn              = False               ; turn on legend title
   r2@lgTitleString          = "SST indices"           ; create legend title
   r2@lgTitleFontHeightF     = .003                ; font of legend title
   r2@lgPerimOn              = False               ; turn off box around
   r2@xyExplicitLegendLabels = (/"SODA3","ORAS5","ORAS4","GODAS","GECCO2","EN4"/)         ; explicit labels 


;
;r2@tiXAxisString = "Calendar month"                  ; yaxis
r2@tiXAxisString = "Calendar month"                  ; yaxis


;----------------------------------------------------------------------------------------------------------------

   



  r2@trYMinF             =    0.                 ; manually set lower limit
  r2@trYMaxF             =   50.           

  r2@pmLegendDisplayMode    = "Always"            ; turn on legend
  ;r2@xyDashPatterns      = (/1.,0.,0.,0.,0.,0./)           ; Dash patterns 
  ;r2@tiYAxisString = "Q~B~net~N~ feedback (W m~S~-2~N~)"                  ; yaxis
  ;r2@tiYAxisString = "~F18~s~F21~T/~F18~s~F21~t~N~ terms (e~S~-6~N~ K s~S~-1~N~)"                  ; yaxis
 
  plot(0)  =  gsn_csm_xy(wks,xlabel,mean_mld_s(lev|:,month|:),r2)

  
 
  r2@trYMaxF             =   8.        
  r2@tiYAxisString = "s.d. [ocean mixed lay depth (m)]"    
  ;r2@pmLegendDisplayMode    = "Never"            ; turn on legend
  plot(1)  =  gsn_csm_xy(wks,xlabel,std_mld_s(lev|:,month|:),r2)
 
 

;=====================================================================
;=====================================================================
nplots = 2
;-----------------------------------------------------------------------------
  txres               = True
  txres@txAngleF      = 0.          ; Rotate string clockwise 
  txres@txFontHeightF = 0.023   ; Use same font height as left axis
  txres@txFont                = "Helvetica-Bold"
 ; Move text string to center/right edge of plot.;
  amres                  = True
  amres@amParallelPosF   =  -0.62 ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.
  amres@amOrthogonalPosF =  -0.52    ; This is the center of the plot.
  amres@amJust           = "TopLeft"
  txres@txFontHeightF = 0.018   ; Use same font height as left axis


 tt = (/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o"/)
create_tt = gsn_create_text(wks,tt,txres);
tt_p0 = new(9,graphic)                ; create graphic array ; only needed if paneling
;---Attach the letters
     do i = 0,nplots-1,1
      tt_p0(i) = gsn_add_annotation(plot(i),create_tt(i),amres)
     end do 
 delete(i)


DATA_NAMES = (/"Mean","Standard deviations","ERA5","ERAI","JRA55","MERRA2"/)
txres@txFontHeightF = 0.017   ; Use same font height as left axis
txres@txFont                = "Helvetica"
 amres@amParallelPosF   =  -0.47 ; 0.5 is the right edge of the plot, so
 amres@amOrthogonalPosF =  -0.40    ; This is the center of the plot.
create_tt2 = gsn_create_text(wks,DATA_NAMES,txres);
tt_p1 = new(9,graphic) 
;---Attach the letters
;     do i = 0,nplots-1,1
;      tt_p1(i) = gsn_add_annotation(plot(i),create_tt2(i),amres)
;     end do 
; delete(i)
;======================================================================== 
;========================================================================
  x_0 = 5.
  x_1 = 5.
  y_0 = -20. 
  y_1 = 50.  

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
;***********************************************
; panel plot only resources
  res_P                     		 = True            ; modify the panel plot
  res_P@gsnDraw            		 = True       ; don't draw yet
  res_P@gsnFrame            		 = False        ; don't advance frame yet
  res_P@gsnPanelYWhiteSpacePercent 	= 0.0    ; default is 1.0
  res_P@gsnPanelXWhiteSpacePercent 	= 2.0    ; default is 1.0
  ;res_P@gsnPanelTop                      = 0.96
  res_P@gsnPanelBottom                   = 0.10
  ;res_P@gsnPanelLeft                     = 0.02
  ;res_P@gsnPanelRight                    = 0.90
  ; res_P@gsnPanelFigureStrings   = (/"(b)","~F0~Peak~F34~;~F0~12 yr~S~-1~N~"/)
  res_P@amJust   			  = "TopLeft"
  res_P@gsnPanelFigureStringsFontHeightF = 0.020
  res_P@gsnPanelFigureStringsPerimOn 	  = True 
  res_P@gsnPanelFigureStringsBackgroundFillColor = "grey"
  res_P@gsnPaperOrientation 		  = "portrait"
  ;res_P@gsnPanelFigureStringsJust    =     "CenterRight"
  ;res_P@amJust               =   "CenterRight"
  res_P@gsnMaximize          = True            ; large format
  gsn_panel(wks,plot(0:),(/2,1/),res_P)     ; now draw as one plot
;======================================================
;======================================================
;===============================================================
;===============================================================
;***********************************************
;==============================================================
frame(wks)
end
