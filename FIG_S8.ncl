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



xy_colors_dheating        = (/"red","greenyellow","deepskyblue","black","darksalmon","black"/)
;xy_colors_thermocline        = (/"blue","yellow3","red","springgreen4","springgreen","black"/)
xy_colors_thermocline        = (/"blue","yellow3","springgreen4","springgreen","black","red"/)

yrStrt = 1984
yrLast = 2009

yyyymm = yyyymm_time(yrStrt,yrLast,"integer")
;print(yyyymm) 

ndata = 6
;=========================================================================================
;=========================================================================================
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


  tmp_annual_mean = month_to_annual_weighted(yyyymm,oisst,1)   ; tmp_annual_mean
  printVarSummary(tmp_annual_mean)
  print(tmp_annual_mean&year)

;=============================================================================================================
;=============================================================================================================

nyears = dimsizes(oisst&time)/nmos
print(nyears)
oisst_s     = onedtond(ndtooned(oisst),(/nyears,nmos,dimsizes(oisst&lat),dimsizes(oisst&lon)/))
copy_VarCoords(oisst(0,:,:),oisst_s(0,0,:,:))                             ; copy coordinate arrays from slp to newslp
oisst_s!0 = "time"
oisst_s!1 = "month"
;
printVarSummary(oisst_s)
;
;
oisst_s = dtrend_n(oisst_s,False,0)
;
;oisst_s = dim_standardize_n(oisst_s,0,0)      ; normalize 
;
;==========================================================================================
aoisst_s = wgt_areaave_Wrap(oisst_s(time |:,month|:,{lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0) 
printVarSummary(aoisst_s)
;
;=================================================
all_years = (/1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,\
		2005,2006,2007,2008,2009/)
print(all_years)
;
;
event_years = (/1984,-999,-999,1987,1988,-999,-999,1991,-999,1993,-999,1995,1996,-999,1998,1999,-999,-999,-999,-999,-999,\
		-999,-999,-999,2008,-999/)
print(event_years)
event_years@missing_value = -999
event_years@_FillValue = event_years@missing_value


;event_years = (/1984, 1987, 1988, 1991, 1993, 1995, 1996, 1998, 1999, 2008/)
;print(event_years)


iz = ind(.not.ismissing(event_years))  ; iz will have indices of non-missing values
                               ; iz(0)=0, iz(1)=2, iz(2)=4
print(iz)

nCases = dim_num(.not.ismissing(event_years)) 
print(nCases)
;

X = all_years(iz)   
print(X)   

xx = aoisst_s(iz,:)
printVarSummary(xx)

aSST = dim_avg_n_Wrap(xx,0)
print(aSST)

;===========
aSSTd = aSST	;transfers metadata
aSSTd = ttest((dim_avg_n_Wrap(xx,0)),(dim_stddev_n_Wrap(xx,0))^2,nCases\				;T-test for sig.
              ,(dim_avg_n_Wrap(aoisst_s,0)),(dim_stddev_n_Wrap(aoisst_s,0))^2,nyears,True, False)

aSSTd = where(abs(aSSTd).le.p95,aSST,event_years@missing_value) ;Mask out significant values
print(aSSTd)
printVarSummary(aSSTd)

delete(X)
delete(xx)
;=============================================================================================================
;===============================================================================================================
;===============================================================================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
f0      = addfile("atl3_merra2_jra55_erai_era5_cfsr_oaflux_1984_2009.nc", "r")
  TIME   = f0->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt.and. YYYY.le.yrLast)
  qnet_atm    = f0->qnet(iYYYY,0:ndata-1)
 
  printVarSummary(qnet_atm) 
  printMinMax(qnet_atm,True)

    
  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set the time dimension
;=============================================================================================================
;=============================================================================================================
;===============================================================================================================
;=========================================================================================
;=========================================================================================
;
qnet_atm_s     = onedtond(ndtooned(qnet_atm),(/nyears,nmos,ndata/))
qnet_atm_s!0 = "time"
qnet_atm_s!1 = "month"
qnet_atm_s!2 = "lev"
printVarSummary(qnet_atm_s)
delete(qnet_atm)


;
qnet_atm_s = dtrend_n(qnet_atm_s,False,0)
;
;qnet_atm_s = dim_standardize_n(qnet_atm_s,0,0)      ; normalize 
;
std_qnet_atm_s = dim_stddev_n_Wrap(qnet_atm_s,0)
printVarSummary(std_qnet_atm_s)
;
print("========================================")
do ii=0,ndata-1,1
print(std_qnet_atm_s(:,ii))
end do
delete(ii)
print("========================================")



xx = qnet_atm_s(iz,:,:)
printVarSummary(xx)

aqnet_atm_s = dim_avg_n_Wrap(xx,0)

print("========================================")
do ii=0,ndata-1,1
print(aqnet_atm_s(:,ii))
end do
delete(ii)
print("========================================")
;===========
aqnet_atm_sd = aqnet_atm_s	;transfers metadata
aqnet_atm_sd = ttest((dim_avg_n_Wrap(xx,0)),(dim_stddev_n_Wrap(xx,0))^2,nCases\				;T-test for sig.
              ,(dim_avg_n_Wrap(qnet_atm_s,0)),(dim_stddev_n_Wrap(qnet_atm_s,0))^2,nyears,True, False)

aqnet_atm_sd = where(abs(aqnet_atm_sd).le.p95,aqnet_atm_s,aqnet_atm_s@_FillValue) ;Mask out significant values
print("====================================================================")
print("====================================================================")
do ii=0,ndata-1,1
print(aqnet_atm_sd(:,ii))
end do
delete(ii)
print("====================================================================")
print("====================================================================")
printVarSummary(aqnet_atm_sd)

delete(xx)

ndata = 5
;=============================================================================================================
;===============================================================================================================
;===============================================================================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
f0      = addfile("atl3_soda342_oras5_godas_gecco2_oaflux_1984_2009.nc", "r")
  TIME   = f0->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt.and. YYYY.le.yrLast)
  qnet_ocn    = f0->qnet(iYYYY,0:ndata-1)
 
  printVarSummary(qnet_ocn) 
  printMinMax(qnet_ocn,True)

    
  delete(f0)
  delete(TIME) 
  delete(YYYY)
  delete(iYYYY)            ; Re-set the time dimension
;=============================================================================================================
;=============================================================================================================
;===============================================================================================================
;=========================================================================================
;=========================================================================================
;
qnet_ocn_s     = onedtond(ndtooned(qnet_ocn),(/nyears,nmos,ndata/))
qnet_ocn_s!0 = "time"
qnet_ocn_s!1 = "month"
qnet_ocn_s!2 = "lev"
printVarSummary(qnet_ocn_s)
delete(qnet_ocn)


;
qnet_ocn_s = dtrend_n(qnet_ocn_s,False,0)
;
;qnet_ocn_s = dim_standardize_n(qnet_ocn_s,0,0)      ; normalize 
;
std_qnet_ocn_s = dim_stddev_n_Wrap(qnet_ocn_s,0)
printVarSummary(std_qnet_ocn_s)
;
print("========================================")
do ii=0,ndata-1,1
print(std_qnet_ocn_s(:,ii))
end do
delete(ii)
print("========================================")



xx = qnet_ocn_s(iz,:,:)
printVarSummary(xx)

aqnet_ocn_s = dim_avg_n_Wrap(xx,0)

print("========================================")
do ii=0,ndata-1,1
print(aqnet_ocn_s(:,ii))
end do
delete(ii)
print("========================================")
;===========
aqnet_ocn_sd = aqnet_ocn_s	;transfers metadata
aqnet_ocn_sd = ttest((dim_avg_n_Wrap(xx,0)),(dim_stddev_n_Wrap(xx,0))^2,nCases\				;T-test for sig.
              ,(dim_avg_n_Wrap(qnet_ocn_s,0)),(dim_stddev_n_Wrap(qnet_ocn_s,0))^2,nyears,True, False)

aqnet_ocn_sd = where(abs(aqnet_ocn_sd).le.p95,aqnet_ocn_s,aqnet_ocn_s@_FillValue) ;Mask out significant values
print("====================================================================")
print("====================================================================")
do ii=0,ndata-1,1
print(aqnet_ocn_sd(:,ii))
end do
delete(ii)
print("====================================================================")
print("====================================================================")
printVarSummary(aqnet_ocn_sd)

delete(xx)

;=============================================================================================================
;===============================================================================================================
;===============================================================================================================
print("==============================================================")
print("==============================================================")
;===============================================================================================================
;======================================
  wks = gsn_open_wks("eps","FIG_S8")         ; send graphics to PNG file
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
   r2@tiXAxisString = "Calendar month"                  ; yaxis
   r2@tiYAxisString = "~F10~Q~B~~F21~net~N~ (W m~S~-2~N~)"                  ; yaxis
  ; r2@tiYAxisString = "SST (~S~o~N~C)"                  ; yaxis
   ;r2@tiYAxisAngleF  =  270.
   r2@tmXBLabelFontHeightF = 0.013
   r2@tmYLLabelFontHeightF = 0.013
   r2@trYLog              = False                ; log scaling
   r2@trXLog              = False               ; log scaling
   r2@trYReverse           = False	            ; log scaling
   ;r2@trXMinF             = 1940.                 ; manually set lower limit
   ;r2@trXMaxF             = 2020.                ;   "          upper
  ;r2@trYMinF             =   0.                 ; manually set lower limit
  ;r2@trYMaxF             =   0.8           

   r2@trXMinF             =  0.                 ; manually set lower limit
   r2@trXMaxF             =  11.                 ;   "          upper
  ; r2@gsnYRefLineThicknesses = (/1/)
  ; r2@gsnYRefLineDashPatterns = (/0./)
   r2@xyLineColors        = xy_colors_dheating
   

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
   r2@xyDashPatterns      = (/0.,0.,0.,0.,0.,1./)           ; Dash patterns 
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
   r2@pmLegendSide           = "Top"               ; Change location of 
   r2@pmLegendParallelPosF   = 0.80                 ; move units right
   r2@pmLegendOrthogonalPosF = -0.98                ; move units down
   r2@pmLegendWidthF         = 0.10                ; Change width and
   r2@pmLegendHeightF        = 0.08                ; height of legend.
   r2@lgLabelFontHeightF     = 0.008                 ; change font height
   r2@lgTitleOn              = False               ; turn on legend title
   r2@lgTitleString          = "SST indices"           ; create legend title
   r2@lgTitleFontHeightF     = .003                ; font of legend title
   r2@lgPerimOn              = False               ; turn off box around
   r2@xyExplicitLegendLabels = (/"MERRA2","JRA55","ERAI","ERA5","CFSR","OAFLUX"/)         ; explicit labels 


;
;r2@tiXAxisString = "Calendar month"                  ; yaxis
r2@tiXAxisString = "Calendar month"                  ; yaxis


;----------------------------------------------------------------------------------------------------------------

   



  r2@trYMinF             =  -30.                 ; manually set lower limit
  r2@trYMaxF             =   10.           

  r2@pmLegendDisplayMode    = "Always"            ; turn on legend
  ;r2@xyDashPatterns      = (/1.,0.,0.,0.,0.,0./)           ; Dash patterns 
  ;r2@tiYAxisString = "Q~B~net~N~ feedback (W m~S~-2~N~)"                  ; yaxis
  ;r2@tiYAxisString = "~F18~s~F21~T/~F18~s~F21~t~N~ terms (e~S~-6~N~ K s~S~-1~N~)"                  ; yaxis
 
  plotA0  =  gsn_csm_xy(wks,xlabel,aqnet_atm_s(lev|:,month|:),r2)

  r2@xyLineColors        = xy_colors_thermocline
  r2@xyDashPatterns      = (/0.,0.,0.,0.,1.,1./)           ; Dash patterns 
  r2@xyExplicitLegendLabels = (/"SODA3","ORAS5","GODAS","GECCO2","OAFLUX","ORAS4"/)         ; oras4 is just a placeholder 

  plotA1  =  gsn_csm_xy(wks,xlabel,aqnet_ocn_s(lev|:,month|:),r2)

  r2@pmLegendDisplayMode    = "Never"            ; turn on legend
   
 

   r2@xyMarkLineMode	 =   "Markers"
   r2@xyMarkers     =  (/16.,16.,16.,16.,16.,16./)
   r2@xyMarkerColors =  r2@xyLineColors
   r2@xyMarkerSizeF   = 0.010
   r2@pmLegendDisplayMode    = "Never"            ; turn on legend
   plotB1  =  gsn_csm_xy(wks,xlabel,aqnet_ocn_sd(lev|:,month|:),r2) 
 


   r2@xyMarkerColors        = xy_colors_dheating
   plotB0  =  gsn_csm_xy(wks,xlabel,aqnet_atm_sd(lev|:,month|:),r2)

  
   


  overlay(plotA0,plotB0)
  plot(0) = plotA0
 
  overlay(plotA1,plotB1)
  plot(1) = plotA1

 
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


DATA_NAMES = (/"Atmospheric reanalyses","Ocean reanalyses","ERA5","ERAI","JRA55","MERRA2"/)
txres@txFontHeightF = 0.017   ; Use same font height as left axis
txres@txFont                = "Helvetica"
 amres@amParallelPosF   =  -0.50 ; 0.5 is the right edge of the plot, so
 amres@amOrthogonalPosF =   -0.65    ; This is the center of the plot.
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
  y_0 = -30. 
  y_1 = 30.  

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
  res_P@gsnPanelXWhiteSpacePercent 	= 0.0    ; default is 1.0
  ;res_P@gsnPanelTop                      = 0.96
  res_P@gsnPanelBottom                   = 0.20
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
