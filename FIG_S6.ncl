;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
begin
; ==============================================================
; User defined parameters that specify region of globe and
; ==============================================================
yrStrt = 1982
yrLast = 2015
 
 
nmos = 12
  
mxlag    = 20   ;Time lag
 
xlabel = ispan(0,nmos-1,1)
print(xlabel)
  
  
p95 = 0.05

latS   = -20.
latN   =  15. 
lonL   = -60.
lonR   =  20.


p95 = 0.05

mxlag = 7
lag = 1	;first-order autocorrelation
;
;===============================================================================================================
;===============================================================================================================
;===============================================================================================================
;===============================================================================================================
;============================================================================================
; ==============================================================
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
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
 oisst_d = new(dimsizes(oisst),typeof(oisst),oisst@_FillValue)
 x   = oisst&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =oisst(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,False,False,0)
    oisst_d(tt::nmos,:,:) = work
 end do

 printVarSummary(oisst_d)
 oisst = oisst_d

;===
delete([/tt,x,x1,oisst_d,work/])
;============================================================================================
avg_OISST = wgt_areaave_Wrap(oisst(time |:, {lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0)
p_OISST = wgt_areaave_Wrap(oisst(time |:, {lon|-150:-90}, {lat | -5:5}),1.0, 1.0, 0)

; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================

; f0      = addfile ("/home/hyacinth/work/datasets/re-rain/gpcp_precip.mon.mean.nc", "r")
  f0      = addfile ("/home/hyacinth/work/datasets/rain/GPCP_precip.mon.mean.nc", "r")
  
  TIME   = f0->time    
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  gpcp    = f0->precip(iYYYY,:,:)
  printVarSummary(gpcp) 
   
   
   printMinMax(gpcp,True)
   
  ; print(xxxxxxxxxxxxxxx)
   
   
   gpcp = lonFlip(gpcp)
   

   delete(f0)
   delete(TIME) 
   delete(YYYY)
   delete(iYYYY)            ; Re-set the time dimension
;=============================================================================================================
;=============================================================================================================
 ;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
 gpcp_d = new(dimsizes(gpcp),typeof(gpcp),gpcp@_FillValue)
 x   = gpcp&time
 x1  = x(::12)
 do tt = 0,nmos-1
    work =gpcp(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,False,False,0)
    gpcp_d(tt::nmos,:,:) = work
 end do

 printVarSummary(gpcp_d)
 gpcp = gpcp_d

;===
delete([/tt,x,x1,gpcp_d,work/])
;---------------------------------------------------------------------------
;=========================================================================================
;=========================================================================================
avg_GPCP = wgt_areaave_Wrap(gpcp(time |:, {lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0)
p_GPCP = wgt_areaave_Wrap(gpcp(time |:, {lon|-150:-90}, {lat | -5:5}),1.0, 1.0, 0)
;---------------------------------------------------------------------------
;=========================================================================================
;=========================================================================================
;Climatological-mean
dimZ = dimsizes(avg_GPCP)
printVarSummary(dimZ)
;
nyears = dimZ(0)/nmos
printVarSummary(nyears)
;
;---------------------------
a_GPCP4d     = onedtond(ndtooned(avg_GPCP(time|:)),(/nyears,nmos/))
printVarSummary(a_GPCP4d)
;
p_GPCP4d     = onedtond(ndtooned(p_GPCP(time|:)),(/nyears,nmos/))
printVarSummary(p_GPCP4d)
;
;---------------------------
z4d     = onedtond(ndtooned(avg_OISST(time|:)),(/nyears,nmos/))
printVarSummary(z4d)



;============================================================================================
;---------------------------

;---------------------------
;Anomalies
a_GPCP = wgt_areaave_Wrap(gpcp(time |:, {lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0)
a_OISST = wgt_areaave_Wrap(oisst(time |:, {lon|-20:0}, {lat | -3:3}),1.0, 1.0, 0)
;
p_GPCP = wgt_areaave_Wrap(gpcp(time |:, {lon|-150:-90}, {lat | -5:5}),1.0, 1.0, 0)
p_OISST = wgt_areaave_Wrap(oisst(time |:, {lon|-150:-90}, {lat | -5:5}),1.0, 1.0, 0)
;-----------------------
;
;---------------------------
a_GPCP4d     = onedtond(ndtooned(a_GPCP(time|:)),(/nyears,nmos/))
a_GPCP4d!0 = "time"
a_GPCP4d!1 = "month"
printVarSummary(a_GPCP4d)
;
a_OISST4d     = onedtond(ndtooned(a_OISST(time|:)),(/nyears,nmos/))
a_OISST4d!0 = "time"
a_OISST4d!1 = "month"
printVarSummary(a_OISST4d)
;
p_OISST4d     = onedtond(ndtooned(p_OISST(time|:)),(/nyears,nmos/))
p_OISST4d!0 = "time"
p_OISST4d!1 = "month"
printVarSummary(p_OISST4d)
;
p_GPCP4d     = onedtond(ndtooned(p_GPCP(time|:)),(/nyears,nmos/))
p_GPCP4d!0 = "time"
p_GPCP4d!1 = "month"
printVarSummary(p_GPCP4d)

;print(xxxxxxxxxxxxxxxxxx)

;========================================================================================================================
;Normalize the SSTs
a_OISST4d = dim_standardize_n_Wrap(a_OISST4d,0,0)  ; Or, overwrite with deviations
p_OISST4d = dim_standardize_n_Wrap(p_OISST4d,0,0)  ; Or, overwrite with deviations
;
;
r_a_months =  new((/nmos,nmos/),float)
printVarSummary(r_a_months)
r_p_months     = r_a_months 


do ii = 0,nmos-1,1
r_a_months(ii,:) = regCoef_n(a_OISST4d(:,ii),a_GPCP4d(:,:),0,0)	;y (latitude) is plotted before x (latitude) in ncl
r_p_months(ii,:) = regCoef_n(p_OISST4d(:,ii),p_GPCP4d(:,:),0,0)
end do
delete(ii)
;print(r_a_months)
;
printVarSummary(r_a_months)
printVarSummary(r_p_months)
;
;
;
d_a_months = r_a_months
d_p_months = r_p_months
do ii = 0,nmos-1,1
d_a_months(ii,:) = escorc_n(a_OISST4d(:,ii),a_GPCP4d(:,:),0,0)	;y (latitude) is plotted before x (latitude) in ncl
d_p_months(ii,:) = escorc_n(p_OISST4d(:,ii),p_GPCP4d(:,:),0,0)
end do
delete(ii)
;
;print(d_a_months)
printVarSummary(d_a_months)
printVarSummary(d_p_months)
;
SS = 0.34
d_a_months = where(abs(d_a_months).le.SS,r_a_months@_FillValue,r_a_months) ; Significant TEMP
d_p_months = where(abs(d_p_months).le.SS,r_p_months@_FillValue,r_p_months) ; Significant TEMP
;print(xxxxxxxxxxxxx)
;===============================================================================
;===============================================================================
;===============================================================================
;   tval = onedtond(r_a_months@tval,dimsizes(r_a_months))
;   df   = onedtond(r_a_months@nptxy,dimsizes(r_a_months)) - 2
;   b = tval    ; b must be same size as tval (and df)
;   b = 0.5
;   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon) 
;   r_a_months@long_name   = "regression coefficient"
;   r_a_months@long_name = "probability"
;   r_a_monthsd = r_a_months
;   r_a_monthsd = where(prob.le.p95,r_a_monthsd@_FillValue,r_a_monthsd) ; Significant TEMP
;   printVarSummary(r_a_monthsd)
;   printMinMax(r_a_monthsd,True)
;=========================================================================================
print("===========================================================")
;-----------------------------------------------------------------------------------------
print("===========================================================")
;=========================================================================================
;=========================================================================================
print("===========================================================")
 wks = gsn_open_wks("eps","FIG_S6")
 plot = new(2,graphic)                ; create graphic array

  res                	   = True          ; plot mods desired
  res@gsnDraw      	   = False                       ; do not draw
  res@gsnFrame     	   = False                       ; do not advance frame
  res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels  
  res@cnMinLevelValF       = -2.            ; set min contour level
  res@cnMaxLevelValF       =  2.             ; set max contour level
  res@cnLevelSpacingF      =  .2            ; set contour spacing   
  res@cnFillOn             = True          ; turn on color
  res@cnLinesOn            = False          ; turn on color
  res@cnLineLabelsOn	   = False
  res@cnInfoLabelOn	   = False
  ;res@cnFillPalette        = "gui_default" ; set color map
  res@tmXBMode = "Explicit"
  res@tmXBValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
  res@tmXBLabels = (/"Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
  res@tmYLMode = "Explicit"
  res@tmYLValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
  res@tmYLLabels = (/"Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
  res@vpWidthF            = 0.8           ; change aspect ratio of plot
  res@vpHeightF           = 0.5

  res@tiXAxisString = "SST"    
  res@tiYAxisString = "Precipitation"
  res@tmXBLabelFontHeightF = 0.020
   
  res@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  ;res@gsnMaximize         = True          ; maximize plot size

  res@lbLabelBarOn         =  True       ; turn off individual lb's
  res@lbOrientation        = "vertical"
  res@lbTitleFontHeightF   = .02
  res@lbLabelFontHeightF   = 0.02
  res@lbLeftMarginF         = 0.10
  res@lbRightMarginF        = 0.55
  res@lbLabelFont   	    = "Helvetica-Bold"     ; label font
  ;res@lbTitleString        = "SST trend (~S~o~N~C per 38-yr)"
  res@lbTitleString        = "Precipitation (mm day~S~-1~N~)"
  res@lbTitleAngleF        = 270.
  res@lbTitlePosition      = "Right"
  res@lbTitleDirection     = "Across"
  
  res@lbLabelStride        = 10
  res@gsnLeftString  = ""  
  res@gsnLeftStringOrthogonalPosF	= -0.15
  res@gsnLeftStringParallelPosF		=  0.08
;-----------------------------------------------------------------------------
  ;-----------------------------------------------------------------------------
  res1                      =   True
  res1@gsnDraw              = False        ; don't draw yet
  res1@gsnFrame             = False        ; don't advance frame yet
  res1@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res1@cnMinLevelValF       =  -1.            ; set min contour level
  res1@cnMaxLevelValF       =   1.             ; set max contour level
  res1@cnLevelSpacingF      =  .01            ; set contour spacing  
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

;plotB0 = gsn_csm_contour(wks,d_OISSTd,res1) 
 ; plotB0 = gsn_contour_shade(plotB0,d_OISSTd@_FillValue,0.90,opt)

  res@gsnLeftString  = "Atl3 region"                       ; title
  ;res@gsnCenterString  = "Atmospheric forcing"  
  plotA0 = gsn_csm_contour(wks,r_a_months,res)       ; contour the variable
  plotB0 = gsn_csm_contour(wks,d_a_months,res1)
  plotB0 = gsn_contour_shade(plotB0,d_a_months@_FillValue,p95,opt)
  
  overlay(plotA0,plotB0)
  plot(0) = plotA0

  ;res@lbLabelBarOn         =  False       ; turn off individual lb's
  ;delete(res@tiXAxisString)
  ;res@tiXAxisString = "Z~B~23~N~"
  res@gsnLeftString  = "Nino3 region"                       ; title
  res@tmYRBorderOn      = False	;remove the border
  plotA3 = gsn_csm_contour(wks,r_p_months,res)       ; contour the variable
  plotB3 = gsn_csm_contour(wks,d_p_months,res1)       ; contour the variable
  plotB3 = gsn_contour_shade(plotB3,d_p_months@_FillValue,p95,opt)
  overlay(plotA3,plotB3)
  plot(1) = plotA3
  
;========================================================================
;========================================================================
;=====================================================
  txres5               = True
  txres5@txAngleF      = 0.          ; Rotate string clockwise 
  txres5@txFontHeightF = 0.030   ; Use same font height as left axis
  txres5@txFont                = "Helvetica-Bold"
 ; txid = gsn_create_text(wks, "SST anomaly (~S~o~N~C)", txres)
   txid5 = gsn_create_text(wks, "a", txres5);
; Move text string to center/right edge of plot.;
  amres5                  = True
  amres5@amParallelPosF   =  -0.65 ; 0.5 is the right edge of the plot, so
                                  ; 0.6 is a little further to the right.
  amres5@amOrthogonalPosF =  -0.50    ; This is the center of the plot.
  amres5@amJust           = "TopLeft"
  annoid5 = gsn_add_annotation(plot(0),txid5,amres5)  ; Attach string to plot


; amres5@amParallelPosF   =  -0.65 ; 0.5 is the right edge of the plot, so

; txres5@txFontHeightF = 0.025   ; Use same font height as left axis
 txid5b = gsn_create_text(wks, "b", txres5);
  annoid5b = gsn_add_annotation(plot(1),txid5b,amres5)  ;

 ;txres5@txFontHeightF = 0.022   ; Use same font height as left axis
 txid5c = gsn_create_text(wks, "c", txres5);
; annoid5c = gsn_add_annotation(plot(2),txid5c,amres5)  ;
 
 txid5d = gsn_create_text(wks, "d", txres5);
; annoid5d = gsn_add_annotation(plot(3),txid5d,amres5)  ; 
   
;======================================================== 
;========================================================
;========================================================
;===============================================================
;==============================================================
;======================================================== 
;========================================================
txr =  True
txr@txAngleF       =  0.
txr@txFontHeightF  =  0.020
txr@txFontColor  = "red"
txr@txFont    =  "Helvetica"

amres        =   True
amres@amParallelPosF    = -0.20    ;left/right pos
amres@amOrthogonalPosF  =  0.29    ;up/down pos
amres@amJust            =  "CenterRight"

txid0 = gsn_create_text(wks,"~F10~d~F21~SST/dx~F18~oq~F10~d~F21~Prec/dx",txr)
;txid1 = gsn_create_text(wks,"~F0~Peak~F0~=~F0~12-13 yr~S~-1~N~",txr)
;anid0 =  gsn_add_annotation(plot(0),txid0,amres)
;
txr@txFontColor  = "springgreen4"
txid1 = gsn_create_text(wks,"~F10~d~F21~Prec/dx~F18~oq~F21~d~F33~t~F21~/dx",txr)
amres@amOrthogonalPosF  =  0.36    ;up/down pos
;anid1 =  gsn_add_annotation(plot(0),txid1,amres)
;
txr@txFontColor  = "springgreen3"
txid2 = gsn_create_text(wks,"~F10~d~F21~SST/dx~F18~oq~F21~d~F33~t~F21~/dx",txr)
amres@amOrthogonalPosF  =  0.43    ;up/down pos
amres@amParallelPosF  =   -0.09    ;left/right pos
;anid2 =  gsn_add_annotation(plot(0),txid2,amres)
;===============================================================
;==============================================================
;-----------------------------------------------------------------------------
;========================================================
;===============================================================
;===============================================================
;---Resources for the polyline
  pres0                  = True           ; polyline mods desired  
  pres0@gsFillColor      = "grey85"         ; color of lines
  pres0@tfPolyDrawOrder       = "PreDraw"    ; Draw any attached primitives
  pres0@gsFillIndex = 0
  pres0@gsFillOpacityF     = 0.40

;-------------------------
;===============================================================
  latEQ = (/ 12.,  12.,  1.,   1., 12./)
  lonEQ = (/12., 1., 1., 12.,12./)
 
  x0 = 0.
  x1 = 11.

  y0 = 0.
  y1 = 11.

;---Resources for the polyline
  pres0                  = True           ; polyline mods desired
  pres0@gsLineThicknessF = 0.5            ; line thickness   
  pres0@gsLineColor      = "black"         ; color of lines
  pres0@gsLineDashPattern = 0
  pres0@tfPolyDrawOrder   = "postDraw"
 dum2 = gsn_add_polyline(wks,plot(0), (/x0,x1/),(/y0,y1/) ,pres0)
 dum3 = gsn_add_polyline(wks,plot(1), (/x0,x1/),(/y0,y1/) ,pres0)


   x_0 = 5.
  x_1 = 11.
  y_0 = -10. 
  y_1 = 40.  
;pres0@gsLineDashPattern = 0
; pres0@gsLineColor      = "grey"         ; color of lines
 pres0@gsLineThicknessF = 3.5            ; line thickness
 pres0@gsLineDashPattern = 16   
 dumL0 = gsn_add_polyline(wks,plot(0), (/x_0,x_0/),(/y_0,y_1/) ,pres0)
 pres0@gsLineThicknessF = 7.0            ; line thickness
 dumL1 = gsn_add_polyline(wks,plot(1), (/x_1,x_1/),(/y_0,y_1/) ,pres0)
;===============================================================
;===============================================================
;===============================================================
; panel plot only resources
  resP                   	  	= True        ; modify the panel plot
  resP@gsnDraw            		 = True       ; don't draw yet
  resP@gsnFrame            		 = False        ; don't advance frame yet
  resP@gsnPanelLabelBar    		= False        ; add common colorbar
  resP@lbLabelAutoStride  		= True         ; auto stride on labels
  resP@gsnPanelYWhiteSpacePercent 	= 1.0    ; default is 1.0
  resP@gsnPanelXWhiteSpacePercent 	= 1.0    ; default is 1.0
  resP@gsnMaximize         		= True         ; large format
  resP@lbTitleFontHeightF		= .010 
  resP@lbLabelFontHeightF               = 0.010
  resP@lbLabelFont   			= "Helvetica-Bold"     ; label font
  resP@lbTitleString    		= "SST anomaly (K)"
  resP@gsnPaperOrientation 		= "portrait"
  resP@gsnPanelFigureStringsBackgroundFillColor = "transparent"
  resP@gsnPanelFigureStringsFontHeightF = 0.010
  resP@gsnPanelFigureStrings = (/"","","","","e","f"/) ; add strings to panel
  resP@amJust   			= "TopLeft"
  resP@gsnPanelTop                      = 0.99
  resP@gsnPanelBottom                   = 0.05
  resP@gsnPanelLeft                     = 0.05
  resP@gsnPanelRight                    = 0.95
  resP@gsnPanelFigureStringsPerimOn 	= False
  resP@gsnMaximize         		= True         ; large format
  ;gsn_panel(wks,plot(0:1),(/1,2/),resP)
  gsn_panel(wks,plot,(/2,1/),resP)
;-------------------------------------------------------
   res_P             			  = True                        ; panel mods desired
   res_P@gsnDraw              		  = True       ; don't draw yet
   res_P@gsnFrame         		  = False        ; don't advance frame yet
   res_P@gsnPanelYWhiteSpacePercent       = 0.0    ; default is 1.0
   res_P@gsnPanelXWhiteSpacePercent       = 2.0    ; default is 1.0
   res_P@gsnPanelYWhiteSpacePercent       = 0.0    ; default is 1.0
   res_P@gsnPanelXWhiteSpacePercent       = 1.0    ; default is 1.0
   res_P@gsnPanelTop                      = 0.75
   res_P@gsnPanelBottom                   = 0.05
   res_P@gsnPanelLeft                     = 0.05
   res_P@gsnPanelRight                    = 0.95
   res_P@gsnPanelFigureStringsPerimOn 	  = True 
   res_P@gsnPaperOrientation 		  = "portrait"
   res_P@gsnPanelDebug      		  = True
   res_P@amJust   			  = "TopLeft"
   res_P@gsnPanelFigureStringsPerimOn     = False
   res_P@gsnPanelFigureStringsFontHeightF = 0.015
   res_P@gsnPanelFigureStringsBackgroundFillColor = "grey"
   res_P@gsnPanelFigureStrings= (/"","","c","d","e","f","g","h","i","j","k","l"/) ; add strings to panel
   res_P@gsnMaximize = True                        ; blow up plot
 ; gsn_panel(wks,plot(2:3),(/1,2/),res_P)               ; create panel plots
;-------------------------------------------------------
;--------------------------
;===============================================================
;***********************************************
;=========================================================================================
;-----------------------------------------------------------------------------------------
frame(wks)

end
