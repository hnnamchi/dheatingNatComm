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
;xy_colors_thermocline        = (/"deeppink4","greenyellow","deepskyblue","red","black","springgreen4"/)
xy_colors_thermocline        = (/"deeppink4","yellow2","deepskyblue","red","grey75","springgreen4"/)

yrStrt = 1984
yrLast = 2009

yyyymm = yyyymm_time(yrStrt,yrLast,"integer")
;print(yyyymm) 

ndata = 4
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
;all_years = (/1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,\
;		2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015/)

all_years = (/1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,\
		2005,2006,2007,2008,2009/)
print(all_years)
;
;
;event_years = (/-999,-999,1984,-999,-999,1987,1988,-999,-999,1991,-999,1993,-999,1995,1996,-999,1998,1999,-999,-999,-999,-999,-999,\
;		-999,-999,-999,2008,-999,2010,-999,-999,-999,-999,-999/)
event_years = (/1984,-999,-999,1987,1988,-999,-999,1991,-999,1993,-999,1995,1996,-999,1998,1999,-999,-999,-999,-999,-999,\
		-999,-999,-999,2008,-999/)
print(event_years)
event_years@missing_value = -999
event_years@_FillValue = event_years@missing_value


;event_years = (/1984, 1987, 1988, 1991, 1993, 1995, 1996, 1998, 1999, 2008, 2010/)
;print(event_years)


iz = ind(.not.ismissing(event_years))  ; iz will have indices of non-missing values
                               ; iz(0)=0, iz(1)=2, iz(2)=4
print(iz)

X = all_years(iz)   
print(X)   

xx = aoisst_s(iz,:)
printVarSummary(xx)

aSST = dim_avg_n_Wrap(xx,0)
print(aSST)

;===========
aSSTd = aSST	;transfers metadata
aSSTd = ttest((dim_avg_n_Wrap(xx,0)),(dim_stddev_n_Wrap(xx,0))^2,11\				;T-test for sig.
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

;f0      = addfile("/home/hyacinth/work/datasets/godas/atl3_qnet_godas_1979_2018.nc", "r")
f0      = addfile("/home/hyacinth/work/datasets/godas/new_atl3_dTdt_and_terms_0_godas_1984_2009.nc", "r")
  TIME   = f0->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt.and. YYYY.le.yrLast)
  godas    = f0->dTdt(iYYYY,0:ndata-1)
 
  printVarSummary(godas) 
  printMinMax(godas,True)

  godas = godas*1000000.
  
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
godas_s     = onedtond(ndtooned(godas),(/nyears,nmos,ndata/))
godas_s!0 = "time"
godas_s!1 = "month"
godas_s!2 = "lev"
printVarSummary(godas_s)
delete(godas)


;
godas_s = dtrend_n(godas_s,False,0)
;
;godas_s = dim_standardize_n(godas_s,0,0)      ; normalize 
;
std_godas_s = dim_stddev_n_Wrap(godas_s,0)
printVarSummary(std_godas_s)
;
print("========================================")
do ii=0,ndata-1,1
print(std_godas_s(:,ii))
end do
delete(ii)
print("========================================")



xx = godas_s(iz,:,:)
printVarSummary(xx)

agodas_s = dim_avg_n_Wrap(xx,0)

print("========================================")
do ii=0,ndata-1,1
print(agodas_s(:,ii))
end do
delete(ii)

print("========================================")
;===========
agodas_sd = agodas_s	;transfers metadata
agodas_sd = ttest((dim_avg_n_Wrap(xx,0)),(dim_stddev_n_Wrap(xx,0))^2,11\				;T-test for sig.
              ,(dim_avg_n_Wrap(godas_s,0)),(dim_stddev_n_Wrap(godas_s,0))^2,nyears,True, False)

agodas_sd = where(abs(agodas_sd).le.p95,agodas_s,agodas_s@_FillValue) ;Mask out significant values
print("====================================================================")
print("====================================================================")
do ii=0,3,1
print(agodas_sd(:,ii))
end do
delete(ii)
print("====================================================================")
print("====================================================================")
printVarSummary(agodas_sd)

delete(xx)
;==========================================================================================
;==========================================================================================
;==========================================================================================
;===============================================================================
x = new((/nmos,5/),typeof(agodas_s),agodas_s@_FillValue)
x(:,0:3) = agodas_s
x(:,4) =   agodas_s(month|:,lev|3) - (agodas_s(month|:,lev|0)+agodas_s(month|:,lev|1)+agodas_s(month|:,lev|2))		;residual
printVarSummary(x)
delete(agodas_s)
agodas_s = x
delete(x)
;-------------------------
x = new((/nmos,5/),typeof(agodas_s),agodas_s@_FillValue)
x(:,0:3) = agodas_sd
x(:,4) =  agodas_sd@_FillValue		;residual
printVarSummary(x)
delete(agodas_sd)
agodas_sd = x
delete(x)
;
;=============================================================================================================
;===============================================================================================================
;===============================================================================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
f0      = addfile("/home/hyacinth/work/datasets/gecco2/new_atl3_dTdt_and_terms_0_gecco2_1984_2009.nc", "r")
  TIME   = f0->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt.and. YYYY.le.yrLast)
  gecco2    = f0->dTdt(iYYYY,0:ndata-1)
 
  printVarSummary(gecco2) 
  printMinMax(gecco2,True)
  gecco2 = gecco2*1000000.
  
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
gecco2_s     = onedtond(ndtooned(gecco2),(/nyears,nmos,ndata/))
gecco2_s!0 = "time"
gecco2_s!1 = "month"
gecco2_s!2 = "lev"
printVarSummary(gecco2_s)
delete(gecco2)


;
gecco2_s = dtrend_n(gecco2_s,False,0)
;
;gecco2_s = dim_standardize_n(gecco2_s,0,0)      ; normalize 
;
std_gecco2_s = dim_stddev_n_Wrap(gecco2_s,0)
printVarSummary(std_gecco2_s)
;
print("========================================")
do ii=0,ndata-1,1
print(std_gecco2_s(:,ii))
end do
delete(ii)
print("========================================")



xx = gecco2_s(iz,:,:)
printVarSummary(xx)

agecco2_s = dim_avg_n_Wrap(xx,0)

print("========================================")
do ii=0,ndata-1,1
print(agecco2_s(:,ii))
end do
delete(ii)
print("========================================")
;===========
agecco2_sd = agecco2_s	;transfers metadata
agecco2_sd = ttest((dim_avg_n_Wrap(xx,0)),(dim_stddev_n_Wrap(xx,0))^2,11\				;T-test for sig.
              ,(dim_avg_n_Wrap(gecco2_s,0)),(dim_stddev_n_Wrap(gecco2_s,0))^2,nyears,True, False)

agecco2_sd = where(abs(agecco2_sd).le.p95,agecco2_s,agecco2_s@_FillValue) ;Mask out significant values
print("====================================================================")
print("====================================================================")
do ii=0,ndata-1,1
print(agecco2_sd(:,ii))
end do
delete(ii)
print("====================================================================")
print("====================================================================")
printVarSummary(agecco2_sd)

delete(xx)
;==========================================================================================
;===============================================================================
x = new((/nmos,5/),typeof(agecco2_s),agecco2_s@_FillValue)
x(:,0:3) = agecco2_s
x(:,4) =   agecco2_s(month|:,lev|3) - (agecco2_s(month|:,lev|0)+agecco2_s(month|:,lev|1)+agecco2_s(month|:,lev|2))		;residual
printVarSummary(x)
delete(agecco2_s)
agecco2_s = x
delete(x)
;-------------------------
x = new((/nmos,5/),typeof(agecco2_s),agecco2_s@_FillValue)
x(:,0:3) = agecco2_sd
x(:,4) =  agecco2_sd@_FillValue		;residual
printVarSummary(x)
delete(agecco2_sd)
agecco2_sd = x
delete(x)
;
;=============================================================================================================
;===============================================================================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
f0      = addfile("/home/hyacinth/work/datasets/oras4/new_atl3_dTdt_and_terms_0_oras4_1984_2009.nc", "r")
  TIME   = f0->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt.and. YYYY.le.yrLast)
  oras4    = f0->dTdt(iYYYY,0:ndata-1)
 
  printVarSummary(oras4) 
  printMinMax(oras4,True)

  oras4 = oras4*1000000.
  
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
oras4_s     = onedtond(ndtooned(oras4),(/nyears,nmos,ndata/))
oras4_s!0 = "time"
oras4_s!1 = "month"
oras4_s!2 = "lev"
printVarSummary(oras4_s)
delete(oras4)


;
oras4_s = dtrend_n(oras4_s,False,0)
;
;oras4_s = dim_standardize_n(oras4_s,0,0)      ; normalize 
;
std_oras4_s = dim_stddev_n_Wrap(oras4_s,0)
printVarSummary(std_oras4_s)
;
print("========================================")
do ii=0,ndata-1,1
print(std_oras4_s(:,ii))
end do
delete(ii)
print("========================================")



xx = oras4_s(iz,:,:)
printVarSummary(xx)

aoras4_s = dim_avg_n_Wrap(xx,0)

print("========================================")
do ii=0,ndata-1,1
print(aoras4_s(:,ii))
end do
delete(ii)
print("========================================")
;===========
aoras4_sd = aoras4_s	;transfers metadata
aoras4_sd = ttest((dim_avg_n_Wrap(xx,0)),(dim_stddev_n_Wrap(xx,0))^2,11\				;T-test for sig.
              ,(dim_avg_n_Wrap(oras4_s,0)),(dim_stddev_n_Wrap(oras4_s,0))^2,nyears,True, False)

aoras4_sd = where(abs(aoras4_sd).le.p95,aoras4_s,aoras4_s@_FillValue) ;Mask out significant values
print("====================================================================")
print("====================================================================")
do ii=0,ndata-1,1
print(aoras4_sd(:,ii))
end do
delete(ii)
print("====================================================================")
print("====================================================================")
printVarSummary(aoras4_sd)

delete(xx)


;==========================================================================================
;===============================================================================
x = new((/nmos,5/),typeof(aoras4_s),aoras4_s@_FillValue)
x(:,0:3) = aoras4_s
x(:,4) =   aoras4_s(month|:,lev|3) - (aoras4_s(month|:,lev|0)+aoras4_s(month|:,lev|1)+aoras4_s(month|:,lev|2))		;residual
printVarSummary(x)
delete(aoras4_s)
aoras4_s = x
delete(x)
;-------------------------
x = new((/nmos,5/),typeof(aoras4_s),aoras4_s@_FillValue)
x(:,0:3) = aoras4_sd
x(:,4) =  aoras4_sd@_FillValue		;residual
printVarSummary(x)
delete(aoras4_sd)
aoras4_sd = x
delete(x)
;
;=============================================================================================================
;===============================================================================================================
;===============================================================================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
f0      = addfile("/home/hyacinth/work/datasets/oras5/new_atl3_dTdt_and_terms_0_oras5_1984_2009.nc", "r")
  TIME   = f0->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt.and. YYYY.le.yrLast)
  oras5    = f0->dTdt(iYYYY,0:ndata-1)
 
  printVarSummary(oras5) 
  printMinMax(oras5,True)

  oras5 = oras5*1000000.
  
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
oras5_s     = onedtond(ndtooned(oras5),(/nyears,nmos,ndata/))
oras5_s!0 = "time"
oras5_s!1 = "month"
oras5_s!2 = "lev"
printVarSummary(oras5_s)
delete(oras5)


;
oras5_s = dtrend_n(oras5_s,False,0)
;
;oras5_s = dim_standardize_n(oras5_s,0,0)      ; normalize 
;
std_oras5_s = dim_stddev_n_Wrap(oras5_s,0)
printVarSummary(std_oras5_s)
;
print("========================================")
do ii=0,ndata-1,1
print(std_oras5_s(:,ii))
end do
delete(ii)
print("========================================")

xx = oras5_s(iz,:,:)
printVarSummary(xx)

aoras5_s = dim_avg_n_Wrap(xx,0)

print("========================================")
do ii=0,ndata-1,1
print(aoras5_s(:,ii))
end do
delete(ii)
print("========================================")
;===========
aoras5_sd = aoras5_s	;transfers metadata
aoras5_sd = ttest((dim_avg_n_Wrap(xx,0)),(dim_stddev_n_Wrap(xx,0))^2,11\				;T-test for sig.
              ,(dim_avg_n_Wrap(oras5_s,0)),(dim_stddev_n_Wrap(oras5_s,0))^2,nyears,True, False)

aoras5_sd = where(abs(aoras5_sd).le.p95,aoras5_s,aoras5_s@_FillValue) ;Mask out significant values
print("====================================================================")
print("====================================================================")
do ii=0,ndata-1,1
print(aoras5_sd(:,ii))
end do
delete(ii)
print("====================================================================")
print("====================================================================")
printVarSummary(aoras5_sd)

delete(xx)
;==========================================================================================
;===============================================================================
x = new((/nmos,5/),typeof(aoras5_s),aoras5_s@_FillValue)
x(:,0:3) = aoras5_s
x(:,4) =   aoras5_s(month|:,lev|3) - (aoras5_s(month|:,lev|0)+aoras5_s(month|:,lev|1)+aoras5_s(month|:,lev|2))		;residual
printVarSummary(x)
delete(aoras5_s)
aoras5_s = x
delete(x)
;-------------------------
x = new((/nmos,5/),typeof(aoras5_s),aoras5_s@_FillValue)
x(:,0:3) = aoras5_sd
x(:,4) =  aoras5_sd@_FillValue		;residual
printVarSummary(x)
delete(aoras5_sd)
aoras5_sd = x
delete(x)
;
;=============================================================================================================
;===============================================================================================================
; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
f0      = addfile("/home/hyacinth/work/datasets/soda342/atl3_dTdt_and_terms_3_soda342_1979_2018.nc", "r")
  TIME   = f0->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt.and. YYYY.le.yrLast)
  soda342    = f0->heat_storage(iYYYY,0:ndata-1)
 
  printVarSummary(soda342) 
  printMinMax(soda342,True)

  soda342 = soda342*1000000.
  
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
soda342_s     = onedtond(ndtooned(soda342),(/nyears,nmos,ndata/))
soda342_s!0 = "time"
soda342_s!1 = "month"
soda342_s!2 = "lev"
printVarSummary(soda342_s)
delete(soda342)


;
soda342_s = dtrend_n(soda342_s,False,0)
;
;soda342_s = dim_standardize_n(soda342_s,0,0)      ; normalize 
;
std_soda342_s = dim_stddev_n_Wrap(soda342_s,0)
printVarSummary(std_soda342_s)
;
print("========================================")
do ii=0,ndata-1,1
print(std_soda342_s(:,ii))
end do
delete(ii)

print("========================================")



xx = soda342_s(iz,:,:)
printVarSummary(xx)

asoda342_s = dim_avg_n_Wrap(xx,0)

print("========================================")
do ii=0,ndata-1,1
print(asoda342_s(:,ii))
end do
delete(ii)
print("========================================")
;===========
asoda342_sd = asoda342_s	;transfers metadata
asoda342_sd = ttest((dim_avg_n_Wrap(xx,0)),(dim_stddev_n_Wrap(xx,0))^2,11\				;T-test for sig.
              ,(dim_avg_n_Wrap(soda342_s,0)),(dim_stddev_n_Wrap(soda342_s,0))^2,nyears,True, False)

asoda342_sd = where(abs(asoda342_sd).le.p95,asoda342_s,asoda342_s@_FillValue) ;Mask out significant values
print("====================================================================")
print("====================================================================")
do ii=0,ndata-1,1
print(asoda342_sd(:,ii))
end do
delete(ii)
print("====================================================================")
print("====================================================================")
printVarSummary(asoda342_sd)

delete(xx)
;==========================================================================================
;===============================================================================
x = new((/nmos,5/),typeof(asoda342_s),asoda342_s@_FillValue)
x(:,0:3) = asoda342_s
x(:,4) =   asoda342_s(month|:,lev|3) - (asoda342_s(month|:,lev|0)+asoda342_s(month|:,lev|1)+asoda342_s(month|:,lev|2))		;residual
printVarSummary(x)
delete(asoda342_s)
asoda342_s = x
delete(x)
;-------------------------
x = new((/nmos,5/),typeof(asoda342_s),asoda342_s@_FillValue)
x(:,0:3) = asoda342_sd
x(:,4) =  asoda342_sd@_FillValue		;residual
printVarSummary(x)
delete(asoda342_sd)
asoda342_sd = x
delete(x)
;
;=============================================================================================================
;===============================================================================================================
;===============================================================
;===============================================================================
print("==============================================================")
print("==============================================================")
;===============================================================================================================
;======================================
  wks = gsn_open_wks("eps","FIG_5")         ; send graphics to PNG file
  plot = new(6,graphic)                ; create graphic array ; only needed if paneling
 
;-----------------------------------------------------------------------------
;=============================================================================
;===============================================================
;=============================================================== 
;***********************************************
   r2               = True                        ; plot mods desired
   r2@gsnDraw       = False                       ; do not draw
   r2@gsnFrame      = False                       ; do not advance frame
   r2@tiMainString  = ""                       ; title
   r2@tiXAxisFontHeightF = 0.012
   r2@tiYAxisFontHeightF = 0.012
   r2@tiYAxisString = "Atl3 sea surface temperature (K)"                  ; yaxis
   r2@tiXAxisString = "Calendar month"                  ; yaxis
  ; r2@tiYAxisString = "SST (~S~o~N~C)"                  ; yaxis
   ;r2@tiYAxisAngleF  =  270.
   r2@tmXBLabelFontHeightF = 0.012
   r2@tmYLLabelFontHeightF = 0.012
   r2@trYLog              = False                ; log scaling
   r2@trXLog              = False               ; log scaling
   ;r2@trYReverse           = True               ; log scaling
   ;r2@trXMinF             = 1940.                 ; manually set lower limit
   ;r2@trXMaxF             = 2020.                ;   "          upper
  r2@trYMinF             =   0.                 ; manually set lower limit
  r2@trYMaxF             =   0.8           

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
   
   r2@vpWidthF                = 0.55
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
   r2@xyDashPatterns      = (/16.,0.,0.,0.,0.,0./)           ; Dash patterns 
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
   r2@pmLegendOrthogonalPosF = -0.48                ; move units down
   r2@pmLegendWidthF         = 0.12                ; Change width and
   r2@pmLegendHeightF        = 0.10                ; height of legend.
   r2@lgLabelFontHeightF     = 0.010                 ; change font height
   r2@lgTitleOn              = False               ; turn on legend title
   r2@lgTitleString          = "SST indices"           ; create legend title
   r2@lgTitleFontHeightF     = .003                ; font of legend title
   r2@lgPerimOn              = False               ; turn off box around
   r2@xyExplicitLegendLabels = (/"~F10~Q~B~~F21~net~N~/~F8~r~F10~C~B~~F21~w~N~~F10~h","1/~F10~h~F21~[~F10~T -T~B~~F21~-h~N~].~F10~w~B~~F21~e~N~",\
				"~F10~u~F18~s~F10~T/~F18~s~F10~x~N~ + v~F18~s~F10~T/~F18~s~F10~y~N~","~F18~s~F10~T/~F18~s~F10~t~N~","~F18~n~N~"/)        


;
;r2@tiXAxisString = "Calendar month"                  ; yaxis
r2@tiXAxisString = "Calendar month"                  ; yaxis


;----------------------------------------------------------------------------------------------------------------

   plotAA  =  gsn_csm_xy(wks,xlabel,aSST,r2)



  r2@trYMinF             =  -0.2                 ; manually set lower limit
  r2@trYMaxF             =   0.2           

  r2@pmLegendDisplayMode    = "Always"            ; turn on legend
  r2@xyDashPatterns      = (/0.,0.,0.,0.,16.,16./)           ; Dash patterns 
  r2@xyLineColors        = xy_colors_thermocline
  ;r2@xyDashPatterns      = (/1.,0.,0.,0.,0.,0./)           ; Dash patterns 
  ;r2@tiYAxisString = "Q~B~net~N~ feedback (W m~S~-2~N~)"                  ; yaxis
  r2@tiYAxisString = "~F18~s~F21~T/~F18~s~F21~t~N~ terms (~F10~e~S~~F21~-6~N~ K s~S~-1~N~)"                  ; yaxis
  plotA0  =  gsn_csm_xy(wks,xlabel,agecco2_s(lev|:,month|:),r2)
 
  plotA2  =  gsn_csm_xy(wks,xlabel,aoras4_s(lev|:,month|:),r2)
  plotA3  =  gsn_csm_xy(wks,xlabel,aoras5_s(lev|:,month|:),r2)
  plotA4  =  gsn_csm_xy(wks,xlabel,asoda342_s(lev|:,month|:),r2)

   r2@trYMinF             =  -0.4                 ; manually set lower limit
   r2@trYMaxF             =   0.4           
   plotA1  =  gsn_csm_xy(wks,xlabel,agodas_s(lev|:,month|:),r2)
  r2@pmLegendDisplayMode    = "Never"            ; turn on legend


   ; r2@trYMinF             =  -0.2                 ; manually set lower limit
   ; r2@trYMaxF             =   0.4           
   ; r2@pmLegendDisplayMode    = "Always"            ; turn on legend
   ; r2@pmLegendOrthogonalPosF = -0.58                ; move units down
    ;r2@xyExplicitLegendLabels = (/"SODA3","ORAS5","ORAS4","GODAS","GECCO2"/)         ; explicit labels 
   ; r2@xyLineColors        =	(/"blue","grey60","red1","seagreen2","yellow1","black"/)
    ;r2@tiYAxisString = "Residual [~F18~n~F21~ (e~S~-6~N~ K s~S~-1~N~)]"                  ; yaxis
    ;plot(6)  =  gsn_csm_xy(wks,xlabel,residu(lev|:,month|:),r2)
    ;r2@pmLegendDisplayMode    = "Never"            ; turn on legend
    r2@tiYAxisString = ""                  ; yaxis
 
  

   delete(r2@xyLineColors)
   r2@xyMarkLineMode	 =   "Markers"
   r2@xyMarkers     =  (/16.,16.,16.,16.,16.,16./)
   r2@xyMarkerColors =  xy_colors_dheating
   r2@xyMarkerSizeF   = 0.010
   r2@pmLegendDisplayMode    = "Never"            ; turn on legend

   plotBB  =  gsn_csm_xy(wks,xlabel,aSSTd,r2)

   r2@xyMarkerColors =  xy_colors_thermocline

   plotB0  =  gsn_csm_xy(wks,xlabel,agecco2_sd(lev|:,month|:),r2)
   plotB1  =  gsn_csm_xy(wks,xlabel,agodas_sd(lev|:,month|:),r2)
   plotB2  =  gsn_csm_xy(wks,xlabel,aoras4_sd(lev|:,month|:),r2)
   plotB3  =  gsn_csm_xy(wks,xlabel,aoras5_sd(lev|:,month|:),r2)
   plotB4  =  gsn_csm_xy(wks,xlabel,asoda342_sd(lev|:,month|:),r2) 
   ;plotB5  =  gsn_csm_xy(wks,xlabel,r_merra2d(lev|:,month|:),r2)
 
   


  overlay(plotAA,plotBB)
  plot(0) = plotAA
 
  overlay(plotA0,plotB0)
  plot(1) = plotA0

  overlay(plotA1,plotB1)
  plot(2) = plotA1

  overlay(plotA2,plotB2)
  plot(3) = plotA2
  
  overlay(plotA3,plotB3)
  plot(4) = plotA3

  overlay(plotA4,plotB4)
  plot(5) = plotA4

  ;overlay(plotA5,plotB5)
  ;plot(5) = plotA5


 

;=====================================================================
;=====================================================================
nplots = 6
;-----------------------------------------------------------------------------
  txres               = True
  txres@txAngleF      = 0.          ; Rotate string clockwise 
  txres@txFontHeightF = 0.020   ; Use same font height as left axis
  txres@txFont                = "Helvetica-Bold"
 ; Move text string to center/right edge of plot.;
  amres                  = True
  amres@amParallelPosF   =  -0.67 ; 0.5 is the right edge of the plot, so
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


DATA_NAMES = (/"GECCO2","GODAS","ORAS4","ORAS5","SODA3","Residual term"/)
txres@txFontHeightF = 0.017   ; Use same font height as left axis
txres@txFont                = "Helvetica"
 amres@amParallelPosF   =  -0.47 ; 0.5 is the right edge of the plot, so
 amres@amOrthogonalPosF =   0.40    ; This is the center of the plot.
create_tt2 = gsn_create_text(wks,DATA_NAMES,txres);
tt_p1 = new(9,graphic) 
;---Attach the letters
     do i = 1,nplots-1,1
      tt_p1(i) = gsn_add_annotation(plot(i),create_tt2(i-1),amres)
     end do 
 delete(i)
;======================================================================== 
;========================================================================
  x_0 = 5.
  x_1 = 5.
  y_0 = -20. 
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
  res_P@gsnPanelXWhiteSpacePercent 	= 1.0    ; default is 1.0
  res_P@gsnPanelTop                      = 0.99
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

;-------------------------------------------------------
  ;res_P@gsnPanelCenter   = False                 ; =True .... must move labels
  ;res_P@gsnPanelRowSpec  = True                  ; tell panel what order to plot
  ;gsn_panel(wks,plot,(/1,2,2,2/),res_P)
;--------------------------------------------------------
  gsn_panel(wks,plot(0:),(/3,2/),res_P)     ; now draw as one plot
;======================================================
;======================================================
;===============================================================
;===============================================================
;***********************************************
;==============================================================
frame(wks)
end
