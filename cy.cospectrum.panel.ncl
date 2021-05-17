; latitude-phase speed spectrum
;TODO ncl -Q cy.cospectrum.panel.ncl
load "./gsn_csm.ncl"
begin
  wcStrt = systemfunc("date")

 ; FMS
  cidCtrl = "b1ctrl"
;  caseid = getenv("CASENAME")
;  if (ismissing(caseid)) then
;    print(caseid)
;    status_exit(1)
;  end if
  cidArr := (/"b2kidston","b12unifit1","b13resi"/)
  cidArr2 = (/"K11","LRF","K11-LRF"/)
  pp = 300.
  FILES1 = "ensemble-wise/cy.cospectrum."+cidArr+"."+pp+".nc"
  fdata = addfiles(FILES1,"r")
  ListSetType(fdata,"join") ;cat")
;  caseid = fdata@caseid
  dt = fdata[0]@dt
;  Ly = fdata[0]@Ly
  Lx = fdata[0]@Lx
  ds = fdata[0]@ds
  phasespeed = fdata[0]->phasespeed
  lat = fdata[0]->lat
  u_y = fdata[:]->u_y  ; repmat
  uEN_y = fdata[:]->uEN_y
  uvClim_yc = fdata[:]->uvClim_yc  ; repmat
  uvEN_yc = fdata[:]->uvEN_yc
  vvClim_yn = fdata[0]->vvClim_yn
  vvEN_yn = fdata[:]->vvEN_yn
  uuClim_yn = fdata[0]->uuClim_yn
  uuEN_yn = fdata[:]->uuEN_yn
  delete(fdata)
;  wallClockElapseTime(wcStrt, "reload nc",0)
;  Lx = 2*get_pi("f")*6371 *1000
  Lx := tofloat(Lx)

;  uvClim_yc(:,{-0.5:0.5}) = uvClim_yc@_FillValue
;  uvEN_yc(:,{-0.5:0.5}) = uvEN_yc@_FillValue

  cos2y_yc = conform(uvClim_yc, cos(lat*get_d2r("f"))^2, 1)
  uvyClim_yc = -1/cos2y_yc*center_finite_diff_n(uvClim_yc*cos2y_yc,lat*get_d2r("f")*6.371e6,False,0,1)*86400  ; Eq.9 in Kidston2011
  uvyEN_yc2 = -1/cos2y_yc*center_finite_diff_n((uvClim_yc+uvEN_yc)*cos2y_yc,lat*get_d2r("f")*6.371e6,False,0,1)*86400
  uvyEN_yc2 = linint1_n(conform(uvyClim_yc,phasespeed,2)-conform(uvyClim_yc,uEN_y-u_y,(/0,1/)), uvyEN_yc2,False,phasespeed,0,2) -uvyClim_yc

  copy_VarCoords(uvClim_yc,uvyClim_yc)
  copy_VarCoords(uvyClim_yc,uvyEN_yc2)

/;
  vvClim_yn1 = vvClim_yn(:,1:)
  vvClim_yn1 = vvClim_yn1 *conform(vvClim_yn1,ispan(1,ds(3)/2,1),1)^2 /2./get_pi("f")/6371/conform(vvClim_yn1,cos(lat*get_d2r("f")),0)
  vvEN_yn1 = vvEN_yn(:,1:) *conform(vvClim_yn1,ispan(1,ds(3)/2,1),1)^2 /2./get_pi("f")/6371/conform(vvClim_yn1,cos(lat*get_d2r("f")),0)
;  vvClim_yn1@long_name  = "v'v' at "+pp+" hPa (m~S~2~N~ s~S~-2~N~ per 1 km bin)"
  copy_VarCoords(vvClim_yn1,vvEN_yn1)
  vvEN_yn1@long_name  = "v'v' at "+pp+" hPa (m~S~2~N~ s~S~-2~N~ per 1 km bin)"
  vvnClim_y = dim_sum_n(vvClim_yn(:,1:)*conform(vvClim_yn1,ispan(1,ds(3)/2,1),1),1)/dim_sum_n(vvClim_yn(:,1:),1)
  vvnEN_y = dim_sum_n((vvClim_yn(:,1:)+vvEN_yn(:,1:))*conform(vvClim_yn1,ispan(1,ds(3)/2,1),1),1)/dim_sum_n(vvClim_yn(:,1:)+vvEN_yn(:,1:),1)
  vvn_ry = dim_sum_n(vv_ryn(:,:,1:)*conform(vv_ryn(:,:,1:),ispan(1,ds(3)/2,1) ,2),2)/dim_sum_n(vv_ryn(:,:,1:),2)
  vvnCI_y = sqrt( (dim_variance_n(vvn_ry(irClim,:),0)*(dimsizes(irClim)-1) + dim_variance_n(vvn_ry(irEN,:),0)*(dimsizes(irEN)-1))/(dimsizes(irClim)+dimsizes(irEN)-2) )*sqrt(1./dimsizes(irClim)+1./dimsizes(irEN)) *cdft_t(0.975,dimsizes(irClim)+dimsizes(irEN)-2)
;  vvnCI_y = mask(vvnCI_y,vvnCI_y.ne.0,True)
  vvnt_y = (vvnEN_y-vvnClim_y) / vvnCI_y
  vvlClim_y = cos(lat*get_d2r("f"))*dim_sum_n(vvClim_yn(:,1:)*conform(vvClim_yn1,2*get_pi("f")*6371/ispan(1,ds(3)/2,1),1),1)/dim_sum_n(vvClim_yn(:,1:),1)
  vvlEN_y = cos(lat*get_d2r("f"))*dim_sum_n((vvClim_yn(:,1:)+vvEN_yn(:,1:))*conform(vvClim_yn1,2*get_pi("f")*6371/ispan(1,ds(3)/2,1),1),1)/dim_sum_n(vvClim_yn(:,1:)+vvEN_yn(:,1:),1)
  vvl_ry = dim_sum_n(vv_ryn(:,:,1:)*conform(vv_ryn(:,:,1:),2*get_pi("f")*6371/ispan(1,ds(3)/2,1) ,2),2)/dim_sum_n(vv_ryn(:,:,1:),2)
  vvl_ry = vvl_ry*conform(vvl_ry,cos(lat*get_d2r("f")),1)
  vvlCI_y = sqrt( (dim_variance_n(vvl_ry(irClim,:),0)*(dimsizes(irClim)-1) + dim_variance_n(vvl_ry(irEN,:),0)*(dimsizes(irEN)-1))/(dimsizes(irClim)+dimsizes(irEN)-2) )*sqrt(1./dimsizes(irClim)+1./dimsizes(irEN)) *cdft_t(0.975,dimsizes(irClim)+dimsizes(irEN)-2)
;  vvlCI_y = mask(vvlCI_y,vvlCI_y.ne.0,True)
  vvlt_y = (vvlEN_y-vvlClim_y) / vvlCI_y
;  vvlClim_y = 2*get_pi("f")*6371*cos(lat*get_d2r("f")) *dim_sum_n(vvClim_yn(:,1:),1)/dim_sum_n(vvClim_yn(:,1:)*conform(vvClim_yn1,ispan(1,ds(3)/2,1),1),1)
;  vvlEN_y = 2*get_pi("f")*6371*cos(lat*get_d2r("f")) *dim_sum_n(vvClim_yn(:,1:)+vvEN_yn(:,1:),1)/dim_sum_n((vvClim_yn(:,1:)+vvEN_yn(:,1:))*conform(vvClim_yn1,ispan(1,ds(3)/2,1),1),1)
;  res@sfXArray        = conform(transpose(vvClim_yn1),lat,1)
;  res@sfYArray        = 2*get_pi("f")*6371*cos(res@sfXArray*get_d2r("f"))/conform(transpose(vvClim_yn1),ispan(1,ds(3)/2,1),0)
  uunClim_y = dim_sum_n(uuClim_yn(:,1:)*conform(vvClim_yn1,ispan(1,ds(3)/2,1),1),1)/dim_sum_n(uuClim_yn(:,1:),1)
  uunEN_y = dim_sum_n((uuClim_yn(:,1:)+uuEN_yn(:,1:))*conform(vvClim_yn1,ispan(1,ds(3)/2,1),1),1)/dim_sum_n(uuClim_yn(:,1:)+uuEN_yn(:,1:),1)
;/

  ekeClim_yn = uuClim_yn
  ekeClim_yn = (uuClim_yn+vvClim_yn)/2
;  ekeClim_yn = mask(ekeClim_yn,ekeClim_yn.ne.0,True)
  ekeEN_yn = uuEN_yn
  ekeEN_yn = (uuEN_yn+vvEN_yn)/2
;  copy_VarCoords(u_y,ekeClim_yn(:,0))
;  copy_VarCoords(ekeClim_yn,ekeEN_yn)
;  ekeEN_yn@long_name  = "u'u'+v'v' at "+pp+" hPa (m~S~2~N~ s~S~-2~N~)"
  ekeEN_yn@long_name  = ""
;  ekenClim_y = dim_sum_n(ekeClim_yn(:,1:)*conform(vvClim_yn1,ispan(1,ds(3)/2,1),1),1)/dim_sum_n(ekeClim_yn(:,1:),1)
;  ekenEN_y = dim_sum_n((ekeClim_yn(:,1:)+ekeEN_yn(:,1:))*conform(vvClim_yn1,ispan(1,ds(3)/2,1),1),1)/dim_sum_n(ekeClim_yn(:,1:)+ekeEN_yn(:,1:),1)

  type = "pdf"
; type@wkOrientation = "portrait" 
  type@wkPaperWidthF  = 190/25.4
  type@wkPaperHeightF = 2.4
  wks = gsn_open_wks(type,get_script_prefix_name+"."+pp)  ; TODO
  gsn_define_colormap(wks,"BlWhRe")  ;"BlueWhiteOrangeRed"

; Set some basic resources
  res = True

  res@cnFillOn          = True
;  res@cnFillMode        = "CellFill"
  res@cnLinesOn         = False
  res@cnLineLabelsOn  = False
  res@cnInfoLabelOn   = False
;  res@cnLevelSelectionMode = "ManualLevels"
;  res@cnMinLevelValF       = -0.08
;  res@cnMaxLevelValF       = 0.1
;  res@cnLevelSpacingF      = 0.02
;  res@cnGridBoundFillColor = "grey50"
;  res@cnMissingValFillColor= "grey50"
  res@lbOrientation      = "Vertical"
;  res@lbLabelPosition    = "Left"
;  res@pmLabelBarOrthogonalPosF = -1.5
;  res@pmLabelBarWidthF   = 0.35
  res@lbBoxSeparatorLinesOn = False
  res@lbLabelFontHeightF    = 0.04

  res@gsnAddCyclic      = False

; res@tmXBMode          = "Explicit"              ; explicit labels
;  res@tmYLMode          = "Manual"
;  res@tmYLMaxTicks      = 4
;  res@tmYLTickStartF    = -100.
;  res@tmYLTickEndF      = -0.
;  res@tmYLTickSpacingF  = 20.
;  res@mpGridLatSpacingF = 20.
; res@tmXBValues        = ispan(0,192,48)         ; location of labels 
; res@tmXBLabels        = lbl                     ; labels themselves
 ;res@tmXBLabelStride   = 24
; res@tmLabelAutoStride = True                    ; nice stride on labels
; res@tmXBFormat        = "0@!;*^se"
  res@tmXBLabelFontHeightF  = 0.04
  res@tmXBLabelDeltaF       = -0.5
  res@tmYLLabelFontHeightF  = 0.04
  res@tmYLLabelDeltaF       = -0.5
; res@xyMonoDashPattern = True
 ;res@xyDashPatterns    = 1
; res@xyMarkLineMode    = "Marklines"
; res@xyLineColors      = (/"Blue","Black"/)
 ;res@gsnXRefLine       = 0.5
 ;res@gsnXRefLineDashPattern = 1
;  res@tiMainFontHeightF = 0.04
  res@tiXAxisFontHeightF= 0.04
  res@tiYAxisFontHeightF= 0.04
  res@tiDeltaF          = 0.
;  res@tiYAxisOffsetXF   = -0.01
;  res@vpWidthF        = 0.82
;  res@vpHeightF       = 0.62
  res@gsnDraw         = False
  res@gsnFrame        = False
 ;res@gsnShape        = True
;  res@gsnMaximize = True  ; panel
;  res@gsnPaperOrientation  = "portrait"
  res@gsnXAxisIrregular2Linear = True
  res@gsnYAxisIrregular2Linear = True
;  res@trYMaxF           = 10*round(max(lat)/10,0)
 ;res@trYMinF           = 0

/;
 ;res@pmLegendDisplayMode    = "Always"            ; turn on legend
 ;res@pmLegendZone           = 1
;;res@pmLegendSide           = "Right"             ; Change location of 
 ;res@pmLegendParallelPosF   = 0.98                ; move
;;res@pmLegendOrthogonalPosF = -0.2                ; move
 ;res@lgJustification        = "BottomRight"
 ;res@pmLegendWidthF         = 0.12                ; Change width and
 ;res@pmLegendHeightF        = 0.1                 ; height of legend.
;;res@lgPerimOn              = False               ; turn off box around
;;res@lgLabelFontHeightF     = .02                 ; label font height
 ;res@xyExplicitLegendLabels = (/"0.1a","a"/)         ; create explicit labels
;/

  ores = True
;  ores@cnFillOn          = True
;  ores@cnFillMode        = "CellFill"
;  ores@cnLinesOn         = False
  ores@cnLineLabelsOn  = False
  ores@cnInfoLabelOn   = False
;  ores@cnInfoLabelFontHeightF= 0.025
;  ores@cnLevelSelectionMode = "AutomaticLevels"  ;"ManualLevels"  ;"ExplicitLevels"  ;
;  ores@cnMaxLevelCount      = 9
;  ores@cnMinLevelValF       = -0.3
;  ores@cnMaxLevelValF       = 0.3
;  ores@cnLevelSpacingF      = 0.05
;  ores@cnLevels             = (/ 0.5 /)
  ores@cnLineThicknessF     = 4.
  ores@cnMonoLineColor      = False
  ores@cnLinePalette        = "BlueWhiteOrangeRed";"BlWhRe"
;  ores@cnLineColors         = (/"blue","blue","blue","green","red","red","red"/)
;  ores@cnConstFEnableFill   = True
;  ores@cnConstFLabelOn      = False
  ores@gsnContourZeroLineThicknessF = 0.
  ores@gsnDraw         = False
  ores@gsnFrame        = False

;  ores@cnLevelSelectionMode = "ManualLevels"
;  ores@cnMinLevelValF       = -0.08
;  ores@cnMaxLevelValF       = 0.1
;  ores@cnLevelSpacingF      = 0.02
  ores@lbOrientation      = "Vertical"
;  ores@lbLabelPosition    = "Left"
;  ores@pmLabelBarOrthogonalPosF = -1.5
  ores@gsnXAxisIrregular2Linear = True
  ores@gsnYAxisIrregular2Linear = True

  gsres = True
  gsres@gsLineThicknessF = 8.
;  gsres@tfPolyDrawOrder  = "Draw"

    resP = True
;   resP@txString    = ""
 ;  resP@txPosYF  = 0.98
    resP@gsnMaximize = True
    resP@gsnPaperOrientation   = "portrait"
    resP@gsnPaperMargin        = 0.
    resP@gsnBoxMargin          = 0.3
;    resP@gsnPanelDebug = True
;   resP@gsnPanelLabelBar = True
;   resP@lbLabelFontHeightF    = 0.012
    resP@gsnPanelFigureStrings = "("+(/"a","b","c"/)+") "+cidArr2
;    resP@gsnPanelFigureStringsFontHeightF  = 0.03
    resP@gsnPanelFigureStringsPerimOn = False
    resP@gsnPanelFigureStringsJust  = "TopLeft"
;   resP@pmLabelBarOrthogonalPosF = -0.02   ; panel
    resP@gsnPanelScalePlotIndex   = 1    ; panel 24
;   resP@gsnPanelYWhiteSpacePercent = 2  ; panel 6
;    resP@gsnPanelBottom = 0.04           ; panel 7
;    resP@gsnPanelTop    = 0.935          ; panel 7
    resP@gsnPanelLeft   = 0.08
    resP@gsnPanelRight  = 0.92           ; panel 7
    resP@gsnFrame    = False

  plot = new(6, graphic)
  dum  = new(6, graphic)

  res@trYReverse      = True
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  ores@cnMonoLineColor      = True
;  ores@cnLinePalette        = "BlueWhiteOrangeRed";"BlWhRe"
  ores@cnLineThicknessF     = 1.
  ores@gsnContourNegLineDashPattern = 2
;  ores@cnLevelSelectionMode = "AutomaticLevels"
;  ores@cnMaxLevelCount      = 11
;  ores@cnLevelSpacingF = 0.1
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
;  res@tiYAxisString     = "Relative angular phase speed (m s~S~-1~N~)"
  res@tiYAxisString     = "(c-u~H-16V23~_~V-23~)/cos~F33~f~F~ (m s~S~-1~N~)"
  res@tiYAxisOffsetXF   = -0.01
  res@trYMaxF         = 20.
  res@sfYArray        = transpose(conform(uvClim_yc(0,:,:),phasespeed,1) - conform(uvClim_yc(0,:,:),u_y(0,:),0))
  res@sfXArray        = transpose(conform(uvClim_yc(0,:,:),lat,0))
 ; d/dy uv
;  uvyEN_yc2@long_name  = "~F18~-s~FB~y~N~u'v' at "+pp+" hPa (m s~S~-1~N~day~S~-1~N~ per m s~S~-1~N~ bin)"
  symMinMaxPlt(uvyEN_yc2,21,False,res)
  res@gsnXAxisIrregular2Linear = False
  res@gsnYAxisIrregular2Linear = False
  res@trGridType      = "TriangularMesh" ;"Curvilinear"
  symMinMaxPlt(uvyClim_yc,21,False,ores)
  ores@trGridType      = res@trGridType
  ores@sfXArray        = res@sfXArray
  ores@sfYArray        = res@sfYArray
do ip = 0,dimsizes(cidArr)-1
    res@tiYAxisOn            = (mod(ip,dimsizes(cidArr)).eq.0)
    res@tmYLLabelsOn         = (mod(ip,dimsizes(cidArr)).eq.0)       ; panel 7
;    res@tmYLOn               = (mod(ip,dimsizes(cidArr)).eq.0)       ; panel 7
    res@lbLabelBarOn         = (mod(ip+1,dimsizes(cidArr)).eq.0)
  plot(ip) = gsn_csm_contour(wks,transpose(uvyEN_yc2(ip,:,:)),res)
  ploto = gsn_csm_contour(wks,transpose(uvyClim_yc(ip,:,:)),ores)
  overlay(plot(ip),ploto)
end do
    resP@gsnPanelFigureStringsJust  = "TopLeft"
    resP@gsnPanelLeft   = 0.08
    resP@gsnPanelRight  = 0.92           ; panel 7
;    gsn_panel(wks,plot(:dimsizes(cidArr)-1),(/1,dimsizes(cidArr)/),resP)
  newplots = gsn_panel_return(wks,plot(:dimsizes(cidArr)-1),(/1,dimsizes(cidArr)/),resP)
;  drawNDCGrid(wks)
  frame(wks)
  getvalues ploto@contour
    "cnMinLevelValF" : wrkcmn
    "cnMaxLevelValF" : wrkcmx
    "cnLevelSpacingF" : wrkciu
  end getvalues
print("CONTOUR FROM "+wrkcmn+" TO "+wrkcmx+" BY "+wrkciu)
;print("Contour LevelSpacing="+ores@cnLevelSpacingF)
  bb_dbg = NhlGetBB(newplots)
      print("min/max x position = " + min(bb_dbg(:,2)) + "/" + max(bb_dbg(:,3)))
;      print("min/max y position = " + min(bb_dbg(:,1)) + "/" + max(bb_dbg(:,0)))
print("wkPaperHeightF="+type@wkPaperWidthF*(max(bb_dbg(:,0))-min(bb_dbg(:,1))))
  delete_VarAtts(ores,(/"trGridType","sfXArray","sfYArray"/))
  delete_VarAtts(res,(/"trGridType","sfXArray","sfYArray"/))
  delete_VarAtts(res,(/"trYReverse"/))

  res@tiXAxisString     = ""
  res@trYReverse      = True
  res@trYMaxF         = 11.
  res@tiYAxisString     = "Wavenumber"
  res@tiYAxisOffsetXF   = -0.
  ores@cnMonoLineColor      = True
;  ores@cnLineThicknessF     = 4.
  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
;  gsres0 = gsres
;  gsres0@gsLineColor      = "Black"
;  gsres0@gsLineDashPattern = 0
;  gsres@gsLineColor      = "Green4"
;  gsres@gsLineDashPattern = 2

  symMinMaxPlt(ekeEN_yn,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 9
do ip = 0,dimsizes(cidArr)-1
    res@tiYAxisOn            = (mod(ip,dimsizes(cidArr)).eq.0)
    res@tmYLLabelsOn         = (mod(ip,dimsizes(cidArr)).eq.0)       ; panel 7
;    res@tmYLOn               = (mod(ip,dimsizes(cidArr)).eq.0)       ; panel 7
    res@lbLabelBarOn         = (mod(ip+1,dimsizes(cidArr)).eq.0)
  plot(ip) = gsn_csm_contour(wks,transpose(ekeEN_yn(ip,:,:)),res)
  ploto = gsn_csm_contour(wks,transpose(ekeClim_yn),ores)
  overlay(plot(ip),ploto)
;  dum(0)  = gsn_add_polyline(wks,plot(0),lat,2*get_pi("f")*6371*cos(lat*get_d2r("f"))/ekelClim_y,gsres0)
;  dum(0)  = gsn_add_polyline(wks,plot(0),lat,ekenClim_y,gsres0)
;  dum(1)  = gsn_add_polyline(wks,plot(0),lat,2*get_pi("f")*6371*cos(lat*get_d2r("f"))/ekelEN_y,gsres)
;  dum(1)  = gsn_add_polyline(wks,plot(0),lat,ekenEN_y,gsres)
end do
    resP@gsnPanelFigureStringsJust  = "BottomLeft"
    resP@gsnPanelLeft   = 0.05
    resP@gsnPanelRight  = 0.93           ; panel 7
  newplots = gsn_panel_return(wks,plot(:dimsizes(cidArr)-1),(/1,dimsizes(cidArr)/),resP)
;  drawNDCGrid(wks)
  frame(wks)
  getvalues ploto@contour
    "cnMinLevelValF" : wrkcmn
    "cnMaxLevelValF" : wrkcmx
    "cnLevelSpacingF" : wrkciu
  end getvalues
print("CONTOUR FROM "+wrkcmn+" TO "+wrkcmx+" BY "+wrkciu)
  bb_dbg = NhlGetBB(newplots)
      print("min/max x position = " + min(bb_dbg(:,2)) + "/" + max(bb_dbg(:,3)))
;      print("min/max y position = " + min(bb_dbg(:,1)) + "/" + max(bb_dbg(:,0)))
print("wkPaperHeightF="+type@wkPaperWidthF*(max(bb_dbg(:,0))-min(bb_dbg(:,1))))

  delete(type)
  type = "png"
; type@wkOrientation = "portrait" 
  type@wkWidth        = 1280
  type@wkHeight       = 1280;640
  wks = gsn_open_wks(type,get_script_prefix_name+"."+pp)
  gsn_define_colormap(wks,"BlWhRe")  ;"BlueWhiteOrangeRed"

do ip = 0,dimsizes(cidArr)-1
    res@tiYAxisOn            = (mod(ip,dimsizes(cidArr)).eq.0)
    res@tmYLLabelsOn         = (mod(ip,dimsizes(cidArr)).eq.0)       ; panel 7
;    res@tmYLOn               = (mod(ip,dimsizes(cidArr)).eq.0)       ; panel 7
    res@lbLabelBarOn         = (mod(ip+1,dimsizes(cidArr)).eq.0)
  plot(ip) = gsn_csm_contour(wks,transpose(ekeEN_yn(ip,:,:)),res)
  ploto = gsn_csm_contour(wks,transpose(ekeClim_yn),ores)
  overlay(plot(ip),ploto)
;  dum(0)  = gsn_add_polyline(wks,plot(0),lat,2*get_pi("f")*6371*cos(lat*get_d2r("f"))/ekelClim_y,gsres0)
;  dum(0)  = gsn_add_polyline(wks,plot(0),lat,ekenClim_y,gsres0)
;  dum(1)  = gsn_add_polyline(wks,plot(0),lat,2*get_pi("f")*6371*cos(lat*get_d2r("f"))/ekelEN_y,gsres)
;  dum(1)  = gsn_add_polyline(wks,plot(0),lat,ekenEN_y,gsres)
end do
;    resP@gsnPanelFigureStringsJust  = "BottomLeft"
    resP@gsnPanelLeft   = 0.11
    resP@gsnPanelRight  = 0.88           ; panel 7
  newplots = gsn_panel_return(wks,plot(:dimsizes(cidArr)-1),(/1,dimsizes(cidArr)/),resP)
;  drawNDCGrid(wks)
  frame(wks)
  bb_dbg = NhlGetBB(newplots)
 ; 80/1280 = 0.0625--0.9375, 0.3125--0.6875,
      print("min/max x position = " + min(bb_dbg(:,2)) + "/" + max(bb_dbg(:,3)))
      print("min/max y position = " + min(bb_dbg(:,1)) + "/" + max(bb_dbg(:,0)))
  system("convert "+wks@name+".png -shave 0x320 "+wks@name+".png")

end
exit

 ;:'<,'> s?filevardef(\([^,]*\), "\([^"]*\)",.*$?filevarattdef(\1,"\2", \2)? |noh
 ;:'<,'> s?filevardef(\([^,]*\), "\([^"]*\)",.*$?\1->\2 = (/\2/)? |noh
 ;:'<,'> s?filevardef(\([^,]*\), "\([^"]*\)",.*$?\2 = fdata->\2? |noh
 ;:'<,'> g/^[^;]*contour/s/wks,\([^,]*\)/wks,transpose(\1)/gce |noh
 ;:'<,'> g/^[^;]*polyline/s/\([^,]*\),\([^,]*\),gsres/\2,\1,gsres/gce |noh

; vim: set fdm=marker foldmarker=/;,;/:

