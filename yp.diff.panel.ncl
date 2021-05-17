; plot difference in climatology
;TODO ncl -Q yp.diff.panel.ncl
load "./gsn_csm.ncl"

begin
  cidCtrl = "b1ctrl"
;  cidArr = getenv("CASENAME")
;  pdfname = "ensemble-wise/"+get_script_prefix_name+"."+getenv("CASENAME")
;  if (ismissing(cidArr)) then
    cidArr := (/"b2kidston","b12unifit1","b13resi"/)
    cidArr2 = (/"K11","LRF","K11-LRF"/)
    pdfname = get_script_prefix_name
;    print(cidArr)
;    status_exit(1)
;  end if

  type = "pdf"
; type@wkWidth        = 1200
; type@wkHeight       = 1000
  type@wkPaperWidthF  = 190/25.4
  type@wkPaperHeightF = 4.9
  wks = gsn_open_wks(type, pdfname)
; gsn_define_colormap(wks,"hotres")
  res = True
  res@tiYAxisString = "Pressure (hPa)"
  res@tiYAxisFontHeightF = 0.04
  res@tiDeltaF           = 0.
  res@tmXBLabelFontHeightF  = 0.04
  res@tmXBLabelDeltaF       = -0.5
  res@tmYLLabelFontHeightF  = 0.04
  res@tmYLLabelDeltaF       = -0.5
  res@gsnRightString  = ""
  res@gsnStringFontHeightF  = 0.04

  res@cnFillOn        = True
; res@cnFillMode      = "RasterFill" ;"AreaFill" ;
  res@cnFillPalette   = "BlWhRe"
  res@cnLinesOn       = False
  res@cnLineLabelsOn  = False
  res@cnInfoLabelOn   = False
;  res@cnLineLabelFormat = "0@.2f";"*+^sg";"*.2f";
 ;res@cnLineLabelFontHeightF = 0.01
 ;res@cnLevelSelectionMode = "EqualSpacedLevels""ExplicitLevels"  ;
 ;res@cnMaxLevelCount      = 99
; res@cnLevelSelectionMode = "ManualLevels"  ;
; res@cnMinLevelValF       = -3*50
; res@cnMaxLevelValF       = 47*50
; res@cnLevelSpacingF      = 0.2*50
; res@cnLevels             = (/ 0.5 /)
; res@cnFillColors         = (/"White","Green"/)
;  res@lbLabelBarOn = False
;  res@lbAutoManage = False
 ;res@lbPerimOn    = False
 ;res@lbPerimDashSegLenF = 0.
  res@lbOrientation      = "Vertical"
;  res@pmLabelBarOrthogonalPosF = 0.1
;  res@pmLabelBarWidthF   = 0.25
;  res@pmLabelBarSide     = "Right"
;  res@lbJustification    = "CenterLeft"
; res@lbLabelStride      = 1
 ;res@lbLabelAutoStride  = False
 ;res@lbBoxSizing        = "ExplicitSizing"
 ;res@lbBoxFractions     = array_append_record(fspan(0,0.30,6),fspan(0.36,1,19),0)
  res@lbBoxSeparatorLinesOn = False
  res@lbLabelFontHeightF = 0.04
  res@lbTitleOn          = True
  res@lbTitleFontHeightF = 0.04
  res@lbTitlePosition    = "Right"
  res@lbTitleDirection   = "Across"
  res@lbTitleAngleF      = 90.
  res@lbTitleOffsetF     = 0.2

;  res@vpXF            = 0.26
;  res@vpWidthF        = 0.82
;  res@vpYF            = 0.51
;  res@vpHeightF       = 0.62
;  res@gsnMaximize     = True  ; panel
;  res@gsnPaperOrientation   = "portrait"
;  res@gsnPaperMargin  = 0.
;  res@gsnBoxMargin    = 0.
  res@gsnFrame        = False        ; Don't draw or advance
  res@gsnDraw         = False        ; frame yet.
; res@pmLegendDisplayMode    = "Always"

; res@tiMainString  = ""
  res@gsnXAxisIrregular2Linear = True
  res@gsnYAxisIrregular2Linear = True
  res@trYReverse    = True
  res@trYMinF       = 0.
  res@trYMaxF       = 1000.
  gsres = True
  gsres@gsLineThicknessF = 3.
  ores = True
;  ores@cnLineLabelsOn  = False
  ores@cnLineLabelFontHeightF = 0.03
  ores@cnInfoLabelOn   = False
;  ores@cnInfoLabelFontHeightF= 0.025
  ores@cnLineThicknessF= 2.
;  ores@cnLevelSelectionMode = "ManualLevels"
;  ores@cnMinLevelValF  = 0.0001
;  ores@cnMaxLevelValF  = 0.001
;  ores@cnLevelSpacingF = 0.0001
  ores@gsnLeftString  = ""
  ores@gsnRightString = ""
  ores@gsnDraw         = False
  ores@gsnFrame        = False
  ores@gsnXAxisIrregular2Linear = True
  ores@gsnYAxisIrregular2Linear = True
  ores@trYReverse    = True
  ores@trYMinF       = 0.
  ores@trYMaxF       = 1000.

    resP = True
;   resP@txString    = ""
 ;  resP@txPosYF  = 0.98
    resP@gsnMaximize = True
    resP@gsnPaperOrientation   = "portrait"
    resP@gsnPaperMargin        = 0.
    resP@gsnBoxMargin          = 0.
;    resP@gsnPanelDebug = True
;    resP@gsnPanelLabelBar = True
;    resP@lbOrientation      = "Vertical"
;    resP@lbBoxSeparatorLinesOn = False
;   resP@lbLabelFontHeightF    = 0.012
;    resP@gsnPanelFigureStrings = (/"N increase"/)
    resP@gsnPanelFigureStringsFontHeightF  = 0.03
    resP@gsnPanelFigureStringsJust  = "TopLeft"
;   resP@pmLabelBarOrthogonalPosF = -0.02   ; panel
    resP@gsnPanelScalePlotIndex   = 1+dimsizes(cidArr)    ; panel 24
    resP@gsnPanelYWhiteSpacePercent = 2  ; panel 6
;    resP@gsnPanelBottom = 0.04           ; panel 7
;    resP@gsnPanelTop    = 0.935          ; panel 7
    resP@gsnPanelLeft   = 0.08
    resP@gsnPanelRight  = 0.89           ; panel 7
    resP@gsnFrame    = False

  plot = new(2*dimsizes(cidArr), graphic)
  dum  = new(2*dimsizes(cidArr), graphic)

  lev_p = 25.*ispan(1,39,1)
  f0  = addfile("ensemble-wise/zyorg_"+cidCtrl+".nc","r")
  ds_org = getfilevardimsizes(f0,"ucomp")
  bk   = f0->bk
  bk_f = f0->pfull / f0->phalf(dimsizes(bk)-1)  ; end
;  ps = f0->ps
;  t0 = vinth2p(f0->temp,bk_f*0,bk_f,lev_p,ps,1,1000.,1,True) ; 1 means linear interpolation
  p_fy = conform_dims(ds_org,bk_f,0)*conform_dims(ds_org,f0->ps/100.,1)
  t0 = linint1_n_Wrap(p_fy,f0->temp,False,lev_p,0,0)
  u0 = linint1_n_Wrap(p_fy,f0->ucomp,False,lev_p,0,0)
  ds = dimsizes(u0)
  u02 = u0(:,0:ds(1)/2-1)
  u02 = (u02 +u0(:,ds(1)-1:ds(1)/2))/2
  t0@units = "degK"  ; orig. deg_k
  s0 = static_stability(lev_p*100,t0,0,0)
  ds = dimsizes(s0)
;  s02 = s0(:,0:ds(1)/2-1)
;  s02 = (s02 +s0(:,ds(1)-1:ds(1)/2))/2
  t0_fy = f0->temp(:,0:ds(1)/2-1)
  t0_fy = (t0_fy +f0->temp(:,ds(1)-1:ds(1)/2))/2
  pTrop0 = trop_wmo(transpose(p_fy(:,0:ds(1)/2-1)+p_fy(:,ds(1)-1:ds(1)/2))/2,transpose(t0_fy),0,False)
  delete(f0)

  f1  = addfile("ensemble-wise/zyorg_"+cidArr(0)+".nc","r")
  p_fy = conform_dims(ds_org,bk_f,0)*conform_dims(ds_org,f1->ps/100.,1)
  t = linint1_n_Wrap(p_fy,f1->temp,False,lev_p,0,0)
  t = t - t0
  symMinMaxPlt(t,31,False,res)
  resu = res
  f1  = addfile("ensemble-wise/zyorg_"+cidArr(2)+".nc","r")
  p_fy = conform_dims(ds_org,bk_f,0)*conform_dims(ds_org,f1->ps/100.,1)
  u = linint1_n_Wrap(p_fy,f1->ucomp,False,lev_p,0,0)
  u = u - u0
  symMinMaxPlt(u,21,False,resu)
do ip = 0,dimsizes(cidArr)-1
  f1  = addfile("ensemble-wise/zyorg_"+cidArr(ip)+".nc","r")
  p_fy = conform_dims(ds_org,bk_f,0)*conform_dims(ds_org,f1->ps/100.,1)
  t = linint1_n_Wrap(p_fy,f1->temp,False,lev_p,0,0)
  t = t - t0
  ds = dimsizes(t)
  t2 = t(:,0:ds(1)/2-1)
  t2 = (t2 +t(:,ds(1)-1:ds(1)/2))/2
  t2@long_name = "("+tochar(toint(tochar("a"))+ip)+") "+cidArr2(ip)
;  symMinMaxPlt(t2,21,False,res)
;  res@gsnDraw         = False
;  res@gsnFrame        = False
  res@lbTitleString      = "~F5~D~F~T (K)"
    res@tiYAxisOn            = (mod(ip,dimsizes(cidArr)).eq.0)
    res@tmYLLabelsOn         = (mod(ip,dimsizes(cidArr)).eq.0)       ; panel 7
;    res@tmYLOn               = (mod(ip,dimsizes(cidArr)).eq.0)       ; panel 7
    res@lbLabelBarOn         = (mod(ip+1,dimsizes(cidArr)).eq.0)
  plot(ip) = gsn_csm_contour(wks, t2, res)
  gsres@gsLineColor      = "grey45";"Black"
  gsres@gsLineDashPattern = 0
;  dum(ip*2)  = gsn_add_polyline(wks,plot(ip),t2&lat,pTrop0,gsres)
  t_fy = f1->temp(:,0:ds(1)/2-1)
  t_fy = (t_fy +f1->temp(:,ds(1)-1:ds(1)/2))/2
  pTrop = trop_wmo(transpose(p_fy(:,0:ds(1)/2-1)+p_fy(:,ds(1)-1:ds(1)/2))/2,transpose(t_fy),0,False)
  gsres@gsLineColor      = "Black";"Green4"
  gsres@gsLineDashPattern = 2
;  dum(ip*2+1)  = gsn_add_polyline(wks,plot(ip),t2&lat,pTrop,gsres)
;  draw(plot(ip))
;  frame(wks)
  delete([/t2,t/])
  u = linint1_n_Wrap(p_fy,f1->ucomp,False,lev_p,0,0)
  u = u - u0
  ds = dimsizes(u)
  u2 = u(:,0:ds(1)/2-1)
  u2 = (u2 +u(:,ds(1)-1:ds(1)/2))/2
  u2@long_name = "("+tochar(toint(tochar("a"))+ip+dimsizes(cidArr))+") "+cidArr2(ip)
  print(cidArr(ip)+"	"+max(abs(u2({200:},:))))
;  symMinMaxPlt(u2,21,False,resu)
  resu@lbTitleString      = "~F5~D~F~u (m s~S~-1~N~)"
    resu@tiYAxisOn            = (mod(ip,dimsizes(cidArr)).eq.0)
    resu@tmYLLabelsOn         = (mod(ip,dimsizes(cidArr)).eq.0)       ; panel 7
;    resu@tmYLOn               = (mod(ip,dimsizes(cidArr)).eq.0)       ; panel 7
    resu@lbLabelBarOn         = (mod(ip+1,dimsizes(cidArr)).eq.0)
  plot(ip+dimsizes(cidArr)) = gsn_csm_contour(wks, u2, resu)
  ploto = gsn_csm_contour(wks, u02, ores)
  overlay(plot(ip+dimsizes(cidArr)),ploto)
;  draw(plot(ip+dimsizes(cidArr)))
;  frame(wks)
  delete([/u2,u,f1/])
end do
 ; gsnPanelBottom causes problem
;    resP@gsnPanelBottom = 0.5            ; panel 7
;    gsn_panel(wks,plot(:dimsizes(cidArr)-1),(/1,dimsizes(cidArr)/),resP)
;    gsn_panel(wks,plot,(/2,dimsizes(cidArr)/),resP)
  newplots = gsn_panel_return(wks,plot,(/2,dimsizes(cidArr)/),resP)
;  drawNDCGrid(wks)
  frame(wks)
;  getvalues ploto@contour
;    "cnMinLevelValF" : wrkcmn
;    "cnMaxLevelValF" : wrkcmx
;    "cnLevelSpacingF" : wrkciu
;  end getvalues
;print("CONTOUR FROM "+wrkcmn+" TO "+wrkcmx+" BY "+wrkciu)
  bb_dbg = NhlGetBB(newplots)
      print("min/max x position = " + min(bb_dbg(:,2)) + "/" + max(bb_dbg(:,3)))
;      print("min/max y position = " + min(bb_dbg(:,1)) + "/" + max(bb_dbg(:,0)))
print("wkPaperHeightF="+type@wkPaperWidthF*(max(bb_dbg(:,0))-min(bb_dbg(:,1))))
exit

 ; uuuuuu
  f1  = addfile("ensemble-wise/zyorg_"+cidArr(2)+".nc","r")
  p_fy = conform_dims(ds_org,bk_f,0)*conform_dims(ds_org,f1->ps/100.,1)
  u = linint1_n_Wrap(p_fy,f1->ucomp,False,lev_p,0,0)
  u = u - u0
  symMinMaxPlt(u,21,False,res)
  delete_VarAtts(ores,(/"cnLevelSelectionMode","cnMinLevelValF","cnMaxLevelValF","cnLevelSpacingF"/))
do ip = 0,dimsizes(cidArr)-1
  f1  = addfile("ensemble-wise/zyorg_"+cidArr(ip)+".nc","r")
  p_fy = conform_dims(ds_org,bk_f,0)*conform_dims(ds_org,f1->ps/100.,1)
  u = linint1_n_Wrap(p_fy,f1->ucomp,False,lev_p,0,0)
  u = u - u0
  ds = dimsizes(u)
  u2 = u(:,0:ds(1)/2-1)
  u2 = (u2 +u(:,ds(1)-1:ds(1)/2))/2
  u2@long_name = "~F5~D~F~u (m/s): "+cidArr2(ip)
  print(cidArr(ip)+"	"+max(abs(u2({200:},:))))
;  symMinMaxPlt(u2,21,False,res)
    res@tiYAxisOn            = (mod(ip,dimsizes(cidArr)).eq.0)
    res@tmYLLabelsOn         = (mod(ip,dimsizes(cidArr)).eq.0)       ; panel 7
;    res@tmYLOn               = (mod(ip,dimsizes(cidArr)).eq.0)       ; panel 7
  plot(ip+dimsizes(cidArr)) = gsn_csm_contour(wks, u2, res)
  ploto = gsn_csm_contour(wks, u02, ores)
  overlay(plot(ip+dimsizes(cidArr)),ploto)
;  draw(plot(ip+dimsizes(cidArr)))
;  frame(wks)
  delete([/u2,u,f1/])
end do
  delete_VarAtts(resP,(/"gsnPanelBottom"/))
    resP@gsnPanelTop    = 0.5            ; panel 7
    gsn_panel(wks,plot(dimsizes(cidArr):),(/1,dimsizes(cidArr)/),resP)
  frame(wks)

;  drawNDCGrid(wks)

end
exit

; vim: set fdm=marker foldmarker=/;,;/:

