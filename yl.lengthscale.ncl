; various lengthscale
;TODO ncl -Q yl.lengthscale.ncl
load "./gsn_csm.ncl"
begin
  pp = 300.
  wcStrt = systemfunc("date")

 ; FMS
;  FILES1 = systemfunc("cat "+get_script_prefix_name+".filelist") ;systemfunc (" ls -1 var_a[14]*.nc")
  cidCtrl = "b1ctrl"
  cidArr = (/cidCtrl,"b12unifit1","b15unifp15","b16unifm15","b14unifm3","b2kidston","b13resi"/)
  mycolor = (/"black","blue","dodgerblue","hotpink","red"/)
  FILES1 = "ensemble-wise/cy.cospectrum."+cidArr+"."+pp+".nc"  ; FILES1(0) dont exist
  fdata = addfile(FILES1(1),"r")
  ds = fdata@ds
  lat = fdata->lat
;  u_y = fdata->u_y
  vv_syn = new((/dimsizes(cidArr),ds(2),ds(3)/2+1/),float)
  uu_syn = new((/dimsizes(cidArr),ds(2),ds(3)/2+1/),float)
  vv_syn(0,:,:) = fdata->vvClim_yn
  uu_syn(0,:,:) = fdata->uuClim_yn
  vvb_syn = new((/dimsizes(cidArr),ds(2),ds(3)/2+1/),float)
  vvb_syn(0,:,:) = fdata->vvbClim_yn
  vvs_syn = new((/dimsizes(cidArr),ds(2),ds(3)/2+1/),float)
  vvs_syn(0,:,:) = fdata->vvsClim_yn
  uus_syn = new((/dimsizes(cidArr),ds(2),ds(3)/2+1/),float)
  uus_syn(0,:,:) = fdata->uusClim_yn
  delete(fdata)
  do ip = 1,dimsizes(cidArr)-1
    fdata = addfile(FILES1(ip),"r")
    vv_syn(ip,:,:) = vv_syn(0,:,:)+ fdata->vvEN_yn
    uu_syn(ip,:,:) = uu_syn(0,:,:)+ fdata->uuEN_yn
    vvb_syn(ip,:,:) = vvb_syn(0,:,:)+ fdata->vvbEN_yn
    vvs_syn(ip,:,:) = vvs_syn(0,:,:)+ fdata->vvsEN_yn
    uus_syn(ip,:,:) = uus_syn(0,:,:)+ fdata->uusEN_yn
    delete(fdata)
  end do  ;ip cidArr

  FILES2 = "ensemble-wise/zyorg_"+cidArr+".nc"
  fdata = addfile(FILES2(0),"r")
  ds_org = getfilevardimsizes(fdata,"temp")
  bk   = fdata->bk
  bk_f = fdata->pfull / fdata->phalf(dimsizes(bk)-1)  ; end
  delete(fdata)
  t_spy = new((/dimsizes(cidArr),ds_org(0),ds_org(1)/2/),float)
  ps_sy = new((/dimsizes(cidArr),ds_org(1)/2/),float)
  do ip = 0,dimsizes(cidArr)-1
    fdata = addfile(FILES2(ip),"r")
    t_spy(ip,:,:) = fdata->temp(:,0:ds_org(1)/2-1)
    t_spy(ip,:,:) = (t_spy(ip,:,:) +fdata->temp(:,ds_org(1)-1:ds_org(1)/2))/2
    ps_sy(ip,:) = fdata->ps(0:ds_org(1)/2-1)
    ps_sy(ip,:) = (ps_sy(ip,:) +fdata->ps(ds_org(1)-1:ds_org(1)/2))/2
    delete(fdata)
  end do  ;ip cidArr
  t_spy!0 = "cidArr"
 ; pres_sigma: rightmost must be lat*lon
  p_syp = conform_dims((/dimsizes(cidArr),ds_org(1)/2,ds_org(0)/),bk_f,2)
  p_syp = p_syp*conform(p_syp,ps_sy,(/0,1/))

  FILES3 = "ensemble-wise/lininst_"+cidArr+".nc"
  lLininst_sy = new((/dimsizes(cidArr),ds_org(1)/2/),float)
  do ip = 0,dimsizes(cidArr)-1
    fdata = addfile(FILES3(ip),"r")
    lLininst_sy(ip,:) = fdata->lLininst
    delete(fdata)
  end do  ;ip cidArr
  wallClockElapseTime(wcStrt, "reload nc",0)



/;
  vvClim_yn1 = vvClim_yn(:,1:)
  vvClim_yn1 = vvClim_yn1 *conform(vvClim_yn1,ispan(1,ds(3)/2,1),1)^2 /2./get_pi("f")/6371/conform(vvClim_yn1,cos(lat*get_d2r("f")),0)
  vvEN_yn1 = vvEN_yn(:,1:) *conform(vvClim_yn1,ispan(1,ds(3)/2,1),1)^2 /2./get_pi("f")/6371/conform(vvClim_yn1,cos(lat*get_d2r("f")),0)
;  vvClim_yn1@long_name  = "v'v' at "+pp+" hPa (m~S~2~N~ s~S~-2~N~ per 1 km bin)"
  copy_VarCoords(vvClim_yn1,vvEN_yn1)
  vvEN_yn1@long_name  = "v'v' at "+pp+" hPa (m~S~2~N~ s~S~-2~N~ per 1 km bin)"
  if isvar("vv_ryn") then
  vvCI_yn = sqrt( (dim_variance_n(vv_ryn(irClim,:,:),0)*(dimsizes(irClim)-1) + dim_variance_n(vv_ryn(irEN,:,:),0)*(dimsizes(irEN)-1))/(dimsizes(irClim)+dimsizes(irEN)-2) )*sqrt(1./dimsizes(irClim)+1./dimsizes(irEN)) *cdft_t(0.975,dimsizes(irClim)+dimsizes(irEN)-2)
  vvCI_yn = mask(vvCI_yn,vvCI_yn.ne.0,True)
  vvt_yn = vvEN_yn / vvCI_yn
  copy_VarCoords(vvClim_yn,vvt_yn)
  vvt_yn@long_name  = "t/t~B~0.975~N~ for change in v'v' at "+pp+" hPa"
  end if
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


;  ekeClim_yn = uuClim_yn+vvClim_yn
;  ekeClim_yn = mask(ekeClim_yn,ekeClim_yn.ne.0,True)
;  ekeEN_yn = uuEN_yn+vvEN_yn
;  copy_VarCoords(u_y,ekeClim_yn(:,0))
;  copy_VarCoords(ekeClim_yn,ekeEN_yn)
;  ekeEN_yn@long_name  = "u'u'+v'v' at "+pp+" hPa (m~S~2~N~ s~S~-2~N~)"
;  ekenClim_y = dim_sum_n(ekeClim_yn(:,1:)*conform(ekeClim_yn(:,1:),ispan(1,ds(3)/2,1),1),1)/dim_sum_n(ekeClim_yn(:,1:),1)
;  ekenEN_y = dim_sum_n((ekeClim_yn(:,1:)+ekeEN_yn(:,1:))*conform(ekeClim_yn(:,1:),ispan(1,ds(3)/2,1),1),1)/dim_sum_n(ekeClim_yn(:,1:)+ekeEN_yn(:,1:),1)
;  lBeta_y = 2*get_pi("f")*sqrt( sqrt(dim_sum_n((ekeClim_yn(:,1:)+ekeEN_yn(:,1:)),1))/(2*7.292e-5*cos(lat*get_d2r("f"))/6371e3) )

  eke_syn = (uu_syn+vv_syn)/2
  lBetap_sy = sqrt( sqrt(dim_sum_n(eke_syn(:,:,1:),2))/conform_dims((/dimsizes(cidArr),ds(2)/),2*7.292e-5*cos(lat*get_d2r("f"))/6371e3,1) )  ;2*get_pi("f")*
  ekes_syn = (uus_syn+vvs_syn)/2
  lBeta_sy = sqrt( sqrt(dim_sum_n(ekes_syn(:,:,1:),2))/conform_dims((/dimsizes(cidArr),ds(2)/),2*7.292e-5*cos(lat*get_d2r("f"))/6371e3,1) )  ;2*get_pi("f")*

 ; brunt_vaisala_atm(z,th), rigrad_bruntv_atm(z,th), static_stability(p)
  t_spy@units = "degK"  ; orig. deg_k
  nsq_syp = 9.80^2/t_spy(cidArr|:,lat|:,pfull|:)^2*p_syp/287.04 *static_stability(p_syp,t_spy(cidArr|:,lat|:,pfull|:),2,0) -0.015^2
  pTrop_sy = new((/dimsizes(cidArr),ds_org(1)/2/),float)
  do kp = 16,3,1
    pTrop_sy = where(ismissing(pTrop_sy).and.nsq_syp(:,:,kp).gt.0., (nsq_syp(:,:,kp)*p_syp(:,:,kp+1)-nsq_syp(:,:,kp+1)*p_syp(:,:,kp))/(nsq_syp(:,:,kp)-nsq_syp(:,:,kp+1)), pTrop_sy)
  end do  ;ip cidArr
;  pTrop_sy = rm_single_dims(linint1_n(nsq_syp(:,:,4:16),p_syp(:,:,4:16),False,0.015^2,0,2))
;print(pTrop_sy)
;  pTrop_sy = trop_wmo(p_syp,t_spy(cidArr|:,lat|:,pfull|:),1,False)
 ; cz2ccm, hydro
 ; ${MiMAROOT}/src/shared/constants/constants.f90
 ; https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pjp/spectral_core.pdf, Eq. 33
  wrk = conform(t_spy,pTrop_sy/ps_sy,(/0,2/))
  zTrop_sy = dim_sum_n(log(dim_max_n((/conform(t_spy,bk(1:),1),wrk/),0)/dim_max_n((/conform(t_spy,bk(:ds_org(0)-1),1),wrk/),0)) *287.04/9.80*t_spy,1)
 ; vinth2p: rightmost must be lev*lat*lon
  lRossby_sy = sqrt(9.80*log( rm_single_dims(linint1_n_Wrap(log(p_syp/conform(p_syp,pTrop_sy,(/0,1/))),t_spy(cidArr|:,lat|:,pfull|:),False,0,0,2)) /t_spy(:,ds_org(0)-1,:)/(pTrop_sy/p_syp(:,:,ds_org(0)-1))^(2./7))*zTrop_sy)/conform(lBeta_sy,2*7.292e-5*sin(abs(lat)*get_d2r("f")),1)  ;2*get_pi("f")*

;  vvl_sy = dim_sum_n(vv_syn(:,:,1:)*conform(vv_syn(:,:,1:),2*get_pi("f")*6371e3/ispan(1,ds(3)/2,1) ,2),2)/dim_sum_n(vv_syn(:,:,1:),2)
;  vvl_sy = vvl_sy*conform(vvl_sy,cos(lat*get_d2r("f")),1)
  vvlm2b_sy = sqrt(dim_sum_n(vvb_syn(:,:,1:)*conform(vvb_syn(:,:,1:),(2*get_pi("f")*6371e3/ispan(1,ds(3)/2,1))^2 ,2),2)/dim_sum_n(vvb_syn(:,:,1:),2))
  vvlm2b_sy = vvlm2b_sy*conform(vvlm2b_sy,cos(lat*get_d2r("f")),1)
  vvlm1b_sy = dim_sum_n(vvb_syn(:,:,1:)*conform(vvb_syn(:,:,1:),(2*get_pi("f")*6371e3/ispan(1,ds(3)/2,1))^1 ,2),2)/dim_sum_n(vvb_syn(:,:,1:),2)
  vvlm1b_sy = vvlm1b_sy*conform(vvlm1b_sy,cos(lat*get_d2r("f")),1)
  vvlm2p_sy = sqrt(dim_sum_n(vv_syn(:,:,1:)*conform(vv_syn(:,:,1:),(2*get_pi("f")*6371e3/ispan(1,ds(3)/2,1))^2 ,2),2)/dim_sum_n(vv_syn(:,:,1:),2))
  vvlm2p_sy = vvlm2p_sy*conform(vvlm2p_sy,cos(lat*get_d2r("f")),1)
  vvlm1p_sy = dim_sum_n(vv_syn(:,:,1:)*conform(vv_syn(:,:,1:),(2*get_pi("f")*6371e3/ispan(1,ds(3)/2,1))^1 ,2),2)/dim_sum_n(vv_syn(:,:,1:),2)
  vvlm1p_sy = vvlm1p_sy*conform(vvlm1p_sy,cos(lat*get_d2r("f")),1)

  type = "pdf"
; type@wkOrientation = "portrait" 
  type@wkPaperWidthF  = 190/25.4
  type@wkPaperHeightF = 190/25.4
  wks = gsn_open_wks(type,get_script_prefix_name+"."+pp)  ; TODO
  gsn_define_colormap(wks,"BlWhRe")  ;"BlueWhiteOrangeRed"

; Set some basic resources
/;
  xyres = True
  xyres@tmXBLabelFontHeightF  = 0.04
  xyres@tmXBLabelDeltaF       = -0.5
  xyres@tmYLLabelFontHeightF  = 0.04
  xyres@tmYLLabelDeltaF       = -0.5
  xyres@tiXAxisFontHeightF= 0.04
  xyres@tiYAxisFontHeightF= 0.04
  xyres@tiDeltaF          = 0.
  xyres@vpWidthF        = 0.82
  xyres@vpHeightF       = 0.62
  xyres@gsnDraw         = False
  xyres@gsnFrame        = False
;/

  res = True

  res@cnFillOn          = True
  res@cnLinesOn         = False
;  res@cnFillMode        = "CellFill"
  res@cnFillPalette     = "ncview_default"
;  res@cnLevelSelectionMode = "ManualLevels"
;  res@cnMinLevelValF       = -0.08
;  res@cnMaxLevelValF       = 0.1
;  res@cnLevelSpacingF      = 0.02
;  res@cnGridBoundFillColor = "grey50"
;  res@cnMissingValFillColor= "grey50"
  res@lbOrientation      = "Vertical"
;  res@lbLabelPosition    = "Left"
;  res@pmLabelBarOrthogonalPosF = -1.5
  res@pmLabelBarWidthF   = 0.35
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
;  res@tiXAxisString     = "U/cos~F33~q~F~ and angular phase speed (m/s)"
;  res@tiMainFontHeightF = 0.04
  res@tiXAxisFontHeightF= 0.04
  res@tiYAxisFontHeightF= 0.04
  res@tiDeltaF          = 0.
  res@vpWidthF        = 0.82
  res@vpHeightF       = 0.62
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
  ores = True
;  ores@cnFillOn          = True
;  ores@cnFillMode        = "CellFill"
;  ores@cnLinesOn         = False
  ores@cnLineLabelsOn  = False
;  ores@cnInfoLabelOn   = False
  ores@cnInfoLabelFontHeightF= 0.025
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
;/

  gsres = True
;  gsres@tfPolyDrawOrder  = "Draw"

    resP = True
;   resP@txString    = ""
 ;  resP@txPosYF  = 0.98
    resP@gsnMaximize = True
    resP@gsnPaperOrientation   = "portrait"
    resP@gsnPaperMargin        = 0.
    resP@gsnBoxMargin          = 0.
;    resP@gsnPanelDebug = True
;   resP@gsnPanelLabelBar = True
;   resP@lbLabelFontHeightF    = 0.012
;    resP@gsnPanelFigureStrings = (/"N increase"/)
    resP@gsnPanelFigureStringsFontHeightF  = 0.02
    resP@gsnPanelFigureStringsPerimOn = False
    resP@gsnPanelFigureStringsJust  = "TopLeft"
;   resP@pmLabelBarOrthogonalPosF = -0.02   ; panel
    resP@gsnPanelScalePlotIndex   = 0    ; panel 24
    resP@gsnPanelYWhiteSpacePercent = 8  ; panel 6
;    resP@gsnPanelBottom = 0.04           ; panel 7
;    resP@gsnPanelTop    = 0.935          ; panel 7
;    resP@gsnPanelLeft   = 0.10
;    resP@gsnPanelRight  = 0.8            ; panel 7
 ;  resP@gsnFrame    = False

  plot = new(6, graphic)
  dum  = new(6, graphic)

/;
  symMinMaxPlt(vvEN_yn1,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@gsnXAxisIrregular2Linear = False
  res@gsnYAxisIrregular2Linear = False
  res@trYMaxF         = 20000.
;  res@trYMinF         = 0.
  res@trGridType      = "TriangularMesh" ;"Curvilinear"
  res@sfXArray        = conform(transpose(vvClim_yn1),lat,1)
  res@sfYArray        = 2*get_pi("f")*6371*cos(res@sfXArray*get_d2r("f"))/conform(transpose(vvClim_yn1),ispan(1,ds(3)/2,1),0)
;  res@sfYArray        = 2*get_pi("f")*6371/ispan(1,ds(3)/2,1)
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavelength (km)"
  res@vpWidthF        = 0.82
  res@vpHeightF       = 0.62
  res@tmYLMaxTicks      = 4
;  res@tmYLTickSpacingF  = 5000.
;  res@tmXBMaxTicks      = 4
;  res@tmXBMinorOn       = False
  plot(0) = gsn_csm_contour(wks,transpose(vvEN_yn1),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 9
  ores@cnLineThicknessF     = 4.
  ores@trGridType      = res@trGridType
  ores@sfXArray        = res@sfXArray
  ores@sfYArray        = res@sfYArray
  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(vvClim_yn1),ores)
  overlay(plot(0),ploto)
  gsres@gsLineThicknessF = 8.
  gsres@gsLineColor      = "Black"
  gsres@gsLineDashPattern = 0
;  dum(0)  = gsn_add_polyline(wks,plot(0),lat,vvlClim_y,gsres)
  dum(0)  = gsn_add_polyline(wks,plot(0),lat,2*get_pi("f")*6371*cos(lat*get_d2r("f"))/vvnClim_y,gsres)
  gsres@gsLineColor      = "Green4"
  gsres@gsLineDashPattern = 2
;  dum(1)  = gsn_add_polyline(wks,plot(0),lat,vvlEN_y,gsres)
  dum(1)  = gsn_add_polyline(wks,plot(0),lat,2*get_pi("f")*6371*cos(lat*get_d2r("f"))/vvnEN_y,gsres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
  delete_VarAtts(res,(/"trGridType","sfXArray","sfYArray" ,"tmYLMaxTicks"/))
  delete_VarAtts(ores,(/"trGridType","sfXArray","sfYArray"/))
  res@gsnXAxisIrregular2Linear = True
  res@gsnYAxisIrregular2Linear = True

;  vvEN_yn = mask(vvEN_yn,vvEN_yn.ne.0,True)
  symMinMaxPlt(vvEN_yn,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  res@vpWidthF        = 0.82
  res@vpHeightF       = 0.62
  plot(0) = gsn_csm_contour(wks,transpose(vvEN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 6
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(vvClim_yn),ores)
  overlay(plot(0),ploto)
;  gsres@gsLineThicknessF = 8.
  gsres@gsLineColor      = "Black"
  gsres@gsLineDashPattern = 0
;  dum(0)  = gsn_add_polyline(wks,plot(0),lat,2*get_pi("f")*6371*cos(lat*get_d2r("f"))/vvlClim_y,gsres)
  dum(0)  = gsn_add_polyline(wks,plot(0),lat,vvnClim_y,gsres)
  gsres@gsLineColor      = "Green4"
  gsres@gsLineDashPattern = 2
;  dum(1)  = gsn_add_polyline(wks,plot(0),lat,2*get_pi("f")*6371*cos(lat*get_d2r("f"))/vvlEN_y,gsres)
  dum(1)  = gsn_add_polyline(wks,plot(0),lat,vvnEN_y,gsres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
;/

;  symMinMaxPlt(t_spy(0,:,:),51,False,res)
  res@gsnXAxisIrregular2Linear = False
  res@gsnYAxisIrregular2Linear = False
  res@trYReverse      = True
;  res@trYMaxF         = 1000.
;  res@trYMinF         = 0.
  res@trGridType      = "TriangularMesh" ;"Curvilinear"
  res@sfXArray        = conform(t_spy(0,:,:),lat,1)
  res@sfYArray        = transpose(p_syp(0,:,:))/100.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Pressure (hPa)"
;  res@tmYLMaxTicks      = 4
;  res@tmYLTickSpacingF  = 5000.
  plot(0) = gsn_csm_contour(wks,t_spy(0,:,:),res)
  gsres@gsLineThicknessF = 8.
  gsres@gsLineColor      = "Black"
  gsres@gsLineDashPattern = 0
  dum(0)  = gsn_add_polyline(wks,plot(0),lat,pTrop_sy(0,:)/100.,gsres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
  delete_VarAtts(res,(/"trGridType","sfXArray","sfYArray"/))
;  delete_VarAtts(ores,(/"trGridType","sfXArray","sfYArray"/))
  res@gsnXAxisIrregular2Linear = True
  res@gsnYAxisIrregular2Linear = True
  res@trYReverse      = False

  xyres = get_res_ne(res,(/"lb","cn","pmLabelBar","gsnAddCyclic","gsnXAxis","gsnYAxis"/))
;  delete_VarAtts(xyres,(/"trYMaxF"/))
  xyres@trXMinF         = min(lat)  ; bug 6105 5794 3445 3345 3106
  xyres@trXMaxF         = max(lat)
;  xyres@xyLineColors  = ndtooned(conform_dims((/2,ds(0)/2/),(/"Black","Blue"/),0))
  xyres@xyLineColors  = ndtooned((/mycolor,mycolor/))
  xyres@tiYAxisOffsetXF   = -0.01
;  lRossby_sy = 2*get_pi("f")*sqrt(9.80*log( rm_single_dims(linint1_n_Wrap(log(p_syp/conform(p_syp,pTrop_sy,(/0,1/))),t_spy(cidArr|:,lat|:,pfull|:),False,0,0,2)) /t_spy(:,ds_org(0)-1,:)/(pTrop_sy/p_syp(:,:,ds_org(0)-1))^(2./7))*zTrop_sy)/conform(lBeta_sy,2*7.292e-5*sin(abs(lat)*get_d2r("f")),1)

  xyres@pmLegendDisplayMode    = "Always"            ; turn on legend
 ;xyres@pmLegendZone           = 1
;;xyres@pmLegendSide           = "Right"             ; Change location of 
 ;xyres@pmLegendParallelPosF   = 0.98                ; move
;;xyres@pmLegendOrthogonalPosF = -0.2                ; move
 ;xyres@lgJustification        = "BottomRight"
 ;xyres@pmLegendWidthF         = 0.12                ; Change width and
 ;xyres@pmLegendHeightF        = 0.1                 ; height of legend.
;;xyres@lgPerimOn              = False               ; turn off box around
;;xyres@lgLabelFontHeightF     = .02                 ; label font height
  xyres@xyExplicitLegendLabels = cidArr         ; create explicit labels
  xyres@tiYAxisString     = "T@bottom (K)"
  plot(0) = gsn_csm_xy(wks,lat,t_spy(:,ds_org(0)-1,:),xyres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
  xyres@pmLegendDisplayMode    = "NoCreate"

  xyres@trYReverse      = True
  xyres@tiYAxisString     = "p Trop"
  plot(0) = gsn_csm_xy(wks,lat,pTrop_sy,xyres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  xyres@trYReverse      = False
  xyres@tiYAxisString     = "theta Trop"
  thetaTrop = rm_single_dims(linint1_n_Wrap(log(p_syp/conform(p_syp,pTrop_sy,(/0,1/))),t_spy(cidArr|:,lat|:,pfull|:),False,0,0,2))/(pTrop_sy/100000.)^(2./7)
  plot(0) = gsn_csm_xy(wks,lat,thetaTrop,xyres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  xyres@tiYAxisString     = "relative theta Trop (K)"
  plot(0) = gsn_csm_xy(wks,lat,thetaTrop-thetaTrop(conform(cidArr,0,-1),:),xyres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  xyres@tiYAxisString     = "h Trop"
  plot(0) = gsn_csm_xy(wks,lat,zTrop_sy,xyres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  xyres@tiYAxisString     = "L~B~D~N~ * f"
  plot(0) = gsn_csm_xy(wks,lat,lRossby_sy*conform(lBeta_sy,2*7.292e-5*sin(abs(lat)*get_d2r("f")),1),xyres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  xyres@tiYAxisString     = "L~B~e"
  plot(0) = gsn_csm_xy(wks,lat,vvlm2b_sy,xyres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
  xyres@tiYAxisString     = "L~BF5~b"
  plot(0) = gsn_csm_xy(wks,lat,lBeta_sy,xyres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
  xyres@tiYAxisString     = "L~B~grow"
  plot(0) = gsn_csm_xy(wks,lat,lLininst_sy,xyres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
  do ip = 0,dimsizes(cidArr)-1
;    do yy = 2,ds_org(1)/2-1
    do yy = 5,ds_org(1)/2-1
      if abs(lLininst_sy(ip,yy)/lLininst_sy(ip,yy-1)-1).gt.0.04 then
        lLininst_sy(ip,yy:) = lLininst_sy@_FillValue
        print(cidArr(ip)+":	"+lat(yy))
        break
      end if
    end do  ;
  end do  ;ip cidArr
  plot(0) = gsn_csm_xy(wks,lat,lLininst_sy,xyres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  xyres@tiYAxisString     = "L~B~e~N~/L~BF5~b"
;  plot(0) = gsn_csm_xy(wks,lat,2*get_pi("f")*6371e3*cos(lat*get_d2r("f"))/ekenEN_y/lBeta_y,xyres)
  plot(0) = gsn_csm_xy(wks,lat,vvlm2b_sy/lBeta_sy,xyres)
;  plot(0) = gsn_csm_xy(wks,lat,array_append_record(vvlm2b_sy/lBeta_sy,vvlm2b_sy/lRossby_sy,0),xyres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  xyres@tiYAxisString     = "L~B~e~N~/L~B~D"
  plot(0) = gsn_csm_xy(wks,lat,vvlm2b_sy/lRossby_sy,xyres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  xyres@tiYAxisString     = "L~B~e~N~/L~B~grow"
  xyres@trYMaxF         = 3.
  plot(0) = gsn_csm_xy(wks,lat,vvlm2b_sy/lLininst_sy,xyres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  xyres@tiYAxisString     = "L~B~forced~N~/L~B~ctrl";"new L/ctrl L"
  xyres@tiYAxisFontHeightF= 0.04
;  xyres@xyLabelMode       = "Custom"
;  xyres@xyExplicitLabels  = (/"m2b","m1b","m2p","m1p","L~BF5~b","300L~BF5~b"/)
;  xyres@xyLineLabelFontHeightF = 0.02
  xyres@xyLineColors := (/"black","hotpink","black","red","blue","blueviolet"/)
  xyres@pmLegendDisplayMode    = "Always"            ; turn on legend
  xyres@xyExplicitLegendLabels:= (/"L~BF5~b","300L~BF5~b","bt-2","300-2","bt-1","300-1"/)    ; create explicit labels
  xyres@trYMaxF         = 1.1
  xyres@trYMinF         = 0.9
  xyres@gsnYRefLine       = 1.
  ip = 1
  xyres@tiYAxisFontColor  = mycolor(ip)
  plot(ip-1) = gsn_csm_xy(wks,lat,(/lBeta_sy(ip,:)/lBeta_sy(0,:),lBetap_sy(ip,:)/lBetap_sy(0,:) ,vvlm2b_sy(ip,:)/vvlm2b_sy(0,:),vvlm2p_sy(ip,:)/vvlm2p_sy(0,:),vvlm1b_sy(ip,:)/vvlm1b_sy(0,:),vvlm1p_sy(ip,:)/vvlm1p_sy(0,:)/),xyres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
  xyres@xyLineColors := ndtooned((/mycolor,mycolor/))
  xyres@pmLegendDisplayMode    = "NoCreate"

;  delete_VarAtts(xyres,(/"vpWidthF","vpHeightF"/))
  xyres@vpWidthF        = 1.
  xyres@vpHeightF       = 0.30
  xyres@tiYAxisString     = "L~B~forced~N~/L~B~ctrl";"new L/ctrl L"
;  xyres@tiYAxisFontHeightF= 0.04
  xyres@xyLabelMode       = "Custom"
  xyres@xyExplicitLabels := (/"L~B~e","L~BF5~b","L~B~D","L~B~grow"/);"L~BF5~b~EF~L~B~D~S1~0.5"
  xyres@xyLineLabelFontHeightF = 0.02
  xyres@trYMaxF         = 1.1
  xyres@trYMinF         = 0.9
  xyres@gsnYRefLine       = 1.
  do ip = 1,dimsizes(mycolor)-1
;  xyres@tiYAxisString     = cidArr(ip)
  xyres@tiYAxisFontColor  = mycolor(ip)
  xyres@tmXBLabelsOn    = (ip.ge.dimsizes(cidArr)-1)
  xyres@tmXBOn          = (ip.ge.dimsizes(cidArr)-1)
  plot(ip-1) = gsn_csm_xy(wks,lat,(/vvlm2b_sy(ip,:)/vvlm2b_sy(0,:),lBeta_sy(ip,:)/lBeta_sy(0,:),lRossby_sy(ip,:)/lRossby_sy(0,:),lLininst_sy(ip,:)/lLininst_sy(0,:)/),xyres) ;,lBeta_sy(ip,:)/lBeta_sy(0,:)*(lRossby_sy(ip,:)/lRossby_sy(0,:))^0.5
  end do  ;ip cidArr
;    resP@txString    = "new L/ ctrl L"
;    resP@gsnPanelFigureStrings = cidArr(1:)
;    resP@gsnPanelFigureStringsFontHeightF  = 0.03
;    resP@gsnPanelFigureStringsJust  = "TopLeft"
    resP@gsnPanelBottom = 0.08
    gsn_panel(wks,plot(:3),(/4,1/),resP)
  delete_VarAtts(resP,(/"gsnPanelBottom"/))

  type@wkOrientation = "portrait" 
  type@wkPaperWidthF  = 39*12/72.27  ; 190/25.4
  type@wkPaperHeightF = 5.8
  wks2 = gsn_open_wks(type,get_script_prefix_name+".pub"+pp)  ; TODO
;  xyres@vpWidthF        = 0.82
;  xyres@vpHeightF       = 0.62
  delete_VarAtts(xyres,(/"vpWidthF","vpHeightF","tmXBLabelsOn","tmXBOn" ,"tiYAxisFontColor"/))
;  ip = 1
  s_i = (/1,5,6/)
  p_i = (/0,2,3/)
;  xyres@tiYAxisFontColor  = "black";mycolor(ip)
  xyres@tiYAxisFontHeightF= 0.04
  do ip = 0,dimsizes(p_i)-1
;    xyres@tiYAxisOn            = (mod(p_i(ip),dimsizes(cidArr)-1).eq.0)
;    xyres@tmYLLabelsOn         = (mod(p_i(ip),dimsizes(cidArr)-1).eq.0)       ; panel 7
  plot(p_i(ip)) = gsn_csm_xy(wks2,lat,(/vvlm2b_sy(s_i(ip),:)/vvlm2b_sy(0,:),lBeta_sy(s_i(ip),:)/lBeta_sy(0,:),lRossby_sy(s_i(ip),:)/lRossby_sy(0,:),lLininst_sy(s_i(ip),:)/lLininst_sy(0,:)/),xyres)
  end do  ;ip cidArr

  delete_VarAtts(xyres,(/"trXMinF","trXMaxF"/))
  xyres@tiXAxisOffsetYF   = -0.01
  xyres@tiXAxisString     = "Forcing applied (~F18~U~F~LRF)"
  xyres@tiYAxisString     = "L~B~forced~N~/L~B~ctrl"
;  xyres@tiYAxisFontHeightF= 0.04
  xyres@xyMarkLineMode    = "MarkLines"
  xyres@xyMarkerColors    = xyres@xyLineColors
  xyres@xyMarkers         = (/12,4,5,7/)
  xyres@xyMarkerSizeF     = 0.02
  xyres@xyLabelMode       = "NoLabels";"Custom"
;  xyres@xyExplicitLabels := (/"L~B~e","L~BF5~b","L~B~D","L~B~grow"/);"L~BF5~b~EF~L~B~D~S1~0.5"
;  xyres@xyLineLabelFontHeightF = 0.02
  xyres@pmLegendDisplayMode    = "Always"            ; turn on legend
  xyres@pmLegendZone           = 1
; xyres@pmLegendSide           = "Right"             ; Change location of 
 ;xyres@pmLegendParallelPosF   = 0.98                ; move
;;xyres@pmLegendOrthogonalPosF = -0.2                ; move
  xyres@lgJustification        = "BottomLeft"
  xyres@pmLegendWidthF         = 0.12                ; Change width and
 ;xyres@pmLegendHeightF        = 0.1                 ; height of legend.
;;xyres@lgPerimOn              = False               ; turn off box around
  xyres@lgLabelFontHeightF     = .03                 ; label font height
  xyres@xyExplicitLegendLabels:= xyres@xyExplicitLabels
  xyres@trYMaxF         = 1.1
  xyres@trYMinF         = 0.9
  xyres@gsnYRefLine       = 1.
  xyres@tmXBMode          = "Manual"
  xyres@tmXBTickSpacingF  = 0.5
;  xyres@tmXBValues        = fspan(1.,-1.,5)
  ip := (/1,2,0,3,4/)
  yy = minind(abs(lat+40)) ;26
;  plot(1) = gsn_csm_xy(wks2,fspan(1.,-1.,5),(/dim_avg(vvlm2b_sy(ip,yy-2:yy+2)/vvlm2b_sy(ip*0,yy-2:yy+2)),dim_avg(lBeta_sy(ip,yy-2:yy+2)/lBeta_sy(ip*0,yy-2:yy+2)),dim_avg(lRossby_sy(ip,yy-2:yy+2)/lRossby_sy(ip*0,yy-2:yy+2)),dim_avg(lLininst_sy(ip,yy-2:yy+2)/lLininst_sy(ip*0,yy-2:yy+2))/),xyres)
  plot(1) = gsn_csm_xy(wks2,fspan(1.,-1.,5),(/vvlm2b_sy(ip,yy)/vvlm2b_sy(0,yy),lBeta_sy(ip,yy)/lBeta_sy(0,yy),lRossby_sy(ip,yy)/lRossby_sy(0,yy),lLininst_sy(ip,yy)/lLininst_sy(0,yy)/),xyres)
  _ls = (/vvlm2b_sy(:,yy), lBeta_sy(:,yy), lRossby_sy(:,yy), lLininst_sy(:,yy)/)
;  _ls(:,1:) = _ls(:,1:)/_ls(:,0)
  alist = [/xyres@xyExplicitLabels, _ls(:,s_i(0))/_ls(:,0), _ls(:,s_i(1))/_ls(:,0), _ls(:,s_i(2))/_ls(:,0)/]
  print_table(alist, "%8s%.3f%.3f%.3f")
;  yy = minind(abs(lat+60)) ;
;  plot(2) = gsn_csm_xy(wks2,fspan(1.,-1.,5),(/vvlm2b_sy(ip,yy)/vvlm2b_sy(0,yy),lBeta_sy(ip,yy)/lBeta_sy(0,yy),lRossby_sy(ip,yy)/lRossby_sy(0,yy),lLininst_sy(ip,yy)/lLininst_sy(0,yy)/),xyres)
;    resP@gsnPanelFigureStrings = (/"(a) LRF","(b) "+get_lat_labels(round(lat(yy),3))/)
    cidArr2 = (/"LRF",get_lat_labels(round(lat(yy),3)),"K11","K11-LRF"/)
    resP@gsnPanelFigureStrings = "("+(/"a","b","c","d"/)+") "+cidArr2
;    resP@gsnPanelScalePlotIndex   = 1    ; panel 24
;    resP@gsnPanelLeft   = 0.01
    resP@gsnPanelRight  = 0.99           ; panel 7
;    gsn_panel(wks2,plot(:1),(/1,2/),resP)
  newplots = gsn_panel_return(wks2,plot(:3),(/2,2/),resP)
;  drawNDCGrid(wks)
;  frame(wks)
  bb_dbg = NhlGetBB(newplots)
      print("min/max x position = " + min(bb_dbg(:,2)) + "/" + max(bb_dbg(:,3)))
;      print("min/max y position = " + min(bb_dbg(:,1)) + "/" + max(bb_dbg(:,0)))
print("wkPaperHeightF="+type@wkPaperWidthF*(max(bb_dbg(:,0))-min(bb_dbg(:,1))))
  delete(ip)

  xyres@trYReverse      = False
;  xyres@trYMaxF         = 20000.
;  xyres@trYMinF         = 0.
  xyres@tiYAxisString     = "~F5~D~F~ Wavenumber"
  xyres@tiYAxisString     = "~F5~D~F~ Wavelength (km)"
;  xyres@tmYLMaxTicks      = 4
;  plot(0) = gsn_csm_xy(wks,lat,vvn_ry-conform(vvn_ry,vvnClim_y,1),xyres)
;  plot(1) = gsn_csm_xy(wks,lat,vvl_ry-conform(vvl_ry,vvlClim_y,1),xyres)
;  xyres@tiYAxisString     = "t/t~B~0.975~N~ for wavenumber"
;  xyres@tiYAxisString     = "t/t~B~0.975~N~ for wavelength"
;  xyres@gsnYRefLine       = (/-1,1/)
;  plot(2) = gsn_csm_xy(wks,lat,vvnt_y,xyres)
;  plot(3) = gsn_csm_xy(wks,lat,vvlt_y,xyres)
;  delete_VarAtts(xyres,(/"gsnYRefLine"/))
;    resP@gsnPanelDebug = True
;    gsn_panel(wks,plot(:3),(/2,2/),resP)

/;
  if isvar("vv_ryn") then
  xyres = get_res_ne(res,(/"lb","cn","pmLabelBar","gsnAddCyclic","gsnXAxis","gsnYAxis"/))
  delete_VarAtts(xyres,(/"trYMaxF"/))
  xyres@trXMinF         = min(lat)  ; bug 6105 5794 3445 3345 3106
  xyres@trXMaxF         = max(lat)
  xyres@xyLineColors  = ndtooned(conform_dims((/2,ds(0)/2/),(/"Black","Blue"/),0))
  xyres@tiYAxisOffsetXF   = -0.01
  xyres@tiYAxisString     = "~F5~D~F~ Wavenumber"
  plot(0) = gsn_csm_xy(wks,lat,vvn_ry-conform(vvn_ry,vvnClim_y,1),xyres)
  xyres@tiYAxisString     = "t/t~B~0.975~N~ for wavenumber"
  xyres@gsnYRefLine       = (/-1,1/)
  plot(2) = gsn_csm_xy(wks,lat,vvnt_y,xyres)
  delete_VarAtts(xyres,(/"gsnYRefLine"/))

  xyres@trYReverse      = False
;  xyres@trYMaxF         = 20000.
;  xyres@trYMinF         = 0.
  xyres@tiYAxisString     = "~F5~D~F~ Wavelength (km)"
;  xyres@tmYLMaxTicks      = 4
  plot(1) = gsn_csm_xy(wks,lat,vvl_ry-conform(vvl_ry,vvlClim_y,1),xyres)
  xyres@tiYAxisString     = "t/t~B~0.975~N~ for wavelength"
  xyres@gsnYRefLine       = (/-1,1/)
  plot(3) = gsn_csm_xy(wks,lat,vvlt_y,xyres)
;    resP@gsnPanelDebug = True
    gsn_panel(wks,plot(:3),(/2,2/),resP)
  end if
;/

end
exit

; vim: set fdm=marker foldmarker=/;,;/:

