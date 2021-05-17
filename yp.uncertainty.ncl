; plot uncertainty of difference in climatology
;TODO ncl -Q yp.uncertainty.ncl
load "./gsn_csm.ncl"

begin
  wcStrt = systemfunc("date")
  cidCtrl = "b1ctrl"
  caseid = "b12unifit1";"d2unifit1";"b2kidston";"b10unifZeroU";"b7nplus3zeroU";"b4nplus3";"b5unif";getenv("CASENAME")
  if (ismissing(caseid)) then
    print(caseid)
    status_exit(1)
  end if
;  nfiles = toint(systemfunc("ls -1dv manyfiles/zyvar_"+caseid+"_???.nc |wc -l"))
;  FILES1 = "manyfiles/zyvar_"+ndtooned((/"c1T42L20_"+ispan(1,nfiles,1),caseid+"_"+ispan(1,nfiles,1)/))+".nc"
  nfiles = toint(systemfunc("ls -1dv fms-output/"+caseid+"ens*/ |wc -l"))
  FILESo = "fms-output/"+ndtooned((/cidCtrl+"ens"+ispan(1,nfiles,1),caseid+"ens"+ispan(1,nfiles,1)/))+"/zyorg_"+ndtooned((/cidCtrl+"ens"+ispan(1,nfiles,1),caseid+"ens"+ispan(1,nfiles,1)/))+".nc"
  lev_p = 25.*ispan(1,39,1)
  fv = addfiles(FILESo,"r")
  ListSetType(fv,"join") ;cat")
  ds_org = array_append_record(2*nfiles, getfilevardimsizes(fv[0],"ucomp"),0)
  bk   = fv[0]->bk
  bk_f = fv[0]->pfull / fv[0]->phalf(dimsizes(bk)-1)  ; end
  po_rzy = conform_dims(ds_org,bk_f,1)*conform_dims(ds_org,fv[:]->ps/100.,(/0,2/))
  to_rzy = linint1_n_Wrap(po_rzy,fv[:]->temp,False,lev_p,0,1)
  uo_rzy = linint1_n_Wrap(po_rzy,fv[:]->ucomp,False,lev_p,0,1)
  ds = dimsizes(to_rzy)  ;yr,z,y
  delete(fv)
;  wallClockElapseTime(wcStrt, "load",0)

  irEN = ispan(ds(0)/2,ds(0)-1,1) ;1
  irClim = ispan(0,ds(0)/2-1,1) ;0

  irEN := ndtooned((/irEN, irEN+ds(0)/))
  irClim := ndtooned((/irClim, irClim+ds(0)/))
  ds = (/ds(0)*2,ds(1),ds(2)/2/)
  t8 = array_append_record(to_rzy(:,:,0:ds(2)-1),to_rzy(:,:,2*ds(2)-1:ds(2)),0)
  to_rzy := t8
  u8 = array_append_record(uo_rzy(:,:,0:ds(2)-1),uo_rzy(:,:,2*ds(2)-1:ds(2)),0)
  uo_rzy := u8
  delete([/t8,u8/])

  tCIo_zy = sqrt( (dim_variance_n(to_rzy(irClim,:,:),0)*(dimsizes(irClim)-1) + dim_variance_n(to_rzy(irEN,:,:),0)*(dimsizes(irEN)-1))/(dimsizes(irClim)+dimsizes(irEN)-2) )*sqrt(1./dimsizes(irClim)+1./dimsizes(irEN)) *cdft_t(0.975,dimsizes(irClim)+dimsizes(irEN)-2)
  tCIo_zy = where(tCIo_zy.ne.0,tCIo_zy,tCIo_zy@_FillValue)
  copy_VarCoords(to_rzy(0,:,:),tCIo_zy)
  tCIo_zy@long_name  = "t~B~0.975~NF18~U~F33~s~F~ for ~F5~D~F~T (K)"
;  vvt_yn = vvEN_yn / tCIo_zy
;  copy_VarCoords(vvClim_yn,vvt_yn)
;  vvt_yn@long_name  = "t/t~B~0.975~N~ for change in EKE at 300 hPa"
  uCIo_zy = sqrt( (dim_variance_n(uo_rzy(irClim,:,:),0)*(dimsizes(irClim)-1) + dim_variance_n(uo_rzy(irEN,:,:),0)*(dimsizes(irEN)-1))/(dimsizes(irClim)+dimsizes(irEN)-2) )*sqrt(1./dimsizes(irClim)+1./dimsizes(irEN)) *cdft_t(0.975,dimsizes(irClim)+dimsizes(irEN)-2)
  uCIo_zy = where(uCIo_zy.ne.0,uCIo_zy,uCIo_zy@_FillValue)
  copy_VarCoords(uo_rzy(0,:,:),uCIo_zy)
  uCIo_zy@long_name  = "t~B~0.975~NF18~U~F33~s~F~ for ~F5~D~F~u (m s~S~-1~N~)"

  FILES1 = "fms-output/"+ndtooned((/cidCtrl+"ens"+ispan(1,nfiles,1),caseid+"ens"+ispan(1,nfiles,1)/))+"/zyvar_"+ndtooned((/cidCtrl+"ens"+ispan(1,nfiles,1),caseid+"ens"+ispan(1,nfiles,1)/))+".nc"
;print(FILES1)
  fv = addfiles(FILES1,"r")
  ListSetType(fv,"join") ;cat")
  t_rzy = fv[:]->t(:,:,:)
  u_rzy = fv[:]->u(:,:,:)

  ds = dimsizes(t_rzy)  ;yr,z,y
  delete(fv)
;  wallClockElapseTime(wcStrt, "load",0)

;  irEN = ispan(ds(0)/2,ds(0)-1,1) ;1
;  irClim = ispan(0,ds(0)/2-1,1) ;0

;  irEN := ndtooned((/irEN, irEN+ds(0)/))
;  irClim := ndtooned((/irClim, irClim+ds(0)/))
  ds = (/ds(0)*2,ds(1),ds(2)/2/)
  t8 = array_append_record(t_rzy(:,:,0:ds(2)-1),t_rzy(:,:,2*ds(2)-1:ds(2)),0)
  t_rzy := t8
  u8 = array_append_record(u_rzy(:,:,0:ds(2)-1),u_rzy(:,:,2*ds(2)-1:ds(2)),0)
  u_rzy := u8
  delete([/t8,u8/])

;  f1  = addfile("diff_zy_"+caseid+".nc","r")
;  t = rm_single_dims(f1->t)
;  ds = dimsizes(t)
;  t2 = t(:,0:ds(1)/2-1)
;  t2 = (t2 +t(:,ds(1)-1:ds(1)/2))/2
;  tc := t2
;  tc@long_name = "~F5~D~F~T: "+caseid+" (K)"
;  delete([/t2,t,f1/])

  tCI_zy = sqrt( (dim_variance_n(t_rzy(irClim,:,:),0)*(dimsizes(irClim)-1) + dim_variance_n(t_rzy(irEN,:,:),0)*(dimsizes(irEN)-1))/(dimsizes(irClim)+dimsizes(irEN)-2) )*sqrt(1./dimsizes(irClim)+1./dimsizes(irEN)) *cdft_t(0.975,dimsizes(irClim)+dimsizes(irEN)-2)
  tCI_zy = where(tCI_zy.ne.0,tCI_zy,tCI_zy@_FillValue)
  copy_VarCoords(t_rzy(0,:,:),tCI_zy)
  tCI_zy@long_name  = "t~B~0.975~NF18~U~F33~s~F~ for ~F5~D~F~T (K)"
;  vvt_yn = vvEN_yn / tCI_zy
;  copy_VarCoords(vvClim_yn,vvt_yn)
;  vvt_yn@long_name  = "t/t~B~0.975~N~ for change in EKE at 300 hPa"
  uCI_zy = sqrt( (dim_variance_n(u_rzy(irClim,:,:),0)*(dimsizes(irClim)-1) + dim_variance_n(u_rzy(irEN,:,:),0)*(dimsizes(irEN)-1))/(dimsizes(irClim)+dimsizes(irEN)-2) )*sqrt(1./dimsizes(irClim)+1./dimsizes(irEN)) *cdft_t(0.975,dimsizes(irClim)+dimsizes(irEN)-2)
  uCI_zy = where(uCI_zy.ne.0,uCI_zy,uCI_zy@_FillValue)
  copy_VarCoords(u_rzy(0,:,:),uCI_zy)
  uCI_zy@long_name  = "t~B~0.975~NF18~U~F33~s~F~ for ~F5~D~F~u (m s~S~-1~N~)"

; v&$v!1$@long_name = v&$v!1$@long_name + " (hPa)"
;  copy_VarMeta(v(0,:,:,0),uv)
;  uv@long_name = "Eddy zonal momentum fluxes"
;  uv@units     = "m~S~2~N~/s~S~2"
; printVarSummary(uv)
; print(uv&$uv!0$)

  type = "pdf"
; type@wkWidth        = 1200
; type@wkHeight       = 1000
; type@wkPaperWidthF  = 11.2
; type@wkPaperHeightF = 6.3
  wks = gsn_open_wks(type, get_script_prefix_name)
; gsn_define_colormap(wks,"hotres")
  res = True
  res@tiYAxisString = "Pressure (hPa)"
  res@tiYAxisFontHeightF = 0.05
  res@tiDeltaF           = 0.
  res@tmXBLabelFontHeightF  = 0.05
  res@tmXBLabelDeltaF       = -0.5
  res@tmYLLabelFontHeightF  = 0.05
  res@tmYLLabelDeltaF       = -0.5

  res@cnFillOn        = True
; res@cnFillMode      = "RasterFill" ;"AreaFill" ;
  res@cnFillPalette   = "WhBlGrYeRe"
  res@cnLinesOn       = False
  res@cnLineLabelsOn  = False
  res@cnInfoLabelOn   = False
 ;res@cnLineLabelFontHeightF = 0.01
 ;res@cnLevelSelectionMode = "EqualSpacedLevels""ExplicitLevels"  ;
 ;res@cnMaxLevelCount      = 99
; res@cnLevelSelectionMode = "ManualLevels"  ;
; res@cnMinLevelValF       = -3*50
; res@cnMaxLevelValF       = 47*50
; res@cnLevelSpacingF      = 0.2*50
; res@cnLevels             = (/ 0.5 /)
; res@cnFillColors         = (/"White","Green"/)
; res@lbLabelBarOn = False
 ;res@lbAutoManage = False
 ;res@lbPerimOn    = False
 ;res@lbPerimDashSegLenF = 0.
  res@lbOrientation      = "Vertical"
;  res@pmLabelBarOrthogonalPosF = 0.1
  res@pmLabelBarWidthF   = 0.27
; res@lbLabelStride      = 1
 ;res@lbLabelAutoStride  = False
 ;res@lbBoxSizing        = "ExplicitSizing"
 ;res@lbBoxFractions     = array_append_record(fspan(0,0.30,6),fspan(0.36,1,19),0)
  res@lbBoxSeparatorLinesOn = False
  res@lbLabelFontHeightF = 0.04

;  res@vpXF            = 0.26
  res@vpWidthF        = 0.82
;  res@vpYF            = 0.51
  res@vpHeightF       = 0.62
  res@gsnMaximize     = True
  res@gsnPaperOrientation   = "portrait"
  res@gsnPaperMargin        = 0.
  res@gsnBoxMargin          = 0.
;  res@gsnFrame        = False        ; Don't draw or advance
;  res@gsnDraw         = False        ; frame yet.
; res@pmLegendDisplayMode    = "Always"

; res@tiMainString  = p&z(k)+"m"
  res@gsnXAxisIrregular2Linear = True
  res@gsnYAxisIrregular2Linear = True
  res@trYReverse    = True
  res@trYMinF       = 0.
  res@trYMaxF       = 1000.
;  symMinMaxPlt(tforcing1,41,False,res)
  mnmxint = nice_mnmxintvl(0.,max((/max(tCIo_zy),max(tCI_zy)/)),16,False)
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = mnmxint(2)  ; skip 0
  res@cnMaxLevelValF       = mnmxint(1)
  res@cnLevelSpacingF      = mnmxint(2)
  print(caseid+"	org	"+max(abs(tCIo_zy)))
  print(caseid+"	var	"+max(abs(tCI_zy)))
  plot = gsn_csm_contour(wks, tCIo_zy, res)
;  plot = gsn_csm_contour(wks, tCI_zy, res)
;  plot = gsn_csm_pres_hgt(wks, tCI_zy, res)
  mnmxint = nice_mnmxintvl(0.,max((/max(uCIo_zy),max(uCI_zy)/)),16,False)
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = mnmxint(2)
  res@cnMaxLevelValF       = mnmxint(1)
  res@cnLevelSpacingF      = mnmxint(2)
  print(caseid+"	org	"+max(abs(uCIo_zy)))
  print(caseid+"	var	"+max(abs(uCI_zy)))
  plot = gsn_csm_contour(wks, uCIo_zy, res)
;  plot = gsn_csm_contour(wks, uCI_zy, res)

end
exit

; vim: set fdm=marker foldmarker=/;,;/:

