; latitude-phase speed spectrum
; sbatch --account=kuang_lab -p huce_intel -J cospectrum --array=1-1 -n 1 -t 1440 --mem=20000 -o "slurm.%a" --wrap='ncl -Q cy.cospectrum.ncl  >& slurm-${SLURM_JOB_NAME}-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}'
;TODO for jj in `seq 2 13`; do  export CASENAME=`basename ~/fms-backup/b${jj}[a-z]*`; echo $CASENAME; ncl -Q cy.cospectrum.ncl pp=800.; done
load "./gsn_csm.ncl"
begin
  wcStrt = systemfunc("date")

 ; ERA-interim
/;
 ; ~/sw/ecmwfapi/ecmwf.djf.sh
 ; ~/sw/ecmwfapi/interim-uv250.py
;  fv = addfile("/n/regal/kuang_lab/pchan/blocking/ERA-interim/uv250_00Z_ERA-interim_1998djf.nc","r")
  yrArr = ispan(1980,2016,1)
;  FILES1 = systemfunc (" ls -1 /n/regal/kuang_lab/pchan/blocking/ERA-interim/uv250_00Z_ERA-interim_????djf.nc")
  FILES1 = "/n/regal/kuang_lab/pchan/blocking/ERA-interim/uv250_00Z_ERA-interim_"+yrArr+"djf.nc"
  fv = addfiles(FILES1,"r")
  ListSetType(fv,"join") ;cat")
  u = short2flt(fv[:]->u)
  v = short2flt(fv[:]->v)
  lat = u&latitude
;/

 ; FMS
;  FILES1 = systemfunc("cat "+get_script_prefix_name+".filelist") ;systemfunc (" ls -1 var_a[14]*.nc")
  cidCtrl = "b1ctrl"
  caseid = getenv("CASENAME")
  if (ismissing(caseid)) then
    print(caseid)
    status_exit(1)
  end if
  if (str_get_cols(caseid,0,0).eq."z") then
    cidCtrl = "zd100"
  end if
;  pp = 300.
  fn_collect = "ensemble-wise/"+get_script_prefix_name+"."+caseid+"."+pp+".nc"
  if (.not. fileexists(fn_collect)) then
;  nfiles = toint(systemfunc("ls -1dv fms-output/"+caseid+"ens*/"+get_script_prefix_name+"."+caseid+"ens*.data.nc |wc -l"))
  nfiles = toint(systemfunc("ls -1dv fms-output/"+caseid+"ens*/nc.2cospectrum.nc |wc -l"))
;  nfiles = 2 ;TODO
;  FILES1 = "fms-output/"+ndtooned((/cidCtrl+"ens"+ispan(1,nfiles,1),caseid+"ens"+ispan(1,nfiles,1)/))+"/"+get_script_prefix_name+"."+ndtooned((/cidCtrl+"ens"+ispan(1,nfiles,1),caseid+"ens"+ispan(1,nfiles,1)/))+".data.nc"
  FILES1 = "fms-output/"+ndtooned((/cidCtrl+"ens"+ispan(1,nfiles,1),caseid+"ens"+ispan(1,nfiles,1)/))+"/nc.2cospectrum.nc"
  fdata = addfiles(FILES1,"r")
  ListSetType(fdata,"join") ;cat")
;  caseid = fdata@caseid
  dt = fdata[0]@dt
;  Ly = fdata[0]@Ly
  Lx = fdata[0]@Lx
  ds = fdata[0]@ds
  ds(0) = ListCount(fdata)
  irEN = ispan(ds(0)/2,ds(0)-1,1) ;1
  irClim = ispan(0,ds(0)/2-1,1) ;0
  lat = fdata[0]->lat
  u_y = dim_avg_n_Wrap(fdata[irClim]->u_py(:,{pp},:),0)
  uEN_y = dim_avg_n_Wrap(fdata[irEN]->u_py(:,{pp},:),0)
  uv_rywn = fdata[:]->uv_pywn(:,{pp},:,:,:)
  uvClim_ywn = dim_avg_n_Wrap(fdata[irClim]->uv_pywn(:,{pp},:,:,:),0)
  uvEN_ywn = dim_avg_n_Wrap(fdata[irEN]->uv_pywn(:,{pp},:,:,:),0) - uvClim_ywn
  uvEN_ywn@long_name  = uvClim_ywn@long_name
  delete_VarAtts(uvClim_ywn,"long_name")
  vvClim_ywn = dim_avg_n_Wrap(fdata[irClim]->vv_pywn(:,{pp},:,:,:),0)
  vvEN_ywn = dim_avg_n_Wrap(fdata[irEN]->vv_pywn(:,{pp},:,:,:),0) - vvClim_ywn
  vvEN_ywn@long_name  = vvClim_ywn@long_name
  delete_VarAtts(vvClim_ywn,"long_name")
  uuClim_ywn = dim_avg_n_Wrap(fdata[irClim]->uu_pywn(:,{pp},:,:,:),0)
  uuEN_ywn = dim_avg_n_Wrap(fdata[irEN]->uu_pywn(:,{pp},:,:,:),0) - uuClim_ywn
  uuEN_ywn@long_name  = uuClim_ywn@long_name
  delete_VarAtts(uuClim_ywn,"long_name")
  vtClim_ywn = dim_avg_n_Wrap(fdata[irClim]->vt_pywn(:,{pp},:,:,:),0)
  vtEN_ywn = dim_avg_n_Wrap(fdata[irEN]->vt_pywn(:,{pp},:,:,:),0) - vtClim_ywn
  vtEN_ywn@long_name  = vtClim_ywn@long_name
  delete_VarAtts(vtClim_ywn,"long_name")
  ttClim_ywn = dim_avg_n_Wrap(fdata[irClim]->tt_pywn(:,{pp},:,:,:),0)
  ttEN_ywn = dim_avg_n_Wrap(fdata[irEN]->tt_pywn(:,{pp},:,:,:),0) - ttClim_ywn
  ttEN_ywn@long_name  = ttClim_ywn@long_name
  delete_VarAtts(ttClim_ywn,"long_name")

  vv_ryn = fdata[:]->vv_pyn(:,{pp},:,:)
  vvClim_yn = dim_avg_n_Wrap(fdata[irClim]->vv_pyn(:,{pp},:,:),0)
  vvEN_yn = dim_avg_n_Wrap(fdata[irEN]->vv_pyn(:,{pp},:,:),0) - vvClim_yn
  copy_VarCoords(u_y,vvEN_yn(:,0))
  vvEN_yn@long_name  = str_sub_str(vvClim_yn@long_name,"x hPa",pp+" hPa")
  delete_VarAtts(vvClim_yn,"long_name")
  uv_ryn = fdata[:]->uv_pyn(:,{pp},:,:)
  uvClim_yn = dim_avg_n_Wrap(fdata[irClim]->uv_pyn(:,{pp},:,:),0)
  uvEN_yn = dim_avg_n_Wrap(fdata[irEN]->uv_pyn(:,{pp},:,:),0) - uvClim_yn
  copy_VarCoords(u_y,uvEN_yn(:,0))
  uvEN_yn@long_name  = uvClim_yn@long_name
  delete_VarAtts(uvClim_yn,"long_name")
  uu_ryn = fdata[:]->uu_pyn(:,{pp},:,:)
  uuClim_yn = dim_avg_n_Wrap(fdata[irClim]->uu_pyn(:,{pp},:,:),0)
  uuEN_yn = dim_avg_n_Wrap(fdata[irEN]->uu_pyn(:,{pp},:,:),0) - uuClim_yn
  copy_VarCoords(u_y,uuEN_yn(:,0))
  uuEN_yn@long_name  = uuClim_yn@long_name
  delete_VarAtts(uuClim_yn,"long_name")
  vt_ryn = fdata[:]->vt_pyn(:,{pp},:,:)
  vtClim_yn = dim_avg_n_Wrap(fdata[irClim]->vt_pyn(:,{pp},:,:),0)
  vtEN_yn = dim_avg_n_Wrap(fdata[irEN]->vt_pyn(:,{pp},:,:),0) - vtClim_yn
  copy_VarCoords(u_y,vtEN_yn(:,0))
  vtEN_yn@long_name  = vtClim_yn@long_name
  delete_VarAtts(vtClim_yn,"long_name")
  tt_ryn = fdata[:]->tt_pyn(:,{pp},:,:)
  ttClim_yn = dim_avg_n_Wrap(fdata[irClim]->tt_pyn(:,{pp},:,:),0)
  ttEN_yn = dim_avg_n_Wrap(fdata[irEN]->tt_pyn(:,{pp},:,:),0) - ttClim_yn
  copy_VarCoords(u_y,ttEN_yn(:,0))
  ttEN_yn@long_name  = ttClim_yn@long_name
  delete_VarAtts(ttClim_yn,"long_name")
  vvb_ryn = fdata[:]->vvb_yn(:,:,:)
  vvbClim_yn = dim_avg_n_Wrap(fdata[irClim]->vvb_yn(:,:,:),0)
  vvbEN_yn = dim_avg_n_Wrap(fdata[irEN]->vvb_yn(:,:,:),0) - vvbClim_yn
  copy_VarCoords(u_y,vvbEN_yn(:,0))
  vvbEN_yn@long_name  = vvbClim_yn@long_name
  delete_VarAtts(vvbClim_yn,"long_name")
  uub_ryn = fdata[:]->uub_yn(:,:,:)
  uubClim_yn = dim_avg_n_Wrap(fdata[irClim]->uub_yn(:,:,:),0)
  uubEN_yn = dim_avg_n_Wrap(fdata[irEN]->uub_yn(:,:,:),0) - uubClim_yn
  copy_VarCoords(u_y,uubEN_yn(:,0))
  uubEN_yn@long_name  = uubClim_yn@long_name
  delete_VarAtts(uubClim_yn,"long_name")
  vvs_ryn = fdata[:]->vvs_yn(:,:,:)
  vvsClim_yn = dim_avg_n_Wrap(fdata[irClim]->vvs_yn(:,:,:),0)
  vvsEN_yn = dim_avg_n_Wrap(fdata[irEN]->vvs_yn(:,:,:),0) - vvsClim_yn
  copy_VarCoords(u_y,vvsEN_yn(:,0))
  vvsEN_yn@long_name  = vvsClim_yn@long_name
  delete_VarAtts(vvsClim_yn,"long_name")
  uus_ryn = fdata[:]->uus_yn(:,:,:)
  uusClim_yn = dim_avg_n_Wrap(fdata[irClim]->uus_yn(:,:,:),0)
  uusEN_yn = dim_avg_n_Wrap(fdata[irEN]->uus_yn(:,:,:),0) - uusClim_yn
  copy_VarCoords(u_y,uusEN_yn(:,0))
  uusEN_yn@long_name  = uusClim_yn@long_name
  delete_VarAtts(uusClim_yn,"long_name")

  delete(fdata)
  wallClockElapseTime(wcStrt, "load ens nc",0)
  vvEN_yn@long_name  = str_sub_str(vvEN_yn@long_name,"x hPa",pp+" hPa")
  uvEN_yn@long_name  = str_sub_str(uvEN_yn@long_name,"x hPa",pp+" hPa")
  uuEN_yn@long_name  = str_sub_str(uuEN_yn@long_name,"x hPa",pp+" hPa")
  vtEN_yn@long_name  = str_sub_str(vtEN_yn@long_name,"x hPa",pp+" hPa")
  ttEN_yn@long_name  = str_sub_str(ttEN_yn@long_name,"x hPa",pp+" hPa")
  vvEN_ywn@long_name  = str_sub_str(vvEN_ywn@long_name,"x hPa",pp+" hPa")
  uvEN_ywn@long_name  = str_sub_str(uvEN_ywn@long_name,"x hPa",pp+" hPa")
  uuEN_ywn@long_name  = str_sub_str(uuEN_ywn@long_name,"x hPa",pp+" hPa")
  vtEN_ywn@long_name  = str_sub_str(vtEN_ywn@long_name,"x hPa",pp+" hPa")
  ttEN_ywn@long_name  = str_sub_str(ttEN_ywn@long_name,"x hPa",pp+" hPa")

  if (.not. fileexists("gridfrac.nc")) then
    system("export infile="+FILES1(0)+"; module load matlab; \matlab -softwareopengl -nosplash -nojvm -nodisplay -singleCompThread -r 'gridfrac;exit'")
  end if
  fgrid = addfile("gridfrac.nc","r")  ; from gridfrac.m
  frac_cwn = fgrid->frac_cwn
  phasespeed = fgrid->phasespeed
  delete(fgrid)
  wallClockElapseTime(wcStrt, "load gridfrac",0)

;  phasespeed = -19.5+ispan(0,74,1)
  uvClim_yc = new((/ds(2),dimsizes(phasespeed)/),float)
  uvClim_yc = 0.
  uvEN_yc = uvClim_yc
  if isvar("uv_rywn") then
  uv_ryc = conform_dims((/ds(0),ds(2),dimsizes(phasespeed)/),uvClim_yc,(/1,2/))
  end if
;  do jw = 1,ds(1)-1
;    do in = 1,ds(3)/2
; ; Unit is correct as long as bin size is 1m/s
;      mc = round(-phasespeed(0)+ ( -ds(1)/2.+mod(ds(1)/2.+jw,ds(1)) )*2*get_pi("f")/ds(1)/86400 /in*6.371e6, 3)
;;      mc = round(-phasespeed(0)+ ( -ds(1)/2.+mod(ds(1)/2.+jw,ds(1)) )*2*get_pi("f")/ds(1)/21600 /in*6.371e6, 3)
;      if (mc.ge.0.and.mc.lt.dimsizes(phasespeed)) then
;        uvClim_yc(:,mc) = uvClim_yc(:,mc)+ uvClim_ywn(:,jw,in)
;        uvEN_yc(:,mc) = uvEN_yc(:,mc)+ uvEN_ywn(:,jw,in)
;  if isvar("uv_rywn") then
;        uv_ryc(:,:,mc) = uv_ryc(:,:,mc)+ uv_rywn(:,:,jw,in)
;  end if
;      end if
;    end do
;  end do
  do mc = 0,dimsizes(phasespeed)-1
    uvClim_yc(:,mc) = dim_sum_n(uvClim_ywn *conform(uvClim_ywn,frac_cwn(mc,:,:),(/1,2/)),(/1,2/))
    uvEN_yc(:,mc) = dim_sum_n(uvEN_ywn *conform(uvEN_ywn,frac_cwn(mc,:,:),(/1,2/)),(/1,2/))
  if isvar("uv_rywn") then
    uv_ryc(:,:,mc) = dim_sum_n(uv_rywn *conform(uv_rywn,frac_cwn(mc,:,:),(/2,3/)),(/2,3/))
  end if
  end do
  delete(frac_cwn)
  if isvar("uv_rywn") then
  delete([/uv_rywn/])
  end if
  copy_VarCoords(u_y,uvClim_yc(:,0))
  uvClim_yc!1 = "phasespeed"
  uvClim_yc&phasespeed = phasespeed
;  uvClim_yc@long_name  = "Cospectra of eddy momentum flux at "+pp+" hPa (m s~S~-1~N~)"
;  uvClim_yc@long_name  = "u'v' at "+pp+" hPa (m~S~2~N~ s~S~-2~N~ per m s~S~-1~N~ bin)"
;  uvClim_yc@units      = "day~S~-1"  ;"m s~S~-1~N~ day~S~-1"
  copy_VarCoords(uvClim_yc,uvEN_yc)
  uvEN_yc@long_name  = "u'v' at "+pp+" hPa (m~S~2~N~ s~S~-2~N~ per m s~S~-1~N~ bin)"
  wallClockElapseTime(wcStrt, "collect wave speed",0)

 ; save ds irEN irClim  phasespeed lat u_y uEN_y  uv{Clim_,EN_,_r}yc {uv,vv,uu,vt,tt,vvb,uub,vvs,uus}{Clim_,EN_,_r}yn {uv,vv,uu,vt,tt}{Clim_,EN_}ywn
setfileoption("nc", "Format", "NetCDF4Classic")
;system("rm "+fn_collect)
fout = addfile(fn_collect, "c")
;system("ln -sf "+get_script_prefix_name+"."+caseid+".data.nc "+get_script_prefix_name+".data.nc")

setfileoption(fout, "DefineMode", True)

fAtt           = True
fAtt@caseid    = caseid
fAtt@creation_date = systemfunc("date")
fAtt@dt            = dt
;fAtt@Ly            = Ly
fAtt@Lx            = Lx
fAtt@ds            = ds
fAtt@irEN          = irEN
fAtt@irClim        = irClim
fileattdef(fout, fAtt)
; copy_VarAtts(f, fout)

 dimNames = (/"slice","lat","freq","kx","phasespeed"/)
 dimUnlim = (/False,False,False,False,False/)
 filedimdef(fout, dimNames, (/ds(0),ds(2),ds(1),ds(3)/2+1,dimsizes(phasespeed)/), dimUnlim)
; filedimdef(fout, dimNames, dimsizes(uv_rywn), dimUnlim)
  system("echo finish def dim")

 filevardef(fout, "phasespeed", typeof(phasespeed), "phasespeed")
 filevardef(fout, "lat", typeof(lat), "lat")
 filevardef(fout, "u_y", typeof(u_y), "lat")
 filevardef(fout, "uEN_y", typeof(uEN_y), "lat")
 filevardef(fout, "uvClim_yc", "float", (/"lat","phasespeed"/))
 filevardef(fout, "uvEN_yc", "float", (/"lat","phasespeed"/))
 filevardef(fout, "uv_ryc", "float", (/"slice","lat","phasespeed"/))
 filevardef(fout, "vvClim_yn", "float", (/"lat","kx"/))
 filevardef(fout, "vvEN_yn", "float", (/"lat","kx"/))
 filevardef(fout, "vv_ryn", "float", (/"slice","lat","kx"/))
 filevardef(fout, "uvClim_yn", "float", (/"lat","kx"/))
 filevardef(fout, "uvEN_yn", "float", (/"lat","kx"/))
 filevardef(fout, "uv_ryn", "float", (/"slice","lat","kx"/))
 filevardef(fout, "uuClim_yn", "float", (/"lat","kx"/))
 filevardef(fout, "uuEN_yn", "float", (/"lat","kx"/))
 filevardef(fout, "uu_ryn", "float", (/"slice","lat","kx"/))
 filevardef(fout, "vtClim_yn", "float", (/"lat","kx"/))
 filevardef(fout, "vtEN_yn", "float", (/"lat","kx"/))
 filevardef(fout, "vt_ryn", "float", (/"slice","lat","kx"/))
 filevardef(fout, "ttClim_yn", "float", (/"lat","kx"/))
 filevardef(fout, "ttEN_yn", "float", (/"lat","kx"/))
 filevardef(fout, "tt_ryn", "float", (/"slice","lat","kx"/))
 filevardef(fout, "vvbClim_yn", "float", (/"lat","kx"/))
 filevardef(fout, "vvbEN_yn", "float", (/"lat","kx"/))
 filevardef(fout, "vvb_ryn", "float", (/"slice","lat","kx"/))
 filevardef(fout, "uubClim_yn", "float", (/"lat","kx"/))
 filevardef(fout, "uubEN_yn", "float", (/"lat","kx"/))
 filevardef(fout, "uub_ryn", "float", (/"slice","lat","kx"/))
 filevardef(fout, "vvsClim_yn", "float", (/"lat","kx"/))
 filevardef(fout, "vvsEN_yn", "float", (/"lat","kx"/))
 filevardef(fout, "vvs_ryn", "float", (/"slice","lat","kx"/))
 filevardef(fout, "uusClim_yn", "float", (/"lat","kx"/))
 filevardef(fout, "uusEN_yn", "float", (/"lat","kx"/))
 filevardef(fout, "uus_ryn", "float", (/"slice","lat","kx"/))
 filevardef(fout, "uvClim_ywn", "float", (/"lat","freq","kx"/))
 filevardef(fout, "uvEN_ywn", "float", (/"lat","freq","kx"/))
 filevardef(fout, "vvClim_ywn", "float", (/"lat","freq","kx"/))
 filevardef(fout, "vvEN_ywn", "float", (/"lat","freq","kx"/))
 filevardef(fout, "uuClim_ywn", "float", (/"lat","freq","kx"/))
 filevardef(fout, "uuEN_ywn", "float", (/"lat","freq","kx"/))
 filevardef(fout, "vtClim_ywn", "float", (/"lat","freq","kx"/))
 filevardef(fout, "vtEN_ywn", "float", (/"lat","freq","kx"/))
 filevardef(fout, "ttClim_ywn", "float", (/"lat","freq","kx"/))
 filevardef(fout, "ttEN_ywn", "float", (/"lat","freq","kx"/))
;  system("echo finish def var1")
; filevarattdef(fout,"phasespeed", phasespeed)
 filevarattdef(fout,"lat", lat)
 filevarattdef(fout,"u_y", u_y)
 filevarattdef(fout,"uEN_y", uEN_y)
 filevarattdef(fout,"uvClim_yc", uvClim_yc)
 filevarattdef(fout,"uvEN_yc", uvEN_yc)
 filevarattdef(fout,"uv_ryc", uv_ryc)
 filevarattdef(fout,"vvClim_yn", vvClim_yn)
 filevarattdef(fout,"vvEN_yn", vvEN_yn)
 filevarattdef(fout,"vv_ryn", vv_ryn)
 filevarattdef(fout,"uvClim_yn", uvClim_yn)
 filevarattdef(fout,"uvEN_yn", uvEN_yn)
 filevarattdef(fout,"uv_ryn", uv_ryn)
 filevarattdef(fout,"uuClim_yn", uuClim_yn)
 filevarattdef(fout,"uuEN_yn", uuEN_yn)
 filevarattdef(fout,"uu_ryn", uu_ryn)
 filevarattdef(fout,"vtClim_yn", vtClim_yn)
 filevarattdef(fout,"vtEN_yn", vtEN_yn)
 filevarattdef(fout,"vt_ryn", vt_ryn)
 filevarattdef(fout,"ttClim_yn", ttClim_yn)
 filevarattdef(fout,"ttEN_yn", ttEN_yn)
 filevarattdef(fout,"tt_ryn", tt_ryn)
 filevarattdef(fout,"vvbClim_yn", vvbClim_yn)
 filevarattdef(fout,"vvbEN_yn", vvbEN_yn)
 filevarattdef(fout,"vvb_ryn", vvb_ryn)
 filevarattdef(fout,"uubClim_yn", uubClim_yn)
 filevarattdef(fout,"uubEN_yn", uubEN_yn)
 filevarattdef(fout,"uub_ryn", uub_ryn)
 filevarattdef(fout,"vvsClim_yn", vvsClim_yn)
 filevarattdef(fout,"vvsEN_yn", vvsEN_yn)
 filevarattdef(fout,"vvs_ryn", vvs_ryn)
 filevarattdef(fout,"uusClim_yn", uusClim_yn)
 filevarattdef(fout,"uusEN_yn", uusEN_yn)
 filevarattdef(fout,"uus_ryn", uus_ryn)
 filevarattdef(fout,"uvClim_ywn", uvClim_ywn)
 filevarattdef(fout,"uvEN_ywn", uvEN_ywn)
 filevarattdef(fout,"vvClim_ywn", vvClim_ywn)
 filevarattdef(fout,"vvEN_ywn", vvEN_ywn)
 filevarattdef(fout,"uuClim_ywn", uuClim_ywn)
 filevarattdef(fout,"uuEN_ywn", uuEN_ywn)
 filevarattdef(fout,"vtClim_ywn", vtClim_ywn)
 filevarattdef(fout,"vtEN_ywn", vtEN_ywn)
 filevarattdef(fout,"ttClim_ywn", ttClim_ywn)
 filevarattdef(fout,"ttEN_ywn", ttEN_ywn)

setfileoption(fout, "DefineMode", False)
  system("echo `date` finish def all")
 wallClockElapseTime(wcStrt, "define nc",0)
; wcStrt = systemfunc("date")

 fout->phasespeed = (/phasespeed/)
 fout->lat = (/lat/)
 fout->u_y = (/u_y/)
 fout->uEN_y = (/uEN_y/)
 fout->uvClim_yc = (/uvClim_yc/)
 fout->uvEN_yc = (/uvEN_yc/)
 fout->uv_ryc = (/uv_ryc/)
 fout->vvClim_yn = (/vvClim_yn/)
 fout->vvEN_yn = (/vvEN_yn/)
 fout->vv_ryn = (/vv_ryn/)
 fout->uvClim_yn = (/uvClim_yn/)
 fout->uvEN_yn = (/uvEN_yn/)
 fout->uv_ryn = (/uv_ryn/)
 fout->uuClim_yn = (/uuClim_yn/)
 fout->uuEN_yn = (/uuEN_yn/)
 fout->uu_ryn = (/uu_ryn/)
 fout->vtClim_yn = (/vtClim_yn/)
 fout->vtEN_yn = (/vtEN_yn/)
 fout->vt_ryn = (/vt_ryn/)
 fout->ttClim_yn = (/ttClim_yn/)
 fout->ttEN_yn = (/ttEN_yn/)
 fout->tt_ryn = (/tt_ryn/)
 fout->vvbClim_yn = (/vvbClim_yn/)
 fout->vvbEN_yn = (/vvbEN_yn/)
 fout->vvb_ryn = (/vvb_ryn/)
 fout->uubClim_yn = (/uubClim_yn/)
 fout->uubEN_yn = (/uubEN_yn/)
 fout->uub_ryn = (/uub_ryn/)
 fout->vvsClim_yn = (/vvsClim_yn/)
 fout->vvsEN_yn = (/vvsEN_yn/)
 fout->vvs_ryn = (/vvs_ryn/)
 fout->uusClim_yn = (/uusClim_yn/)
 fout->uusEN_yn = (/uusEN_yn/)
 fout->uus_ryn = (/uus_ryn/)
 fout->uvClim_ywn = (/uvClim_ywn/)
 fout->uvEN_ywn = (/uvEN_ywn/)
 fout->vvClim_ywn = (/vvClim_ywn/)
 fout->vvEN_ywn = (/vvEN_ywn/)
 fout->uuClim_ywn = (/uuClim_ywn/)
 fout->uuEN_ywn = (/uuEN_ywn/)
 fout->vtClim_ywn = (/vtClim_ywn/)
 fout->vtEN_ywn = (/vtEN_ywn/)
 fout->ttClim_ywn = (/ttClim_ywn/)
 fout->ttEN_ywn = (/ttEN_ywn/)
delete(fout) ;close the file
 wallClockElapseTime(wcStrt, "write nc",0)

  else  ;(.not. fileexists(fn_collect))
  fdata = addfile(fn_collect,"r")
;  caseid = fdata@caseid
  dt = fdata@dt
;  Ly = fdata@Ly
  Lx = fdata@Lx
  ds = fdata@ds
  irEN = fdata@irEN
  irClim = fdata@irClim
  phasespeed = fdata->phasespeed
  lat = fdata->lat
  u_y = fdata->u_y
  uEN_y = fdata->uEN_y
  uvClim_yc = fdata->uvClim_yc
  uvEN_yc = fdata->uvEN_yc
;  uv_ryc = fdata->uv_ryc  ;TODO speed concern
  vvClim_yn = fdata->vvClim_yn
  vvEN_yn = fdata->vvEN_yn
  vv_ryn = fdata->vv_ryn
  uvClim_yn = fdata->uvClim_yn
  uvEN_yn = fdata->uvEN_yn
  uv_ryn = fdata->uv_ryn
  uuClim_yn = fdata->uuClim_yn
  uuEN_yn = fdata->uuEN_yn
  uu_ryn = fdata->uu_ryn
  vtClim_yn = fdata->vtClim_yn
  vtEN_yn = fdata->vtEN_yn
  vt_ryn = fdata->vt_ryn
  ttClim_yn = fdata->ttClim_yn
  ttEN_yn = fdata->ttEN_yn
  tt_ryn = fdata->tt_ryn
  vvbClim_yn = fdata->vvbClim_yn
  vvbEN_yn = fdata->vvbEN_yn
  vvb_ryn = fdata->vvb_ryn
  uubClim_yn = fdata->uubClim_yn
  uubEN_yn = fdata->uubEN_yn
  uub_ryn = fdata->uub_ryn
  vvsClim_yn = fdata->vvsClim_yn
  vvsEN_yn = fdata->vvsEN_yn
  vvs_ryn = fdata->vvs_ryn
  uusClim_yn = fdata->uusClim_yn
  uusEN_yn = fdata->uusEN_yn
  uus_ryn = fdata->uus_ryn
  uvClim_ywn = fdata->uvClim_ywn
  uvEN_ywn = fdata->uvEN_ywn
  vvClim_ywn = fdata->vvClim_ywn
  vvEN_ywn = fdata->vvEN_ywn
  uuClim_ywn = fdata->uuClim_ywn
  uuEN_ywn = fdata->uuEN_ywn
  vtClim_ywn = fdata->vtClim_ywn
  vtEN_ywn = fdata->vtEN_ywn
  ttClim_ywn = fdata->ttClim_ywn
  ttEN_ywn = fdata->ttEN_ywn
  delete(fdata)
  wallClockElapseTime(wcStrt, "reload nc",0)
  end if  ;(.not. fileexists(fn_collect))
;  Lx = 2*get_pi("f")*6371 *1000
;  Lx = 32.e6
  Lx := tofloat(Lx)
  fn_lininst = "ensemble-wise/lininst_"+cidCtrl+".nc"
;  if (.not. fileexists(fn_lininst)) then
;    system("export CASENAME="+cidCtrl+"; module load matlab; \matlab -softwareopengl -nosplash -nojvm -nodisplay -singleCompThread -r 'lininst_pchan;exit'")
;  end if
  fdata = addfile(fn_lininst,"r")
  rGrowthClim_yn = fdata->rGrowth_yn
  iGrowthClim_yn = fdata->iGrowth_yn
  delete(fdata)
  fn_lininst = "ensemble-wise/lininst_"+caseid+".nc"
  fdata = addfile(fn_lininst,"r")
  rGrowthEN_yn = fdata->rGrowth_yn
  iGrowthEN_yn = fdata->iGrowth_yn
  delete(fdata)

;TODO  uvClim_yc(:,{-0.5:0.5}) = uvClim_yc@_FillValue
;TODO  uvEN_yc(:,{-0.5:0.5}) = uvEN_yc@_FillValue

  cos2y_yc = conform(uvClim_yc, cos(lat*get_d2r("f"))^2, 0)
  uvyClim_yc = -1/cos2y_yc*center_finite_diff_n(uvClim_yc*cos2y_yc,lat*get_d2r("f")*6.371e6,False,0,0)*86400  ; Eq.9 in Kidston2011
;  uvyClim_yc = wgt_runave_n(uvyClim_yc,(/1./3,1./3,1./3/), 1,1)
  uvyEN_yc = -1/cos2y_yc*center_finite_diff_n(uvEN_yc*cos2y_yc,lat*get_d2r("f")*6.371e6,False,0,0)*86400
;  uvyEN_yc = wgt_runave_n(uvyEN_yc,(/1./3,1./3,1./3/), 1,1)
;  uvyEN_yc = smth9(uvyEN_yc,0.5,0.25,False)
;  uvyEN_yc = smth9(uvyEN_yc,0.5,0.25,False)
;  uvyEN_yc = smth9(uvyEN_yc,0.5,0.25,False)
  if isvar("uv_ryc") then
  uvy_ryc = -1/conform(uv_ryc,cos2y_yc,(/1,2/))*center_finite_diff_n(uv_ryc*conform(uv_ryc,cos2y_yc,(/1,2/)),lat*get_d2r("f")*6.371e6,False,0,1)*86400
;  uvy_ryc = smth9(uvy_ryc,0.5,0.25,False)
;  uvy_ryc = smth9(uvy_ryc,0.5,0.25,False)
;  uvy_ryc = smth9(uvy_ryc,0.5,0.25,False)
  uvyCI_yc = sqrt( (dim_variance_n(uvy_ryc(irClim,:,:),0)*(dimsizes(irClim)-1) + dim_variance_n(uvy_ryc(irEN,:,:),0)*(dimsizes(irEN)-1))/(dimsizes(irClim)+dimsizes(irEN)-2) )*sqrt(1./dimsizes(irClim)+1./dimsizes(irEN)) *cdft_t(0.975,dimsizes(irClim)+dimsizes(irEN)-2)
  end if
  uvEN_yc2 = uvClim_yc+uvEN_yc
  uvEN_yc2 = linint1_n(conform(uvClim_yc,phasespeed,1)-conform(uvClim_yc,uEN_y-u_y,0), uvEN_yc2,False,phasespeed,0,1) -uvClim_yc
  uvyEN_yc2 = -1/cos2y_yc*center_finite_diff_n((uvClim_yc+uvEN_yc)*cos2y_yc,lat*get_d2r("f")*6.371e6,False,0,0)*86400
  uvyEN_yc2 = linint1_n(conform(uvyClim_yc,phasespeed,1)-conform(uvyClim_yc,uEN_y-u_y,0), uvyEN_yc2,False,phasespeed,0,1) -uvyClim_yc
;  uvyEN_yc2 = smth9(uvyEN_yc2,0.5,0.25,False)
;  uvyEN_yc2 = smth9(uvyEN_yc2,0.5,0.25,False)
;  uvyEN_yc2 = smth9(uvyEN_yc2,0.5,0.25,False)
;  wallClockElapseTime(wcStrt, "collect wave speed",0)

  copy_VarCoords(uvClim_yc,uvyClim_yc)
;  uvyClim_yc@long_name  = "Cospectra of eddy momentum flux convergence at "+pp+" hPa (day~S~-1~N~)"
;  uvyClim_yc@long_name  = "Cospectra of eddy momentum flux convergence at "+pp+" hPa (s~S~-1~N~)"
;  uvyClim_yc@long_name  = "~F18~-s~FB~y~N~u'v' at "+pp+" hPa (m s~S~-1~N~day~S~-1~N~ per m s~S~-1~N~ bin)"
  copy_VarCoords(uvyClim_yc,uvyEN_yc)
  uvyEN_yc@long_name  = "~F18~-s~FB~y~N~u'v' at "+pp+" hPa (m s~S~-1~N~day~S~-1~N~ per m s~S~-1~N~ bin)"
;  uvyEN_yc(:,{-0.5:0.5}) = uvyEN_yc@_FillValue
  if isvar("uv_ryc") then
  copy_VarCoords(uvyClim_yc,uvyCI_yc)
;  uvyCI_yc(:,{-0.5:0.5}) = uvyCI_yc@_FillValue
  uvyCI_yc = mask(uvyCI_yc,uvyCI_yc.ne.0,True)
  uvyt_yc = uvyEN_yc/uvyCI_yc
  copy_VarCoords(uvyClim_yc,uvyt_yc)
  uvyt_yc@long_name  = "t/t~B~0.975~N~ for change in ~F18~-s~FB~y~N~u'v' at "+pp+" hPa"
  end if
  copy_VarCoords(uvClim_yc,uvEN_yc2)
  copy_VarCoords(uvyClim_yc,uvyEN_yc2)

/;
  wavelength = 500*(0.5+ispan(0,39,1))
  vvClim_yl = new((/ds(2),dimsizes(wavelength)/),float)
  vvClim_yl = 0.
  vvEN_yl = vvClim_yl
;  vv_ryc = conform_dims((/ds(0),ds(2),dimsizes(wavelength)/),vvClim_yl,(/1,2/))
  jl_yn = round(-0.5+ (2*get_pi("f")*6371*conform(vvClim_yn(:,1:),cos(lat*get_d2r("f")),0)/conform(vvClim_yn(:,1:),ispan(1,ds(3)/2,1),1))/500., 3)
  do jl = 0,dimsizes(wavelength)-1
    vvClim_yl(:,jl) = vvClim_yl(:,jl) + dim_sum_n( where(jl_yn.eq.jl,vvClim_yn(:,1:),0), 1)
    vvEN_yl(:,jl) = vvEN_yl(:,jl) + dim_sum_n( where(jl_yn.eq.jl,vvEN_yn(:,1:),0), 1)
  end do
  vvClim_yl = vvClim_yl/500
  vvEN_yl = vvEN_yl/500
;  do jw = 1,ds(1)-1
;    do in = 1,ds(3)/2
 ; Unit is correct as long as bin size is 1m/s
;      mc = round(-wavelength(0)+ ( -ds(1)/2.+mod(ds(1)/2.+jw,ds(1)) )*2*get_pi("f")/ds(1)/86400 /in*6.371e6, 3)
;      mc = round(-wavelength(0)+ ( -ds(1)/2.+mod(ds(1)/2.+jw,ds(1)) )*2*get_pi("f")/ds(1)/21600 /in*6.371e6, 3)
;      if (mc.ge.0.and.mc.lt.dimsizes(wavelength)) then
;        vvClim_yl(:,mc) = vvClim_yl(:,mc)+ vvClim_ywn(:,jw,in)
;        vvEN_yl(:,mc) = vvEN_yl(:,mc)+ vvEN_ywn(:,jw,in)
;        vv_ryc(:,:,mc) = vv_ryc(:,:,mc)+ vv_rywn(:,:,jw,in)
;      end if
;    end do
;  end do
;  delete([/vv_rywn/])
  copy_VarCoords(u_y,vvClim_yl(:,0))
  vvClim_yl!1 = "wavelength"
  vvClim_yl&wavelength = wavelength
  vvClim_yl@long_name  = "v'v' at "+pp+" hPa (m~S~2~N~ s~S~-2~N~ per 1 km bin)"
;  vvClim_yl@units      = "day~S~-1"  ;"m s~S~-1~N~ day~S~-1"
  copy_VarCoords(vvClim_yl,vvEN_yl)
  vvClim_yl = wgt_runave_n(vvClim_yl,(/1./3,1./3,1./3/), 1,1)
  vvEN_yl = wgt_runave_n(vvEN_yl,(/1./3,1./3,1./3/), 1,1)
;/

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

 ; Ellipse
 ; Ellipse is always well defined. Angle is only defined when a metric is defined
 ; Multivariate Gaussian - [u v]'*inv(cov)*[u v] = const
 ; Eigenvector (of cov or inv(cov)) is major axes
 ; let tan = v/u = (lambda - uu)/uv; lambda = uv*tan+uu
 ; (uv*tan+uu)^2 -(uu+vv)*(uv*tan+uu) +uu*vv -uv^2 =0
 ; tan^2 +(uu-vv)/uv*tan -1 =0
 ; tan2 = 2t/(1-t^2) = 2*uv/(uu-vv)
  angClim_yn = atan2(2*uvClim_yn, uuClim_yn-vvClim_yn)/2/get_d2r("f")
  angEN_yn = -90.+mod(270.+atan2(2*(uvClim_yn+uvEN_yn), uuClim_yn+uuEN_yn-vvClim_yn-vvEN_yn)/2/get_d2r("f") -angClim_yn, 180.)
  copy_VarCoords(u_y,angClim_yn(:,0))
  angClim_yn@long_name  = "Angle from x-axis to major axis at "+pp+" hPa (deg)"
  copy_VarCoords(angClim_yn,angEN_yn)
  angEN_yn@long_name  = "Angle from x-axis to major axis at "+pp+" hPa (deg)"

  ekeClim_yn = uuClim_yn+vvClim_yn
;  ekeClim_yn = mask(ekeClim_yn,ekeClim_yn.ne.0,True)
  ekeEN_yn = uuEN_yn+vvEN_yn
  copy_VarCoords(u_y,ekeClim_yn(:,0))
  copy_VarCoords(ekeClim_yn,ekeEN_yn)
  ekeEN_yn@long_name  = "u'u'+v'v' at "+pp+" hPa (m~S~2~N~ s~S~-2~N~)"
  ekenClim_y = dim_sum_n(ekeClim_yn(:,1:)*conform(vvClim_yn1,ispan(1,ds(3)/2,1),1),1)/dim_sum_n(ekeClim_yn(:,1:),1)
  ekenEN_y = dim_sum_n((ekeClim_yn(:,1:)+ekeEN_yn(:,1:))*conform(vvClim_yn1,ispan(1,ds(3)/2,1),1),1)/dim_sum_n(ekeClim_yn(:,1:)+ekeEN_yn(:,1:),1)

  nn = 8
  yy = minind(uvEN_yn(:,nn))  ; 27, -38.2
  yy = 27
  print("k="+nn+", uvEN_yn most neg: "+yy+" "+lat(yy))
  print("k="+nn+",lat="+lat(yy)+", uvEN_ywn most neg: "+minind(uvEN_ywn(yy,:,nn)))  ; 25
  print(rGrowthClim_yn&kx(maxind(rGrowthClim_yn(yy,:)))+"	"+86400*max(rGrowthClim_yn(yy,:)))
  print(rGrowthClim_yn&kx(maxind(rGrowthClim_yn(yy,0:5)))+"	"+86400*max(rGrowthClim_yn(yy,0:5)))
  print(-86400*(iGrowthEN_yn(yy,:)-iGrowthClim_yn(yy,:))/2/get_pi("f")+"	"+86400*(rGrowthEN_yn(yy,:)-rGrowthClim_yn(yy,:))+"	"+86400*rGrowthClim_yn(yy,:))
  e1Clim_yn = uvClim_yn/mask(ekeClim_yn,ekeClim_yn.ne.0,True)
  e1EN_yn = (uvClim_yn+uvEN_yn)/mask(ekeClim_yn+ekeEN_yn,ekeClim_yn.ne.0,True) -e1Clim_yn
  uv0EN_yn = e1Clim_yn*ekeEN_yn
  uv1EN_yn = e1EN_yn*ekeClim_yn
;  e1EN_yn = log((uvClim_yn+uvEN_yn)/mask(ekeClim_yn+ekeEN_yn,ekeClim_yn.ne.0,True) /e1Clim_yn)
;  replace_ieeenan(e1EN_yn,e1EN_yn@_FillValue,0)
  copy_VarCoords(u_y,e1Clim_yn(:,0))
  copy_VarCoords(e1Clim_yn,e1EN_yn)
  copy_VarCoords(e1Clim_yn,uv0EN_yn)
  copy_VarCoords(e1Clim_yn,uv1EN_yn)
  e1EN_yn@long_name  = "u'v'/(u'u'+v'v')"
  uv0EN_yn@long_name  = "~F5~D~F~ u'v' due to ~F5~D~F~ (u'u'+v'v') (m~S~2~N~ s~S~-2~N~)"
  uv1EN_yn@long_name  = "~F5~D~F~ u'v' due to ~F5~D~F~ u'v'/(u'u'+v'v') (m~S~2~N~ s~S~-2~N~)"

  e2Clim_yn = uvClim_yn/mask(sqrt(uuClim_yn*vvClim_yn),vvClim_yn.ne.0,True)
  e2EN_yn = (uvClim_yn+uvEN_yn)/mask(sqrt((uuClim_yn+uuEN_yn)*(vvClim_yn+vvEN_yn)),vvClim_yn.ne.0,True) -e2Clim_yn
  uv2EN_yn = e2EN_yn*sqrt(uuClim_yn*vvClim_yn)
  copy_VarCoords(u_y,e2Clim_yn(:,0))
  copy_VarCoords(e2Clim_yn,e2EN_yn)
  copy_VarCoords(e2Clim_yn,uv2EN_yn)
  e2EN_yn@long_name  = "u'v'/sqrt(u'u'*v'v')"
  uv2EN_yn@long_name  = "~F5~D~F~ u'v' due to ~F5~D~F~ u'v'/sqrt(u'u'*v'v') (m~S~2~N~ s~S~-2~N~)"
;  e2Clim_ywn = uvClim_ywn/mask(sqrt(uuClim_ywn*vvClim_ywn),vvClim_ywn.ne.0,True)
;  e2EN_ywn = (uvClim_ywn+uvEN_ywn)/mask(sqrt((uuClim_ywn+uuEN_ywn)*(vvClim_ywn+vvEN_ywn)),vvClim_ywn.ne.0,True) -e2Clim_ywn
;  uv2EN_ywn = e2EN_ywn*sqrt(uuClim_ywn*vvClim_ywn)
  vvClim_ywn(:,0,0) = 0.
  vvEN_ywn(:,0,0) = 0.
  uuClim_ywn(:,0,0) = 0.
  uuEN_ywn(:,0,0) = 0.
  uvClim_ywn(:,0,0) = 0.
  uvEN_ywn(:,0,0) = 0.
  ttClim_ywn(:,0,0) = 0.
  ttEN_ywn(:,0,0) = 0.
  vtClim_ywn(:,0,0) = 0.
  vtEN_ywn(:,0,0) = 0.

  et2Clim_yn = vtClim_yn/mask(sqrt(ttClim_yn*vvClim_yn),vvClim_yn.ne.0,True)
  et2EN_yn = (vtClim_yn+vtEN_yn)/mask(sqrt((ttClim_yn+ttEN_yn)*(vvClim_yn+vvEN_yn)),vvClim_yn.ne.0,True) -et2Clim_yn
  vt3EN_yn = et2Clim_yn*(sqrt((ttClim_yn+ttEN_yn)*(vvClim_yn+vvEN_yn))-sqrt(ttClim_yn*vvClim_yn))
  vt2EN_yn = et2EN_yn*sqrt(ttClim_yn*vvClim_yn)
  copy_VarCoords(u_y,et2Clim_yn(:,0))
  copy_VarCoords(et2Clim_yn,et2EN_yn)
  copy_VarCoords(et2Clim_yn,vt3EN_yn)
  copy_VarCoords(et2Clim_yn,vt2EN_yn)
  et2EN_yn@long_name  = "v'T'/sqrt(v'v'*T'T')"
  vt3EN_yn@long_name  = "~F5~D~F~ v'T' due to ~F5~D~F~ sqrt(v'v'*T'T') (K m s~S~-1~N~)"
  vt2EN_yn@long_name  = "~F5~D~F~ v'T' due to ~F5~D~F~ v'T'/sqrt(v'v'*T'T') (K m s~S~-1~N~)"

  type = "pdf"
; type@wkOrientation = "portrait" 
; type@wkPaperWidthF  = 11.2
; type@wkPaperHeightF = 6.3
  wks = gsn_open_wks(type,"ensemble-wise/"+get_script_prefix_name+"."+caseid+"."+pp)  ; TODO
  gsn_define_colormap(wks,"BlWhRe")  ;"BlueWhiteOrangeRed"

; Set some basic resources
  res = True

  res@cnFillOn          = True
  res@cnLinesOn         = False
;  res@cnFillMode        = "CellFill"
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
  res@tmXBLabelFontHeightF  = 0.05
  res@tmXBLabelDeltaF       = -0.5
  res@tmYLLabelFontHeightF  = 0.05
  res@tmYLLabelDeltaF       = -0.5
; res@xyMonoDashPattern = True
 ;res@xyDashPatterns    = 1
; res@xyMarkLineMode    = "Marklines"
; res@xyLineColors      = (/"Blue","Black"/)
 ;res@gsnXRefLine       = 0.5
 ;res@gsnXRefLineDashPattern = 1
;  res@tiXAxisString     = "u~H-16V23~_~V-23~/cos~F33~q~F~ and angular phase speed (m s~S~-1~N~)"
;  res@tiMainFontHeightF = 0.04
  res@tiXAxisFontHeightF= 0.05
  res@tiYAxisFontHeightF= 0.05
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

  gsres = True
;  gsres@tfPolyDrawOrder  = "Draw"

    resP = True
;   resP@txString    = "Time-step precip rate (mm/hr) at "+ times(it)
 ;  resP@txPosYF  = 0.98
    resP@gsnMaximize = True
    resP@gsnPaperOrientation   = "portrait"
    resP@gsnPaperMargin        = 0.
    resP@gsnBoxMargin          = 0.
;    resP@gsnPanelDebug = True
;   resP@gsnPanelLabelBar = True
;   resP@lbLabelFontHeightF    = 0.012
;    resP@gsnPanelFigureStrings = (/"El Nin~H-13V2F35~D~FV-2H3~o"/)
;    resP@gsnPanelFigureStrings = (/"La Nin~H-13V2F35~D~FV-2H3~a"/)
;    resP@gsnPanelFigureStrings = "2009/10" ;(/"2004/05"/) ;"82/83,97/98" ;
;    resP@gsnPanelFigureStrings = (/"N increase"/)
    resP@gsnPanelFigureStringsFontHeightF  = 0.03
    resP@amJust  = "TopRight"
;   resP@pmLabelBarOrthogonalPosF = -0.02   ; panel
    resP@gsnPanelScalePlotIndex   = 0    ; panel 24
;   resP@gsnPanelYWhiteSpacePercent = 2  ; panel 6
;    resP@gsnPanelBottom = 0.04           ; panel 7
;    resP@gsnPanelTop    = 0.935          ; panel 7
;    resP@gsnPanelLeft   = 0.10
;    resP@gsnPanelRight  = 0.8            ; panel 7
 ;  resP@gsnFrame    = False

  plot = new(6, graphic)
  dum  = new(6, graphic)

 ; uv without d/dy
;  res@tiXAxisString     = "u~H-16V23~_~V-23~/cos~F33~q~F~ and angular phase speed (m s~S~-1~N~)"
  res@trYReverse      = True
  res@tiYAxisString     = "u~H-16V23~_~V-23~/cos~F33~q~F~, c/cos~F33~q~F~ (m s~S~-1~N~)"
;  symMinMaxPlt(uvClim_yc,41,False,res)
  symMinMaxPlt(uvEN_yc,41,False,res)
;  plot(0) = gsn_csm_contour(wks,uvClim_yc,res)
;  plot(0) = gsn_csm_contour(wks,uvEN_yc,res)
  plot(0) = gsn_csm_contour(wks,transpose(uvEN_yc),res)
;  symMinMaxPlt(uvEN_yc(1:ds(2)-2,1:dimsizes(phasespeed)-2),9,False,ores)
;  symMinMaxPlt(uvEN_yc,7,False,ores)
  symMinMaxPlt(uvClim_yc,21,False,ores)
  ores@cnMonoLineColor      = True
;  ores@cnLinePalette        = "BlWhRe"
  ores@cnLineThicknessF     = 4.
  ores@gsnContourNegLineDashPattern = 2
;  ores@cnLevelSelectionMode = "AutomaticLevels"
;  ores@cnMaxLevelCount      = 11
;  ores@cnLevelSpacingF = 0.1
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
;  ploto = gsn_csm_contour(wks,uvClim_yc,ores)
  ploto = gsn_csm_contour(wks,transpose(uvClim_yc),ores)
  overlay(plot(0),ploto)
;  draw(plot(0))
;  frame(wks)
  gsres@gsLineThicknessF = 8.
;  dum(0)  = gsn_add_polyline(wks,plot(0),u_y,lat,gsres)
  dum(0)  = gsn_add_polyline(wks,plot(0),lat,u_y,gsres)
  gsres@gsLineColor      = "Green4"
  gsres@gsLineDashPattern = 2
;  dum(1)  = gsn_add_polyline(wks,plot(0),uEN_y,lat,gsres)
  dum(1)  = gsn_add_polyline(wks,plot(0),lat,uEN_y,gsres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

 ; d/dy uv
;  if isvar("uv_ryc") then
;  symMinMaxPlt(uvyEN_yc*where(abs(uvyt_yc).gt.1,1,0),7,False,res)
;  else
  symMinMaxPlt(uvyEN_yc(1:ds(2)-2,1:dimsizes(phasespeed)-2),21,False,res)
;  end if
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
;  plot(0) = gsn_csm_contour(wks,uvyEN_yc,res)
  plot(0) = gsn_csm_contour(wks,transpose(uvyEN_yc),res)
  symMinMaxPlt(uvyClim_yc,21,False,ores)
;  ores@cnLinePalette        = "BlWhRe"
;  ores@cnLineThicknessF     = 2.
  ores@cnMonoLineColor      = True
  ores@gsnContourNegLineDashPattern = 2
;  ores@cnLevelSelectionMode = "AutomaticLevels"
;  ores@cnMaxLevelCount      = 11
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
;  ploto = gsn_csm_contour(wks,uvyClim_yc,ores)
  ploto = gsn_csm_contour(wks,transpose(uvyClim_yc),ores)
  overlay(plot(0),ploto)
  gsres@gsLineThicknessF = 8.
  gsres@gsLineColor      = "Black"
  gsres@gsLineDashPattern = 0
;  dum(0)  = gsn_add_polyline(wks,plot(0),u_y,lat,gsres)
  dum(0)  = gsn_add_polyline(wks,plot(0),lat,u_y,gsres)
  gsres@gsLineColor      = "Green4"
  gsres@gsLineDashPattern = 2
;  dum(1)  = gsn_add_polyline(wks,plot(0),uEN_y,lat,gsres)
  dum(1)  = gsn_add_polyline(wks,plot(0),lat,uEN_y,gsres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

/;
 ; switch contour/shading?
  symMinMaxPlt(uvyClim_yc,81,False,res)
  plot(0) = gsn_csm_contour(wks,uvyClim_yc,res)
;  if isvar("uv_ryc") then
;  symMinMaxPlt(uvyEN_yc*where(abs(uvyt_yc).gt.1,1,0),7,False,ores)
;  else
  symMinMaxPlt(uvyEN_yc(1:ds(2)-2,1:dimsizes(phasespeed)-2),11,False,ores)
;  end if
  ores@cnLinePalette        = "BlueWhiteOrangeRed";"BlWhRe"
;  ores@cnLevelSpacingF = 0.04 ; TODO La Nina
  ores@cnLineThicknessF     = 4.
  ploto = gsn_csm_contour(wks,uvyEN_yc,ores)
  overlay(plot(0),ploto)
  gsres@gsLineColor      = "Black"
  gsres@gsLineDashPattern = 0
  dum(0)  = gsn_add_polyline(wks,plot(0),u_y,lat,gsres)
  gsres@gsLineColor      = "Green4"
  gsres@gsLineDashPattern = 2
  dum(1)  = gsn_add_polyline(wks,plot(0),uEN_y,lat,gsres)
;    resP@gsnPanelFigureStrings = (/"N increase"/)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
;/

  if isvar("uv_ryc") then
;  symMinMaxPlt(uvyCI_yc,41,False,res)
;  res@cnLevelSelectionMode = "AutomaticLevels"
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -1.5
  res@cnMaxLevelValF       = 1.5
  res@cnLevelSpacingF      = 0.5
  plot(0) = gsn_csm_contour(wks,transpose(uvyt_yc),res)
  gsres@gsLineColor      = "Black"
  gsres@gsLineDashPattern = 0
  dum(0)  = gsn_add_polyline(wks,plot(0),lat,u_y,gsres)
  gsres@gsLineColor      = "Green4"
  gsres@gsLineDashPattern = 2
  dum(1)  = gsn_add_polyline(wks,plot(0),lat,uEN_y,gsres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
  end if

;  res@tiXAxisString     = "Relative angular phase speed (m s~S~-1~N~)"
;  res@trXMaxF         = 30.
;  res@sfXArray        = conform(uvClim_yc,phasespeed,1) - conform(uvClim_yc,u_y,0)
;  res@sfYArray        = conform(uvClim_yc,lat,0)
  res@tiYAxisString     = "(c-u~H-16V23~_~V-23~)/cos~F33~q~F~ (m s~S~-1~N~)";"La Nin~H-13V2F35~D~FV-2H3~a"
  res@trYMaxF         = 30.
  res@sfYArray        = transpose(conform(uvClim_yc,phasespeed,1) - conform(uvClim_yc,u_y,0))
  res@sfXArray        = transpose(conform(uvClim_yc,lat,0))
  uvEN_yc2@long_name  = "u'v' at "+pp+" hPa (m~S~2~N~ s~S~-2~N~ per m s~S~-1~N~ bin)"
  symMinMaxPlt(uvEN_yc2,21,False,res)
  res@gsnXAxisIrregular2Linear = False
  res@gsnYAxisIrregular2Linear = False
  res@trGridType      = "TriangularMesh" ;"Curvilinear"
  plot(0) = gsn_csm_contour(wks,transpose(uvEN_yc2),res)
  symMinMaxPlt(uvClim_yc,21,False,ores)
;  ores@cnLinePalette        = "BlueWhiteOrangeRed";"BlWhRe"
;  ores@cnLineThicknessF     = 4.
  ores@trGridType      = res@trGridType
  ores@sfXArray        = res@sfXArray
  ores@sfYArray        = res@sfYArray
  ploto = gsn_csm_contour(wks,transpose(uvClim_yc),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
  delete_VarAtts(ores,(/"trGridType","sfXArray","sfYArray"/))

;  uvyClim_yc@long_name  = "~F18~-s~FB~y~N~u'v' at "+pp+" hPa (m s~S~-1~N~day~S~-1~N~ per m s~S~-1~N~ bin)"
  uvyEN_yc2@long_name  = "~F18~-s~FB~y~N~u'v' at "+pp+" hPa (m s~S~-1~N~day~S~-1~N~ per m s~S~-1~N~ bin)"
;  symMinMaxPlt(uvyClim_yc,81,False,res)
  symMinMaxPlt(uvyEN_yc2,21,False,res)
  res@gsnXAxisIrregular2Linear = False
  res@gsnYAxisIrregular2Linear = False
  res@trGridType      = "TriangularMesh" ;"Curvilinear"
;  plot(0) = gsn_csm_contour(wks,uvyClim_yc,res)
  plot(0) = gsn_csm_contour(wks,transpose(uvyEN_yc2),res)
;  symMinMaxPlt(uvyEN_yc2(1:ds(2)-2,1:dimsizes(phasespeed)-2),11,False,ores)
  symMinMaxPlt(uvyClim_yc,21,False,ores)
;  ores@cnLinePalette        = "BlueWhiteOrangeRed";"BlWhRe"
;  ores@cnLineThicknessF     = 4.
  ores@trGridType      = res@trGridType
  ores@sfXArray        = res@sfXArray
  ores@sfYArray        = res@sfYArray
;  ploto = gsn_csm_contour(wks,uvyEN_yc2,ores)
  ploto = gsn_csm_contour(wks,transpose(uvyClim_yc),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
  delete_VarAtts(ores,(/"trGridType","sfXArray","sfYArray"/))

;  uvyEN_yc2@long_name  = "Change in ~F18~-s~FB~y~N~u'v' at "+pp+" hPa (m s~S~-1~N~day~S~-1~N~ per m s~S~-1~N~ bin)"
;  symMinMaxPlt(uvyEN_yc2,21,False,res)
;  plot(0) = gsn_csm_contour(wks,uvyEN_yc2,res)
;    gsn_panel(wks,plot(:0),(/1,1/),resP)
  delete_VarAtts(res,(/"trGridType","sfXArray","sfYArray"/))
;  delete_VarAtts(res,(/"trXMaxF"/))
  delete_VarAtts(res,(/"trYReverse"/))


/;
 ; AOFD 2019
  vvEN_yl = mask(vvEN_yl,vvEN_yl.ne.0,True)
  symMinMaxPlt(vvEN_yl,9,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYMaxF         = max(wavelength)
  res@trYMinF         = min(wavelength)
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavelength (km)"
  plot(0) = gsn_csm_contour(wks,transpose(vvEN_yl),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 6
  ores@cnLineThicknessF     = 4.
  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(vvClim_yl),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
;/

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
;/

;  vvEN_yn = mask(vvEN_yn,vvEN_yn.ne.0,True)
  symMinMaxPlt(vvEN_yn,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(vvEN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 6
  ores@cnLineThicknessF     = 4.
  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
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

  symMinMaxPlt(uuEN_yn,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(uuEN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 6
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(uuClim_yn),ores)
  overlay(plot(0),ploto)
;  gsres@gsLineThicknessF = 8.
  gsres@gsLineColor      = "Black"
  gsres@gsLineDashPattern = 0
;  dum(0)  = gsn_add_polyline(wks,plot(0),lat,2*get_pi("f")*6371*cos(lat*get_d2r("f"))/uulClim_y,gsres)
  dum(0)  = gsn_add_polyline(wks,plot(0),lat,uunClim_y,gsres)
  gsres@gsLineColor      = "Green4"
  gsres@gsLineDashPattern = 2
;  dum(1)  = gsn_add_polyline(wks,plot(0),lat,2*get_pi("f")*6371*cos(lat*get_d2r("f"))/uulEN_y,gsres)
  dum(1)  = gsn_add_polyline(wks,plot(0),lat,uunEN_y,gsres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  symMinMaxPlt(ekeEN_yn,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(ekeEN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 7
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(ekeClim_yn),ores)
  overlay(plot(0),ploto)
;  gsres@gsLineThicknessF = 8.
  gsres@gsLineColor      = "Black"
  gsres@gsLineDashPattern = 0
;  dum(0)  = gsn_add_polyline(wks,plot(0),lat,2*get_pi("f")*6371*cos(lat*get_d2r("f"))/ekelClim_y,gsres)
  dum(0)  = gsn_add_polyline(wks,plot(0),lat,ekenClim_y,gsres)
  gsres@gsLineColor      = "Green4"
  gsres@gsLineDashPattern = 2
;  dum(1)  = gsn_add_polyline(wks,plot(0),lat,2*get_pi("f")*6371*cos(lat*get_d2r("f"))/ekelEN_y,gsres)
  dum(1)  = gsn_add_polyline(wks,plot(0),lat,ekenEN_y,gsres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  symMinMaxPlt(uvEN_yn,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(uvEN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 21
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(uvClim_yn),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  symMinMaxPlt(uv0EN_yn,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(uv0EN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 21
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(uvClim_yn),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

/;
  symMinMaxPlt(uv1EN_yn,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(uv1EN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 21
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(uvClim_yn),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
;/

  symMinMaxPlt(uv2EN_yn,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(uv2EN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 21
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(uvClim_yn),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

/;
  symMinMaxPlt(e1EN_yn(2:ds(2)-2,:15),51,False,res)
;  res@cnLevelSelectionMode = "ManualLevels"
;  res@cnMinLevelValF       = -0.75
;  res@cnMaxLevelValF       = 0.75
;  res@cnLevelSpacingF      = 0.03
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(e1EN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 21
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(e1Clim_yn),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  symMinMaxPlt(e2EN_yn(2:ds(2)-2,:15),51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(e2EN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 21
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(e2Clim_yn),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -75.
  res@cnMaxLevelValF       = 75.
  res@cnLevelSpacingF      = 15.
  res@cnFillPalette     = "MPL_hsv"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(angClim_yn),res)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
;/
/;
  symMinMaxPlt(angEN_yn,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(angEN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 6
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(angClim_yn),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
;/

  if isvar("vv_ryn") then
/;
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -1.5
  res@cnMaxLevelValF       = 1.5
  res@cnLevelSpacingF      = 0.5
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  plot(0) = gsn_csm_contour(wks,transpose(vvt_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLineColor          = "Green4"
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 6
  ores@cnLineThicknessF     = 4.
  ores@gsnLeftString        = ""
  ploto = gsn_csm_contour(wks,transpose(vvClim_yn),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
;/

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

 ; switch contour/shading?
/;
  res@cnFillPalette     = "WhBlGrYeRe"
  res@cnLevelSelectionMode = "AutomaticLevels"
  delete_VarAtts(res,(/"cnLevelSpacingF"/))
;  res@trGridType      = "Curvilinear"
;  res@sfYArray        = conform(vvClim_yn(:,1:),lat,0)
;  res@sfXArray        = 2*get_pi("f")*6371*cos(res@sfYArray*get_d2r("f"))/conform(vvClim_yn(:,1:),ispan(1,ds(3)/2,1),1)
  res@tiXAxisString     = "Wavelength (km)"
;  plot(0) = gsn_csm_contour(wks,(/vvClim_yn(:,1:)/),res)
  plot(0) = gsn_csm_contour(wks,vvClim_yl,res)

;  symMinMaxPlt(vvEN_yl(1:ds(2)-2,1:dimsizes(phasespeed)-2),9,False,ores)
;  symMinMaxPlt(vvEN_yl*where(abs(uvyt_yc).gt.1,1,0),7,False,ores)
  symMinMaxPlt(vvEN_yl,7,False,ores)
  ores@cnMonoLineColor      = False
  ores@cnLinePalette        = "BlWhRe"
;  ores@cnLevelSpacingF = 0.04
;  ores@sfXArray  = res@sfXArray
;  ores@sfYArray  = res@sfYArray
  ploto = gsn_csm_contour(wks,vvEN_yl,ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
;/

  res@sfXArray        = -0.5+tofloat(ispan(0,ds(1),1))/ds(1)
  ores@sfXArray        = -0.5+tofloat(ispan(0,ds(1),1))/ds(1)
  res@gsnAddCyclic      = True
  ores@gsnAddCyclic      = True
  res@gsnRightString    = get_lat_labels(round(lat(yy),3))
  ores@cnInfoLabelJust       = "TopLeft"
  ores@cnInfoLabelParallelPosF= 0.
  res@pmLabelBarOrthogonalPosF = 0.3

  symMinMaxPlt(vvEN_ywn(yy,:,:),51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = "Freq (cycle/day)"
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(vvEN_ywn(yy,ndtooned((/ispan(ds(1)/2,ds(1)-1,1),ispan(0,ds(1)/2-1,1)/)),:)),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 6
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(vvClim_ywn(yy,ndtooned((/ispan(ds(1)/2,ds(1)-1,1),ispan(0,ds(1)/2-1,1)/)),:)),ores)
  overlay(plot(0),ploto)
  gsres@gsLineColor      = "Black"
  gsres@gsLineDashPattern = 0
  dum(0)  = gsn_add_polyline(wks,plot(0),-iGrowthClim_yn(yy,:)/2/get_pi("f")*86400,iGrowthClim_yn&kx,gsres)
  gsres@gsLineColor      = "Green4"
  gsres@gsLineDashPattern = 2
  dum(1)  = gsn_add_polyline(wks,plot(0),-iGrowthEN_yn(yy,:)/2/get_pi("f")*86400,iGrowthEN_yn&kx,gsres)
  yxres  = True
;  yxres@vpXF            = 0.02
  yxres@vpWidthF        = res@vpHeightF*0.4
;  yxres@vpYF            = res@vpYF
  yxres@vpHeightF       = res@vpHeightF
  yxres@gsnDraw   = False
  yxres@gsnFrame  = False
;  yxres@tmXBLabelsOn    = False
  yxres@tmXBMaxTicks    = 4
  yxres@tmXBFormat      = "#f"
  yxres@tmXBLabelDeltaF = -0.5
  yxres@tmYLLabelsOn    = False
  yxres@tiYAxisString = ""
  yxres@tiXAxisString = "GrowthRate ~C~(1/day)"
;  yxres@tiXAxisString = "xxx m~S~2~N~ s~S~-2~N~"
  yxres@tiDeltaF      = 0.
;  yxres@trYLog        = True
  yxres@trYReverse    = res@trYReverse
  yxres@trYMinF       = 0.
  yxres@trYMaxF       = res@trYMaxF
;  yxres@gsnXRefLine   = 0
;  yxres@xyLabelMode      = "Custom"
;  yxres@xyExplicitLabels = (/"uv","vv","uu"/)
  ploty = gsn_csm_xy(wks,(/rGrowthClim_yn(yy,:),rGrowthEN_yn(yy,:)/)*86400,rGrowthClim_yn&kx,yxres)
  dum4 = gsn_attach_plots(plot(0),ploty,res,False)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  symMinMaxPlt(uuEN_ywn(yy,:,:),51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = "Freq (cycle/day)"
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(uuEN_ywn(yy,ndtooned((/ispan(ds(1)/2,ds(1)-1,1),ispan(0,ds(1)/2-1,1)/)),:)),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 6
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(uuClim_ywn(yy,ndtooned((/ispan(ds(1)/2,ds(1)-1,1),ispan(0,ds(1)/2-1,1)/)),:)),ores)
  overlay(plot(0),ploto)
  gsres@gsLineColor      = "Black"
  gsres@gsLineDashPattern = 0
  dum(0)  = gsn_add_polyline(wks,plot(0),-iGrowthClim_yn(yy,:)/2/get_pi("f")*86400,iGrowthClim_yn&kx,gsres)
  gsres@gsLineColor      = "Green4"
  gsres@gsLineDashPattern = 2
  dum(1)  = gsn_add_polyline(wks,plot(0),-iGrowthEN_yn(yy,:)/2/get_pi("f")*86400,iGrowthEN_yn&kx,gsres)
  dum4 = gsn_attach_plots(plot(0),ploty,res,False)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  symMinMaxPlt(uvEN_ywn(yy,:,:),51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = "Freq (cycle/day)"
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(uvEN_ywn(yy,ndtooned((/ispan(ds(1)/2,ds(1)-1,1),ispan(0,ds(1)/2-1,1)/)),:)),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 6
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(uvClim_ywn(yy,ndtooned((/ispan(ds(1)/2,ds(1)-1,1),ispan(0,ds(1)/2-1,1)/)),:)),ores)
  overlay(plot(0),ploto)
  gsres@gsLineColor      = "Black"
  gsres@gsLineDashPattern = 0
  dum(0)  = gsn_add_polyline(wks,plot(0),-iGrowthClim_yn(yy,:)/2/get_pi("f")*86400,iGrowthClim_yn&kx,gsres)
  gsres@gsLineColor      = "Green4"
  gsres@gsLineDashPattern = 2
  dum(1)  = gsn_add_polyline(wks,plot(0),-iGrowthEN_yn(yy,:)/2/get_pi("f")*86400,iGrowthEN_yn&kx,gsres)
  dum4 = gsn_attach_plots(plot(0),ploty,res,False)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
  delete_VarAtts(res,(/"sfXArray","gsnAddCyclic","gsnRightString" ,"pmLabelBarOrthogonalPosF"/))
  delete_VarAtts(ores,(/"sfXArray","gsnAddCyclic"/))

  symMinMaxPlt(ttEN_yn,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(ttEN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 6
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(ttClim_yn),ores)
  overlay(plot(0),ploto)
;  gsres@gsLineThicknessF = 8.
;  gsres@gsLineColor      = "Black"
;  gsres@gsLineDashPattern = 0
;  dum(0)  = gsn_add_polyline(wks,plot(0),lat,2*get_pi("f")*6371*cos(lat*get_d2r("f"))/ttlClim_y,gsres)
;  dum(0)  = gsn_add_polyline(wks,plot(0),lat,ttnClim_y,gsres)
;  gsres@gsLineColor      = "Green4"
;  gsres@gsLineDashPattern = 2
;  dum(1)  = gsn_add_polyline(wks,plot(0),lat,2*get_pi("f")*6371*cos(lat*get_d2r("f"))/ttlEN_y,gsres)
;  dum(1)  = gsn_add_polyline(wks,plot(0),lat,ttnEN_y,gsres)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  symMinMaxPlt(vtEN_yn,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(vtEN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 21
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(vtClim_yn),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  symMinMaxPlt(vt3EN_yn,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(vt3EN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 21
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(vtClim_yn),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  symMinMaxPlt(vt2EN_yn,51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(vt2EN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 21
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(vtClim_yn),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  symMinMaxPlt(et2EN_yn(2:ds(2)-2,:15),51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = ""
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(et2EN_yn),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 21
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(et2Clim_yn),ores)
  overlay(plot(0),ploto)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  res@sfXArray        = -0.5+tofloat(ispan(0,ds(1),1))/ds(1)
  ores@sfXArray        = -0.5+tofloat(ispan(0,ds(1),1))/ds(1)
  res@gsnAddCyclic      = True
  ores@gsnAddCyclic      = True
  res@gsnRightString    = get_lat_labels(round(lat(yy),3))
  ores@cnInfoLabelJust       = "TopLeft"
  ores@cnInfoLabelParallelPosF= 0.
  res@pmLabelBarOrthogonalPosF = 0.3

  symMinMaxPlt(ttEN_ywn(yy,:,:),51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = "Freq (cycle/day)"
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(ttEN_ywn(yy,ndtooned((/ispan(ds(1)/2,ds(1)-1,1),ispan(0,ds(1)/2-1,1)/)),:)),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 6
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(ttClim_ywn(yy,ndtooned((/ispan(ds(1)/2,ds(1)-1,1),ispan(0,ds(1)/2-1,1)/)),:)),ores)
  overlay(plot(0),ploto)
  gsres@gsLineColor      = "Black"
  gsres@gsLineDashPattern = 0
  dum(0)  = gsn_add_polyline(wks,plot(0),-iGrowthClim_yn(yy,:)/2/get_pi("f")*86400,iGrowthClim_yn&kx,gsres)
  gsres@gsLineColor      = "Green4"
  gsres@gsLineDashPattern = 2
  dum(1)  = gsn_add_polyline(wks,plot(0),-iGrowthEN_yn(yy,:)/2/get_pi("f")*86400,iGrowthEN_yn&kx,gsres)
  yxres  = True
;  yxres@vpXF            = 0.02
  yxres@vpWidthF        = res@vpHeightF*0.4
;  yxres@vpYF            = res@vpYF
  yxres@vpHeightF       = res@vpHeightF
  yxres@gsnDraw   = False
  yxres@gsnFrame  = False
;  yxres@tmXBLabelsOn    = False
  yxres@tmXBMaxTicks    = 4
  yxres@tmXBFormat      = "#f"
  yxres@tmXBLabelDeltaF = -0.5
  yxres@tmYLLabelsOn    = False
  yxres@tiYAxisString = ""
  yxres@tiXAxisString = "GrowthRate ~C~(1/day)"
;  yxres@tiXAxisString = "xxx m~S~2~N~ s~S~-2~N~"
  yxres@tiDeltaF      = 0.
;  yxres@trYLog        = True
  yxres@trYReverse    = res@trYReverse
  yxres@trYMinF       = 0.
  yxres@trYMaxF       = res@trYMaxF
;  yxres@gsnXRefLine   = 0
;  yxres@xyLabelMode      = "Custom"
;  yxres@xyExplicitLabels = (/"uv","vv","uu"/)
  ploty = gsn_csm_xy(wks,(/rGrowthClim_yn(yy,:),rGrowthEN_yn(yy,:)/)*86400,rGrowthClim_yn&kx,yxres)
  dum4 = gsn_attach_plots(plot(0),ploty,res,False)
    gsn_panel(wks,plot(:0),(/1,1/),resP)

  symMinMaxPlt(vtEN_ywn(yy,:,:),51,False,res)
  res@cnFillPalette     = "BlWhRe";"BlueWhiteOrangeRed"
  res@trYReverse      = True
  res@trYMaxF         = 15.
  res@tiXAxisString     = "Freq (cycle/day)"
  res@tiYAxisString     = "Wavenumber"
  plot(0) = gsn_csm_contour(wks,transpose(vtEN_ywn(yy,ndtooned((/ispan(ds(1)/2,ds(1)-1,1),ispan(0,ds(1)/2-1,1)/)),:)),res)
  ores@cnMonoLineColor      = True
  ores@cnLevelSelectionMode = "AutomaticLevels"
  ores@cnMaxLevelCount      = 6
  ores@cnLineThicknessF     = 4.
;  delete_VarAtts(ores,(/"cnLevelSpacingF"/))
  ploto = gsn_csm_contour(wks,transpose(vtClim_ywn(yy,ndtooned((/ispan(ds(1)/2,ds(1)-1,1),ispan(0,ds(1)/2-1,1)/)),:)),ores)
  overlay(plot(0),ploto)
  gsres@gsLineColor      = "Black"
  gsres@gsLineDashPattern = 0
  dum(0)  = gsn_add_polyline(wks,plot(0),-iGrowthClim_yn(yy,:)/2/get_pi("f")*86400,iGrowthClim_yn&kx,gsres)
  gsres@gsLineColor      = "Green4"
  gsres@gsLineDashPattern = 2
  dum(1)  = gsn_add_polyline(wks,plot(0),-iGrowthEN_yn(yy,:)/2/get_pi("f")*86400,iGrowthEN_yn&kx,gsres)
  dum4 = gsn_attach_plots(plot(0),ploty,res,False)
    gsn_panel(wks,plot(:0),(/1,1/),resP)
  delete_VarAtts(res,(/"sfXArray","gsnAddCyclic","gsnRightString" ,"pmLabelBarOrthogonalPosF"/))
  delete_VarAtts(ores,(/"sfXArray","gsnAddCyclic"/))

end
exit

 ;:'<,'> s?filevardef(\([^,]*\), "\([^"]*\)",.*$?filevarattdef(\1,"\2", \2)? |noh
 ;:'<,'> s?filevardef(\([^,]*\), "\([^"]*\)",.*$?\1->\2 = (/\2/)? |noh
 ;:'<,'> s?filevardef(\([^,]*\), "\([^"]*\)",.*$?\2 = fdata->\2? |noh
 ;:'<,'> g/^[^;]*contour/s/wks,\([^,]*\)/wks,transpose(\1)/gce |noh
 ;:'<,'> g/^[^;]*polyline/s/\([^,]*\),\([^,]*\),gsres/\2,\1,gsres/gce |noh

; vim: set fdm=marker foldmarker=/;,;/:

