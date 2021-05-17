; latitude-phase speed spectrum
;TODO sbatch --account=kuang_lab -p huce_intel -J cospectrum --array=1-1 -n 1 -t 1440 --mem=20000 -o "slurm.%a" --wrap='ncl -Q nc.cospectrum.ncl  >& slurm-${SLURM_JOB_NAME}-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}'
begin
;  pp = (/800.,300./)
  wcStrt = systemfunc("date")
  tid = stringtoint(getenv("tid")) ;TODO

 ; FMS
 ; \ls -1v analysis/var_a[15]*_*.nc > cy.cospectrum.filelist
;  FILES1 = systemfunc("cat "+get_script_prefix_name+".filelist") ;systemfunc (" ls -1 var_a[14]*.nc")
  caseid = getenv("CASENAME")
  if (ismissing(caseid)) then
    print(caseid)
    status_exit(1)
  end if
;  fn_save = "cy.cospectrum."+caseid+".data.nc"
;  system("mkdir -p analysis")
 ; NCO -n option requires fixed width numbers. http://nco.sourceforge.net/nco.html#Specifying-Input-Files
  fn_save = "analysis/"+get_script_prefix_name+"_"+sprinti("%03d",tid)+".nc"  ;tostring(tid)
;  if (.not.ismissing(systemfunc("find "+fn_save1+" -newer "+get_script_name))) then
;    print("File exists: "+fn_save1)
;    status_exit(0)
;  end if
  if (systemfunc("cmp -s ~/script/fms/nc.1cospectrum.ncl "+get_script_name+"; echo $?").ne."0") then
    print("ncl file different!")
    status_exit(1)
  end if

 ; obtain grid only
  Run = ispan((tid-1)+1,tid,1)
  fv = addfile("analysis/var_"+caseid+"_"+sprinti("%03d",Run(0))+".nc","r")
  ds = getfilevardimsizes(fv,"u")
  delete(fv)  ; no list
  fn_fms = addfile("history/day"+sprinti("%05d",1000+100*Run(0))+"h00.1xday.nc","r")  ; no list
  ds_org = toint(getfilevardimsizes(fn_fms,"ucomp"))
  dt = (fn_fms->time(1)-fn_fms->time(0))*86400.
;  Ly = fn_fms@Ly
  Lx = 6.371e6*2*get_pi("d")
;  lev_p = fn_fms&lev_p
;  lat = fn_fms&lat
  delete(fn_fms)  ; no list

;  nfiles = toint(systemfunc("ls -1dv analysis/var_"+caseid+"_???.nc |wc -l"))
;  FILES1 = "analysis/var_"+caseid+"_"+sprinti("%03d",ispan(1,nfiles,1))+".nc"
;print(FILES1)
;  ds = (/dimsizes(Run)*ds_org(0),dimsizes(lev_p),ds_org(2),ds_org(3)/)
  ds = ds*(/dimsizes(Run),1,1,1/)
  ds_org = ds_org*(/dimsizes(Run),1,1,1/)

  u = new(ds,float)
  v = new(ds,float)
  t = new(ds,float)
  uo = new(ds_org,float)
  vo = new(ds_org,float)
;  z = new(ds,float)
;  z1 = new(ds_org+(/0,0,0,0/),float)
; wallClockElapseTime(wcStrt, "metadata",0)
; wcStrt = systemfunc("date")

do n = 0,dimsizes(Run)-1
  system("echo -ne "+str_get_dq+"\r"+Run(n)+str_get_dq)

  fv = addfile("analysis/var_"+caseid+"_"+sprinti("%03d",Run(n))+".nc","r")
  fn_fms = addfile("history/day"+sprinti("%05d",1000+100*Run(n))+"h00.1xday.nc","r")  ; no list
  u (n*ds_org(0) : (n+1)*ds_org(0)-1,:,:,:) = fv->u
  v (n*ds_org(0) : (n+1)*ds_org(0)-1,:,:,:) = fv->v
  t (n*ds_org(0) : (n+1)*ds_org(0)-1,:,:,:) = fv->t
  uo(n*ds_org(0) : (n+1)*ds_org(0)-1,:,:,:) = fn_fms->ucomp
  vo(n*ds_org(0) : (n+1)*ds_org(0)-1,:,:,:) = fn_fms->vcomp
  delete(fv)
  delete(fn_fms)
end do

  lev_p = u&lev_p
  lat = u&lat
  wallClockElapseTime(wcStrt, "load",0)

  u := reshape(u,(/1,ds(0),ds(1),ds(2),ds(3)/))
;  u := conform_dims((/1,ds(0),ds(1),ds(2),ds(3)/),u,(/1,2,3,4/))
  copy_VarMeta(v,u(0,:,:,:,:))
;printVarSummary(u)
  v := reshape(v,(/1,ds(0),ds(1),ds(2),ds(3)/))
  t := reshape(t,(/1,ds(0),ds(1),ds(2),ds(3)/))
  ds := toint(dimsizes(u))  ;yr,t,p,y,x
  ds := (/ds(0)*2,ds(1),ds(2),ds(3)/2,ds(4)/)
  u2 = array_append_record(u(:,:,:,0:ds(3)-1,:),u(:,:,:,2*ds(3)-1:ds(3),:),0)
  v2 = array_append_record(v(:,:,:,0:ds(3)-1,:),-v(:,:,:,2*ds(3)-1:ds(3),:),0)
  t2 = array_append_record(t(:,:,:,0:ds(3)-1,:),t(:,:,:,2*ds(3)-1:ds(3),:),0)
  u := u2
  v := v2
  t := t2
  delete([/u2,v2,t2/])
  lat := u&lat

  u_py = dim_avg_n_Wrap(dim_avg_n_Wrap(u,4),(/0,1/))/conform_dims(ds(2:3),cos(lat*get_d2r("f")),1)
;  u_py = dim_avg_n_Wrap(dim_avg_n_Wrap( reshape(u(irClim,:,:,:),(/dimsizes(irClim),ds(1),ds(3),ds(4)/)) ,3),(/0,1/))/cos(lat*get_d2r("f"))
;  uEN_y = dim_avg_n_Wrap(dim_avg_n_Wrap( reshape(u(irEN,:,:,:),(/dimsizes(irEN),ds(1),ds(3),ds(4)/)) ,3),(/0,1/))/cos(lat*get_d2r("f"))
  copy_VarMeta(u(0,0,:,:,0),u_py)

;  u = dtrend_quadratic_msg_n(u,True,False,1)  ; why?
;  v = dtrend_quadratic_msg_n(v,True,False,1)
  u_rpyiwn = new((/ds(0),ds(2),ds(3),2,ds(1),ds(4)/2+1/),float)  ; yr,p,y,Re/Im,w,kx
  v_rpyiwn = new((/ds(0),ds(2),ds(3),2,ds(1),ds(4)/2+1/),float)
  t_rpyiwn = new((/ds(0),ds(2),ds(3),2,ds(1),ds(4)/2+1/),float)
  do ir = 0,ds(0)-1
    do kp = 0,ds(2)-1
      do jy = 0,ds(3)-1
        u_rpyiwn(ir,kp,jy,:,:,:) = fft2df(u(ir,:,kp,jy,:))
        v_rpyiwn(ir,kp,jy,:,:,:) = fft2df(v(ir,:,kp,jy,:))
        t_rpyiwn(ir,kp,jy,:,:,:) = fft2df(t(ir,:,kp,jy,:))
      end do
    end do
  end do
;  vv_rpywn = 0.5*dim_sum_n(v_rpyiwn*v_rpyiwn,3)  ; i sum
;  vv_pywn  = dim_avg_n(vv_rpywn,0)  ; yr avg
;print("global max: "+max(vv_pywn(1,26,:,:)))
;print("k=6 max: "+max(vv_pywn(1,26,:,6)))
;print("k=6 maxind: "+maxind(vv_pywn(1,26,:,6)))  ; 5
;exit
  uv_rpywn = 0.5*dim_sum_n(u_rpyiwn*v_rpyiwn,3)  ; i sum
  uv_pywn  = dim_avg_n(uv_rpywn,0)  ; yr avg
;  uvClim_pywn  = dim_avg_n( reshape(uv_rpywn(irClim,:,:,:),(/dimsizes(irClim),ds(3),ds(1),ds(4)/2+1/)) ,0)
;  uvEN_pywn  = dim_avg_n( reshape(uv_rpywn(irEN,:,:,:),(/dimsizes(irEN),ds(3),ds(1),ds(4)/2+1/)) ,0) - uvClim_pywn
  vv_rpywn = 0.5*dim_sum_n(v_rpyiwn*v_rpyiwn,3)  ; i sum
  vv_pywn  = dim_avg_n(vv_rpywn,0)  ; yr avg
  uu_rpywn = 0.5*dim_sum_n(u_rpyiwn*u_rpyiwn,3)  ; i sum
  uu_pywn  = dim_avg_n(uu_rpywn,0)  ; yr avg
  vt_rpywn = 0.5*dim_sum_n(v_rpyiwn*t_rpyiwn,3)  ; i sum
  vt_pywn  = dim_avg_n(vt_rpywn,0)  ; yr avg
  tt_rpywn = 0.5*dim_sum_n(t_rpyiwn*t_rpyiwn,3)  ; i sum
  tt_pywn  = dim_avg_n(tt_rpywn,0)  ; yr avg
  delete([/u_rpyiwn,v_rpyiwn,uv_rpywn,vv_rpywn,uu_rpywn,vt_rpywn,tt_rpywn/])
  vv_pywn@long_name  = "v'v' at x hPa (m~S~2~N~ s~S~-2~N~)"
  uv_pywn@long_name  = "u'v' at x hPa (m~S~2~N~ s~S~-2~N~)"
  uu_pywn@long_name  = "u'u' at x hPa (m~S~2~N~ s~S~-2~N~)"
  vt_pywn@long_name  = "v'T' at x hPa (K m s~S~-1~N~)"
  tt_pywn@long_name  = "T'T' at x hPa (K~S~2~N~)"
;  vv_rpyn = 0.5*dim_sum_n(v_rpyiwn^2,(/2,3/))
  v_irtpyn = new((/2,ds(0),ds(1),ds(2),ds(3),ds(4)/2+1/),float)
  v_irtpyn = 0.
  u_irtpyn = v_irtpyn
  t_irtpyn = v_irtpyn
  v_irtpyn(:,:,:,:,:,1:) = ezfftf(v)
  u_irtpyn(:,:,:,:,:,1:) = ezfftf(u)
  t_irtpyn(:,:,:,:,:,1:) = ezfftf(t)
  vv_rpyn = 0.5*dim_sum_n(dim_avg_n(v_irtpyn^2,2),0)  ; t avg, i sum
  vv_pyn  = dim_avg_n(vv_rpyn,0)
  uv_rpyn = 0.5*dim_sum_n(dim_avg_n(u_irtpyn*v_irtpyn,2),0)  ; t avg, i sum
  uv_pyn  = dim_avg_n(uv_rpyn,0)
  uu_rpyn = 0.5*dim_sum_n(dim_avg_n(u_irtpyn^2,2),0)  ; t avg, i sum
  uu_pyn  = dim_avg_n(uu_rpyn,0)
  vt_rpyn = 0.5*dim_sum_n(dim_avg_n(v_irtpyn*t_irtpyn,2),0)  ; t avg, i sum
  vt_pyn  = dim_avg_n(vt_rpyn,0)
  tt_rpyn = 0.5*dim_sum_n(dim_avg_n(t_irtpyn^2,2),0)  ; t avg, i sum
  tt_pyn  = dim_avg_n(tt_rpyn,0)
;  vvClim_pyn  = dim_avg_n( reshape(vv_rpyn(irClim,:,:),(/dimsizes(irClim),ds(3),ds(4)/2+1/)) ,0)
;  vvEN_pyn  = dim_avg_n( reshape(vv_rpyn(irEN,:,:),(/dimsizes(irEN),ds(3),ds(4)/2+1/)) ,0) - vvClim_pyn
  copy_VarCoords(u_py,vv_pyn(:,:,0))
  vv_pyn@long_name  = "v'v' at x hPa (m~S~2~N~ s~S~-2~N~)"
  uv_pyn@long_name  = "u'v' at x hPa (m~S~2~N~ s~S~-2~N~)"
  uu_pyn@long_name  = "u'u' at x hPa (m~S~2~N~ s~S~-2~N~)"
  vt_pyn@long_name  = "v'T' at x hPa (K m s~S~-1~N~)"
  tt_pyn@long_name  = "T'T' at x hPa (K~S~2~N~)"
;  copy_VarCoords(u_py,vvEN_pyn(:,0))
;  delete([/u,v,u_rpyiwn,v_rpyiwn ,v_irtpyn/])
  delete([/u,v,u_irtpyn,v_irtpyn,vv_rpyn,uv_rpyn,uu_rpyn,vt_rpyn,tt_rpyn/])
  wallClockElapseTime(wcStrt, "fft",0)

  uo := reshape(uo,(/1,ds_org(0),ds_org(1),ds_org(2),ds_org(3)/))
  vo := reshape(vo,(/1,ds_org(0),ds_org(1),ds_org(2),ds_org(3)/))
  ds_org := toint(dimsizes(uo))  ;yr,t,p,y,x
  ds_org := (/ds_org(0)*2,ds_org(1),ds_org(2),ds_org(3)/2,ds_org(4)/)
  uo2 = array_append_record(uo(:,:,:,0:ds(3)-1,:),uo(:,:,:,2*ds(3)-1:ds(3),:),0)
  vo2 = array_append_record(vo(:,:,:,0:ds(3)-1,:),-vo(:,:,:,2*ds(3)-1:ds(3),:),0)
  uo := uo2
  vo := vo2
  delete([/uo2,vo2/])

  vo_irtpyn = new((/2,ds_org(0),ds_org(1),ds_org(2),ds_org(3),ds_org(4)/2+1/),float)
  vo_irtpyn = 0.
  uo_irtpyn = vo_irtpyn
  vo_irtpyn(:,:,:,:,:,1:) = ezfftf(vo)
  uo_irtpyn(:,:,:,:,:,1:) = ezfftf(uo)
  vvb_yn = 0.5*dim_sum_n(dim_avg_n(dim_avg_n(vo_irtpyn,3)^2,(/1,2/)),0)  ; p avg, ^2, rt avg, i sum
  uub_yn = 0.5*dim_sum_n(dim_avg_n(dim_avg_n(uo_irtpyn,3)^2,(/1,2/)),0)  ; p avg, ^2, rt avg, i sum
  vvs_yn = 0.5*dim_sum_n(dim_avg_n(dim_avg_n(vo_irtpyn^2,3),(/1,2/)),0)  ; ^2, p avg, rt avg, i sum
  uus_yn = 0.5*dim_sum_n(dim_avg_n(dim_avg_n(uo_irtpyn^2,3),(/1,2/)),0)  ; ^2, p avg, rt avg, i sum
  vvb_yn@long_name  = "Barotropic v'v' (m~S~2~N~ s~S~-2~N~)"
  uub_yn@long_name  = "Barotropic u'u' (m~S~2~N~ s~S~-2~N~)"
  vvs_yn@long_name  = "Vert. avg. v'v' (m~S~2~N~ s~S~-2~N~)"
  uus_yn@long_name  = "Vert. avg. u'u' (m~S~2~N~ s~S~-2~N~)"
  delete([/uo,vo,uo_irtpyn,vo_irtpyn/])
  wallClockElapseTime(wcStrt, "fft",0)

; test fft unit
/;
  u_yiwn = new((/ds(3),2,ds(1),ds(4)/2+1/),float)
;  v_yiwn = new((/ds(3),2,ds(1),ds(4)/2+1/),float)
  u_yiwn = 0.
  u_yiwn(0,0,3,3) = 100.
  u_yiwn@M = ds(1)
  u_yiwn@N = ds(4)
  v_yiwn = u_yiwn
  uu = fft2db(u_yiwn(0,:,:,:))
  vv = fft2db(v_yiwn(0,:,:,:))

    uu_iwn = fft2df(uu)
    vv_iwn = fft2df(vv)
;print(uu_iwn)
    uv_wn = 0.5*dim_sum_n(uu_iwn*vv_iwn,0)
print(sum(uv_wn))  ; 5000
print(dim_avg_n(dim_avg_n(uu*vv,1),0))  ; 5000
;/

 ; save ds  lev_p lat u_py {uv,vv,uu,vt,tt}_{pywn,pyn} {vv,uu}{b,s}_yn
setfileoption("nc", "Format", "NetCDF4Classic")
system("rm "+fn_save)
fout = addfile(fn_save, "c")
;system("ln -sf "+get_script_prefix_name+"."+caseid+".data.nc "+get_script_prefix_name+".data.nc")

setfileoption(fout, "DefineMode", True)

fAtt           = True
fAtt@caseid    = caseid
fAtt@creation_date = systemfunc("date")
fAtt@dt        = dt
;fAtt@Ly        = Ly
fAtt@Lx        = Lx
fAtt@ds            = ds((/0,1,3,4/))
;fAtt@irEN          = irEN
;fAtt@irClim        = irClim
fileattdef(fout, fAtt)
; copy_VarAtts(f, fout)

 dimNames = (/"lev_p","lat","freq","kx"/)
 dimUnlim = (/False,False,False,False/)
; filedimdef(fout, dimNames, (/ds(3),ds(1),ds(4)/2+1/), dimUnlim)
; filedimdef(fout, dimNames, dimsizes(uv_rpywn), dimUnlim)
 filedimdef(fout, dimNames, dimsizes(uv_pywn), dimUnlim)
;  system("echo finish def dim")

 filevardef(fout, "lev_p", typeof(lev_p), "lev_p")
 filevardef(fout, "lat", typeof(lat), "lat")
 filevardef(fout, "u_py", typeof(u_py), dimNames(:1))
 filevardef(fout, "uv_pywn", "float", dimNames)
 filevardef(fout, "vv_pywn", "float", dimNames)
 filevardef(fout, "uu_pywn", "float", dimNames)
 filevardef(fout, "vt_pywn", "float", dimNames)
 filevardef(fout, "tt_pywn", "float", dimNames)
 filevardef(fout, "uv_pyn", "float", dimNames((/0,1,3/)))
 filevardef(fout, "vv_pyn", "float", dimNames((/0,1,3/)))
 filevardef(fout, "uu_pyn", "float", dimNames((/0,1,3/)))
 filevardef(fout, "vt_pyn", "float", dimNames((/0,1,3/)))
 filevardef(fout, "tt_pyn", "float", dimNames((/0,1,3/)))
 filevardef(fout, "vvb_yn", "float", dimNames((/1,3/)))
 filevardef(fout, "uub_yn", "float", dimNames((/1,3/)))
 filevardef(fout, "vvs_yn", "float", dimNames((/1,3/)))
 filevardef(fout, "uus_yn", "float", dimNames((/1,3/)))
;  system("echo finish def var1")
 filevarattdef(fout,"lev_p", lev_p)
 filevarattdef(fout,"lat", lat)
 filevarattdef(fout,"u_py", u_py)
 filevarattdef(fout,"uv_pywn", uv_pywn)
 filevarattdef(fout,"vv_pywn", vv_pywn)
 filevarattdef(fout,"uu_pywn", uu_pywn)
 filevarattdef(fout,"vt_pywn", vt_pywn)
 filevarattdef(fout,"tt_pywn", tt_pywn)
 filevarattdef(fout,"uv_pyn", uv_pyn)
 filevarattdef(fout,"vv_pyn", vv_pyn)
 filevarattdef(fout,"uu_pyn", uu_pyn)
 filevarattdef(fout,"vt_pyn", vt_pyn)
 filevarattdef(fout,"tt_pyn", tt_pyn)
 filevarattdef(fout,"vvb_yn", vvb_yn)
 filevarattdef(fout,"uub_yn", uub_yn)
 filevarattdef(fout,"vvs_yn", vvs_yn)
 filevarattdef(fout,"uus_yn", uus_yn)

setfileoption(fout, "DefineMode", False)
;  system("echo `date` finish def all")
; wallClockElapseTime(wcStrt, "define nc",0)

 fout->lev_p = (/lev_p/)
 fout->lat = (/lat/)
 fout->u_py = (/u_py/)
 fout->uv_pywn = (/uv_pywn/)
 fout->vv_pywn = (/vv_pywn/)
 fout->uu_pywn = (/uu_pywn/)
 fout->vt_pywn = (/vt_pywn/)
 fout->tt_pywn = (/tt_pywn/)
 fout->uv_pyn = (/uv_pyn/)
 fout->vv_pyn = (/vv_pyn/)
 fout->uu_pyn = (/uu_pyn/)
 fout->vt_pyn = (/vt_pyn/)
 fout->tt_pyn = (/tt_pyn/)
 fout->vvb_yn = (/vvb_yn/)
 fout->uub_yn = (/uub_yn/)
 fout->vvs_yn = (/vvs_yn/)
 fout->uus_yn = (/uus_yn/)
delete(fout) ;close the file
 wallClockElapseTime(wcStrt, "write nc",0)

end
exit

 ;:'<,'> s?filevardef(\([^,]*\), "\([^"]*\)",.*$?filevarattdef(\1,"\2", \2)? |noh
 ;:'<,'> s?filevardef(\([^,]*\), "\([^"]*\)",.*$?\1->\2 = (/\2/)? |noh
 ;:'<,'> s?filevardef(\([^,]*\), "\([^"]*\)",.*$?\2 = fdata->\2? |noh

; vim: set fdm=marker foldmarker=/;,;/:

