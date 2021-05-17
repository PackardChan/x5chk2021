; vertically interpolate FMS results to pressure level
; see PPZTT63_pchan.m for an older, more comprehensive version that includes geopotential height
; sbatch --account=kuang_lab -p huce_intel -J interp_HS94 -n 1 -t 1440 --mem=126820 -o "slurm" --wrap="ncl -Q nc.interp_HS94.ncl >& slurm-\${SLURM_JOB_NAME}-\${SLURM_JOBID}" --mail-type=END
; sbatch --account=kuang_lab -p huce_intel -J vint --array=1-44 -n 1 -t 1440 --mem=6000 -o "slurm.%a" --wrap='ncl -Q nc.fms.vint.ncl tid=${SLURM_ARRAY_TASK_ID}  >& slurm-${SLURM_JOB_NAME}-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}'
;TODO sbatch --account=kuang_lab -p huce_intel -J vint --array=1-10 -n 1 -t 1440 --mem=6000 -o "logfile/slurm.%a" --wrap='export CASENAME=a1ctrl; ncl -Q nc.fms.vint.ncl  >& logfile/slurm-${SLURM_JOB_NAME}-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}'
 ; 2 mins in huce_intel for 2*100
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

 ; Flexible Modeling System (FMS)
 ; MiMAROOT=/n/home05/pchan/model/f1-fms_TestsV1
 ; /n/kuanglfs/${USER}/jetshift/fms-output/${CASENAME}/fms-runscript.sbatch
 ; ${MiMAROOT}/src/atmos_spectral
 ; ${MiMAROOT}/src/shared/constants/constants.f90
 ; https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pjp/spectral_core.pdf, Eq. 33

begin
 wcStrt = systemfunc("date")
  tid = stringtoint(getenv("tid")) ;TODO
;  fileAd = "/n/holylfs/LABS/kuang_lab/pchan/fms_output/"  ; copied 20161118
;  text = "AOefbCV1HF"
;  caseid=text+"T63x1day"
;  fileAd = "fms-output/"
 ; \ls -1 fms-output/ > zc-cidarr_txt
;  cidArr = systemfunc("cat zc-cidarr_txt")
  caseid = getenv("CASENAME") ;TODO cidArr(tid-1)

;  nfiles = toint(systemfunc("ls -1dv "+fileAd+caseid+"/history/day*00h00/ |wc -l"))
;  Run = 100*ispan(11,nfiles,1)  ; discarding first 1000 days TODO
;  Run = 100*ispan(11,20,1)  ; discarding first 1000 days
;  Run = 100*ispan(tid*10+1,tid*10+10,1)  ; discarding first 1000 days
  Run = 1000+100*ispan((tid-1)+1,tid,1)  ; discarding first 1000 days

  name = sprinti("%05d",Run(0))  ; no list: obtain grid only
;  fn_t42   = addfile("../sks/int_z500_zg_day_MIROC-ESM-CHEM_historical_r1i1p1_19660101-20051231.nc","r")
;  fn_fms   = addfile(fileAd+text+"/history/day"+name+"h00/day"+name+"h00.1xday.nc","r")  ; no list
;  fn_fms   = addfile(fileAd+caseid+"/history/day"+name+"h00.1xday.nc","r")  ; no list
  fn_fms   = addfile("history/day"+name+"h00.1xday.nc","r")  ; no list
  ds_org = getfilevardimsizes(fn_fms,"ucomp")
  bk   = fn_fms->bk
  bk_f = fn_fms->pfull / fn_fms->phalf(dimsizes(bk)-1)  ; end
;  delete(fn_t42)
  delete(fn_fms)  ; no list

;  system("mkdir -p analysis")
 ; NCO -n option requires fixed width numbers. http://nco.sourceforge.net/nco.html#Specifying-Input-Files
  fn_save1 = "analysis/var_"+caseid+"_"+sprinti("%03d",tid)+".nc"  ;tostring(tid)
  fn_save2 = "analysis/zyvar_"+caseid+"_"+sprinti("%03d",tid)+".nc"
;  fn_save1 = "uv_"+caseid+".nc"
;  if (systemfunc("find "+fn_save1+" -newer "+get_script_name).eq.fn_save1) then
  if (.not.ismissing(systemfunc("find "+fn_save1+" -newer "+get_script_name))) then
    print("File exists: "+fn_save1)
    status_exit(0)
  end if

 ; vinth2p requires the lev_p to be expressed in mb [hPa]
;  lev_p = (/500., 250./)
;  lev_p = fspan(1000.,100.,10)
  lev_p = (/1000.,900.,800.,700.,600.,500.,400.,300.,240.,200.,100./)
  lev_p!0         = "lev_p"                  ; variable/dim name 
  lev_p&lev_p     =  lev_p                   ; create coordinate variable  
  lev_p@long_name = "pressure"               ; attach some attributes
  lev_p@units     = "hPa"
  lev_p@positive  = "up"
;  k850min = 31  ;find( bk_f>85000/max(ps(:)), 1)-1;  %33
;  k850max = 38  ;find( bk_f>85000/min(ps(:)), 1);    %38
  ds = (/dimsizes(Run)*100,dimsizes(lev_p),ds_org(2),ds_org(3)/)

;  U3 = new(ds*(/4,1,1,1/),float)
  u = new(ds,float)
  v = new(ds,float)
  t = new(ds,float)
  z = new(ds,float)
  z1 = new(ds_org+(/0,0,0,0/),float)
 wallClockElapseTime(wcStrt, "metadata",0)
 wcStrt = systemfunc("date")

do n = 0,dimsizes(Run)-1
  system("echo -ne "+str_get_dq+"\r"+Run(n)+str_get_dq)
  name = sprinti("%05d",Run(n))

  fn_fms = addfile("history/day"+name+"h00.1xday.nc","r")  ; no list

  ps    = fn_fms->ps   ;1xday(3::4  ,:,:)   ; [n] 3::4
  u1    = fn_fms->ucomp;1xday(3::4,:,:,:)
  v1    = fn_fms->vcomp;1xday(3::4,:,:,:)
  t1    = fn_fms->temp ;1xday(3::4,:,:,:)
;  W    = fn_fms->omega(:,k850min:k850max,:,:)

  delete(fn_fms)  ; no list
; wallClockElapseTime(wcStrt, "read nc",0)
; wcStrt = systemfunc("date")

 ; cz2ccm, hydro
;  wrk = conform(t_spy,pTrop_sy/ps_sy,(/0,2/))
;  z1 = dim_cumsum_n(log(dim_max_n((/conform(t1,bk(1:),1),wrk/),0)/dim_max_n((/conform(t1,bk(:ds_org(0)-1),1),wrk/),0)) *287.04/9.80*t1,1)
  z1 = 0.
  z1(:,1:,:,:) = dim_cumsum_n(conform(t1(:,1:,:,:),log(bk(2::-1)/bk(1:ds_org(1)-1:-1)),1)*287.04/9.80*t1(:,1::-1,:,:),0,1) ; avoid div 0
  z (n*100 : (n+1)*100-1,:,:,:) = vinth2p(z1(:,::-1,:,:),bk_f*0,bk(1:),lev_p,ps,2,1000.,1,True) ; 2 means log interpolation
;  U850 = rm_single_dims(vinth2p(u1,bk_f(k850min:k850max)*0,bk_f(k850min:k850max),850,ps,1,1000,1,False))
  u (n*100 : (n+1)*100-1,:,:,:) = vinth2p(u1,bk_f*0,bk_f,lev_p,ps,2,1000.,1,True) ; 2 means log interpolation
  v (n*100 : (n+1)*100-1,:,:,:) = vinth2p(v1,bk_f*0,bk_f,lev_p,ps,2,1000.,1,True) ; 2 means log interpolation
  t (n*100 : (n+1)*100-1,:,:,:) = vinth2p(t1,bk_f*0,bk_f,lev_p,ps,2,1000.,1,True) ; 2 means log interpolation
;  T2 = vinth2p(T,bk_f*0,bk_f,lev_p,ps,2,1000.,1,True) ; 2 means log interpolation
; wallClockElapseTime(wcStrt, "v intp",0)
; wcStrt = systemfunc("date")

;  g2gshv(U850,V850,U3,V3,0)
;  g2fshv_Wrap(U2,U2*0,U3,V3)
;  T3 = g2fsh_Wrap(T2,ds(2:3))
; wallClockElapseTime(wcStrt, "g2f uv",0)
; wcStrt = systemfunc("date")

;  u (n*100 : (n+1)*100-1,:,:) = dim_avg_n_Wrap(U2,3)
;  v (n*100 : (n+1)*100-1,:,:) = dim_avg_n_Wrap(V2,3)
;  V850Daily (n*100 : (n+1)*100-1,:,:) = V3
; wallClockElapseTime(wcStrt, "x mean",0)
; wcStrt = systemfunc("date")
end do

;  delete(fn_fms)  ; list
 wcStrt = systemfunc("date")
  lat = u1&lat
  lat_wgt = latGauWgt(dimsizes(lat),lat!0,"gaussian weights","")
  lon = u1&lon
;delete([/ps,u1,v1/])
delete([/ps,u1,v1,t1,z1/])


;===============================================
;setfileoption("nc", "Format", "64BitOffset")
setfileoption("nc", "Format", "NetCDF4Classic")
system("rm "+fn_save1)
fout = addfile(fn_save1, "c")

setfileoption(fout, "DefineMode", True)

fAtt           = True
fAtt@caseid    = caseid
fAtt@creation_date = systemfunc("date")
fileattdef(fout, fAtt)
; copy_VarAtts(f, fout)

 dimNames = (/"time","lev_p","lat","lon"/)
; ds = (/dimsizes(Run)*100,64,128/)
 dimUnlim = (/True,False,False,False/)
 filedimdef(fout, dimNames, ds, dimUnlim)
  system("echo finish def dim")

;filedimdef(fout,"time",-1,True) 
 filevardef(fout, "lev_p", typeof(lev_p), "lev_p")
 filevardef(fout, "lat", typeof(lat), "lat")
 filevardef(fout, "lat_wgt", typeof(lat_wgt), "lat")
 filevardef(fout, "lon", typeof(lon), "lon")
 filevardef(fout, "u", "float", dimNames)
 filevardef(fout, "v", "float", dimNames)
 filevardef(fout, "t", "float", dimNames)
 filevardef(fout, "z", "float", dimNames)
  system("echo finish def var1")
 filevarattdef(fout,"lev_p", lev_p)
 filevarattdef(fout,"lat", lat)
 filevarattdef(fout,"lat_wgt", lat_wgt)
 filevarattdef(fout,"lon", lon)
 filevarattdef(fout,"u", u)
 filevarattdef(fout,"v", v)
 filevarattdef(fout,"t", t)
 filevarattdef(fout,"z", z)

setfileoption(fout, "DefineMode", False)
  system("echo `date` finish def all")
 wallClockElapseTime(wcStrt, "define nc",0)
 wcStrt = systemfunc("date")

 fout->u = (/u/)
 fout->v = (/v/)
 fout->t = (/t/)
 fout->z = (/z/)
 fout->lev_p = (/lev_p/)
 fout->lat = (/lat/)
 fout->lat_wgt = (/lat_wgt/)
 fout->lon = (/lon/)
;fout->latb=latb
;fout->lonb=lonb
delete(fout) ;close the file
delete([/u,v,t,z/])
system("ncwa -O -a time,lon "+fn_save1+" "+fn_save2)
 ; fms-runscript.sbatch: nco commands to follow
 ; CASENAME=x2srcmod
 ; ncra `\ls -1v analysis/var_${CASENAME}_*.nc` tmp.ncra.${CASENAME}.nc
 ; ncwa -a lon tmp.ncra.${CASENAME}.nc zyvar_${CASENAME}.nc
 ; ncdiff zyvar_${CASENAME}.nc zyvar_a1ctrl.nc diff_zy_${CASENAME}.nc

/;

filename = "day45100h00.1xday.nc"
F = addfiles(readin_folder+filename, "r")


do ii=1,1 ;the first time snapshot. here, for the purpose of illustration, we only calculate the first snapshot of this netcdf data. you need to calculate ii=1,400

tstart=ii-1
uwind0 = F[:]->ucomp(tstart,:,:,:)
pfull = F[:]->pfull
ps0 = F[:]->ps(tstart,:,:)
latb=F[:]->latb
lonb=F[:]->lonb
printVarSummary(uwind0)

;------------------------PART I: regrid--------------------------------
;now do the regridding
;horizontal regridding
uwind = g2fsh_Wrap (uwind0, ds(2:3))
ps = g2fsh_Wrap (ps0, ds(2:3))

lev_p!0         = "lev_p"                  ; variable/dim name 
lev_p&lev_p     =  lev_p                   ; create coordinate variable  
lev_p@long_name = "pressure"               ; attach some attributes
lev_p@units     = "hPa"
lev_p@positive  = "up"

;some other variables
hybm            = pfull/1.e3
hyam            = hybm            ; create a bogus hyam 1D array
hyam            = 0.0             ; set all values to zero
P0mb            = 1000.
kxtrp = True
tbot=ps*0.0
PHIS=ps*0.0
;vertical regridding
uwind2  = vinth2p(uwind,hyam,hybm,lev_p,ps,2,P0mb, 1, kxtrp) ; 2 means log interpolation

;------------------------PART X: calculate the zonal mean  --------------------------------
tbar = dim_avg_n_Wrap(uwind2,2)

;------------------------PART II: write out-----------------------------------
;now, write out the variables
filename_out="output_interp_HS94_"+sprinti("%0.3i", ii)+".nc"

system("rm -f "+writeout_folder+filename_out)
fout = addfile(writeout_folder+filename_out,"c")
fout@title = "simple sample"
fout@creation_date = systemfunc("date")
filedimdef(fout,"time",-1,True) 
fout->tbar = tbar
fout->latb=latb
fout->lonb=lonb
delete(fout) ;close the file
delete(uwind)
delete(uwind2)

end do
;/

end

; :'<,'>s?filevardef(\([^,]*\), "\([^"]*\)",.*$?\1->\2 = (/\2/)? | noh
exit

; vim: set fdm=marker foldmarker=/;,;/:

