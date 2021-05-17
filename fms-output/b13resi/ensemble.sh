#!/bin/sh

 # TODO Usage:
 # Prepare in current directory:
 #   ForcingM2EOF.txt
 #   fms-runscript.sbatch
 #   srcmods/  (opt)
 # Then, ./ensemble.sh

CASEprefix=$(basename $PWD)  # TODO
MiMAROOT=/n/home05/pchan/model/f1-fms_TestsV1   # original source code
oldpwd=$PWD

EnsSize=20
for ii in `seq 1 $EnsSize`; do
  cd $oldpwd
  CASENAME=${CASEprefix}ens$ii
  if [ ! -f ../${CASENAME}/history/day26000h00.1xday.nc ] && [ -z "$(\sacct -S `date +%H:%M` --name=$CASENAME -XnP)" ]; then
#  if [ $(\ls ../${CASENAME}/history/*.nc |wc -l) -lt 260 ] && [ -z "$(squeue |awk -v cn=$CASENAME '$3==cn {print}')" ]; then
    mkdir -p ../$CASENAME
    cp -a ForcingM2EOF.txt fms-runscript.sbatch ../$CASENAME/
    [ -d srcmods ] && cp -a srcmods ../$CASENAME/
    mkdir -p ../$CASENAME/srcmods
    cp -a $MiMAROOT/src/atmos_spectral/init/spectral_initialize_fields.f90 ../$CASENAME/srcmods/
    sed -i "/real :: initial_perturbation/s/=.*$/= ${ii}.e-7/" ../$CASENAME/srcmods/spectral_initialize_fields.f90
    cd ../$CASENAME
#    cp -a commands_reload commands_reload-kuanglfs
#    sed -i "s?kuanglfs?holylfs/LABS/kuang_lab?" commands_reload
    sbatch -J $CASENAME -p huce_intel -t 1690 fms-runscript.sbatch
#    rm $oldpwd/../../finish-${CASENAME}
    rm finish
    sleep 1s
  fi
done

cd $oldpwd/../..
echo $CASEprefix $EnsSize b1ctrl >> crontab-running-ensemble
exit

 # Jobs to run:
 #   fms (fms-runscript.sbatch)
 #   1vint, per nc file (fms-runscript.sbatch)
 #   2nces, per ensemble member (fms-runscript.sbatch)
 #   ncdiff, per ensemble (crontab-ensemble-fms-monitor.sh)
 # crontab is used to run ncdiff because 2nces is not yet submitted.

# JOB=`sbatch -J $CASENAME -p huce_intel -t 720 fms-runscript.sbatch | egrep -o -e "\b[0-9]+$"`
# JOBlist="$JOBlist:$JOB"
# echo $JOBlist
#done
#sbatch -p huce_amd -J "summarize" -n 1 -t 10 --mem=2000 "--dependency=afterok$JOBlist" -o /dev/null --wrap='sleep 30s; exit'
#sbatch --account=kuang_lab -p huce_intel,test -J 2nces_${CASENAME} -n 1 -t 60 --mem=$(( $iMax * 100 +5000 )) --dependency=afterok:${JOB} -o slurm --mail-type=END <<'EOF'
#!/bin/sh
#set -x
#EOF

( for ii in `seq 1 20`; do \ls fms-output/b10unifZeroUens$ii/history/*nc -1 |wc -l ; done )|uniq -c

 # postprocess (1vint & 2nces) for each ensemble member -> fms-runscript.sbatch

 # check each ensemble member
CASEprefix=b5unif
#rm zc-cidarr_txt
echo ${CASEprefix} > zc-cidarr_txt
for ii in `seq 10`; do
  CASENAME=${CASEprefix}ens$ii
  ncdiff -O zyvar_${CASENAME}.nc zyvar_${CASEprefix}.nc diff_zy_${CASENAME}.nc
  echo ${CASENAME} >> zc-cidarr_txt
done
ncl -Q yp.diff.ncl

 # clean fms
CASEprefix=$(basename $PWD)
oldpwd=$PWD
for ii in `seq 1 20`; do
  CASENAME=${CASEprefix}ens$ii
  cd ../$CASENAME
#  rm $oldpwd/../../finish-${CASENAME}
#  rm -r history/ logfile/ restart/ obj/ commands_reload commands_reload_prev fms.x
  \ls -1vr restart/*.cpio |sed '1d' |xargs rm
done

 # clean manyfiles
CASEprefix=$(basename $PWD)
#oldpwd=$PWD
#for ii in `seq 1 20`; do
#  CASENAME=${CASEprefix}ens$ii
#  rm manyfiles/var_${CASENAME}_[0-9][0-9][0-9].nc manyfiles/zyvar_${CASENAME}_[0-9][0-9][0-9].nc manyfiles/zyorg_${CASENAME}_[0-9][0-9][0-9].nc
#  rm manyfiles/slurm.vint_${CASENAME}.[0-9] manyfiles/slurm.vint_${CASENAME}.[0-9][0-9] manyfiles/slurm.vint_${CASENAME}.[0-9][0-9][0-9]
#  rm cy.cospectrum.${CASENAME}.data.nc zyvar_${CASENAME}.nc zyorg_${CASENAME}.nc finish-${CASENAME}
#  rm slurm-ncra_${CASENAME}
#done
#rm manyfiles/*nc.pid*.tmp

 # clean 1vint 2nces
find ./ -path '*/analysis/nc.1cospectrum_[0-9][0-9][0-9].nc' -delete
CASEprefix=$(basename $PWD)
for ii in `seq 1 20`; do
  CASENAME=${CASEprefix}ens$ii
#  rm -r ../${CASENAME}/analysis
  rm ../${CASENAME}/analysis/nc.1cospectrum_[0-9][0-9][0-9].nc
#  rm ../${CASENAME}/logfile/slurm.vint_${CASENAME}.[0-9]*

#  rm ../${CASENAME}/cy.cospectrum.${CASENAME}.data.nc ../${CASENAME}/zyvar_${CASENAME}.nc ../${CASENAME}/zyorg_${CASENAME}.nc ../${CASENAME}/finish
#  rm ../${CASENAME}/slurm-ncra_${CASENAME}
#  rm ../${CASENAME}/slurm-2nces_${CASENAME}
done
#rm manyfiles/*nc.pid*.tmp

 # clean ncdiff
CASEprefix=$(basename $PWD)
rm ../../ensemble-wise/zyvar_${CASEprefix}.nc ../../ensemble-wise/zyorg_${CASEprefix}.nc ../../ensemble-wise/diff_zy_${CASEprefix}.nc

