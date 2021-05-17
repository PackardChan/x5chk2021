#!/bin/bash

#SBATCH --account=kuang_lab
#SBATCH -p huce_intel  #,test
#SBATCH -J fms
#SBATCH -N 1
#SBATCH -n 32  # num_lat_rows(96)/num_pes=int
##SBATCH --ntasks-per-node=8
#SBATCH --hint=compute_bound
##SBATCH --contiguous  #TODO
##SBATCH -w holy2a18106   #.rc.fas.harvard.edu
#SBATCH -x /n/home05/pchan/sw/crontab/node-exclude
#SBATCH -t 0-8:0:0  #999999  # minutes
#SBATCH --mem-per-cpu=1000  #4025
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH -o slurm.out-%j
##SBATCH -e eo  #-%j
#SBATCH --mail-type=FAIL,REQUEUE  #TODO
##SBATCH --mail-user=name@example.com  # default $USER@holy-slurm01.rc.fas.harvard.edu

 # ody: /n/home11/cwalker/FMS/atm_dycores/exp/spectral/fms_runscript
 # che: /glade/u/home/wykang/build_mima_combine.sh
 # ody: ~pedram/fms_lrf/scripts/NewRuns/P900L90p15/P900L90p15
 # ody: ~pedram/fms_TestsV1/scripts/AObaroefbV1HF/AObaroefbV1HF
 # ody: ~pchan/script/fms/fms-runscript.sbatch

 # TODO Usage:
 # cp -ai ~pedram/fms_TestsV1 $MiMAROOT
 # use baroclinicity.m to:
 #   mkdir fms-output/${CASENAME}
 #   generate ForcingM2EOF.txt
 # copy this script to fms-output/${CASENAME}/
 # make changes to this script, e.g. runs_per_script, num_script_runs
 # (opt) srcmods/
 # sbatch fms-runscript.sbatch

VERSION=1.0.0
DATE="May 24, 2018"
#USAGE="Usage: $0 CASENAME"
#OPTIONS="casename [default:test], ncores [default:32]"
#----------------------------------------------------------------
## get in arguments ##
#----------------------------------------------------------------
#if [ $# == 0 ] ; then
#    echo $USAGE
#    echo $OPTIONS
#    exit 1;
#fi
CASENAME=$(basename $PWD)  # TODO
MiMAROOT=/n/home05/pchan/model/f1-fms_TestsV1	# original source code
#MiMAHOME=${HOME}/fms_caseroot    # srcmod, input-mod, caseroot
scriptdir=/n/holylfs/LABS/kuang_lab/${USER}/jetshift/fms-output/${CASENAME}  # location of this sbatch script
#execdir=${MiMAHOME}/${CASENAME}				# exeroot: compile and generate executable
execdir=/n/holylfs/LABS/kuang_lab/${USER}/jetshift/fms-output/${CASENAME}	# exeroot: compile and generate executable
#workdir=/scratch/${USER}/mima_output/${CASENAME}	# rundir, lfs: running model and output
workdir=/scratch/${USER}/${SLURM_JOBID}-${SLURM_RESTART_COUNT}	# rundir, lfs: running model and output
#workdir=/n/holylfs/LABS/kuang_lab/${USER}/temp/${SLURM_JOBID}	# rundir, lfs: running model and output
#scratchdir=/scratch/${CASENAME}
output_dir=/n/holylfs/LABS/kuang_lab/${USER}/jetshift/fms-output/${CASENAME}	# short-term archive, not node local
runs_per_script=999                           # TODO number of runs
num_script_runs=1
days=100
daystop=26000
daynow=0

irun=1                                             # loop counter
ireload=1                                          # resubmit counter
init_cpio=""

# if exists, load reload file ( irun ireload init_cpio )

reload_file=$output_dir/commands_reload

if [ -f $reload_file ] ; then
   echo "Found commands_reload"
   source $reload_file
fi

echo "ireload =" $ireload

#export PLATFORM=impi
#export NCORES=${2:-"32"}							# how many cores to use
#NCORES=${SLURM_NTASKS}							# how many cores to use
#export NCORE_NODE=36								# machine attribution
#export NNODE=$(((NCORES - 1) / NCORE_NODE + 1))		# should be integer
#export NSUBMIT=1									# how many times to rerun

#----------------------------------------------------------------
## set PATH and environment variables##
#----------------------------------------------------------------
## new library ....................................
#export compiler_type=intel
#export compiler_name=intel/17.0.1
#export compiler_fmulti=mpiifort
#export compiler_cmulti=mpiicc
#export compiler_fsingle=ifort
compiler_csingle=icc
#export mpi_type=mpt
#export mpi_name=mpt/2.15f # support both MPI and OpenMP 
#export omp_cmd=omplace
#export mpi_cmd='mpiexec -n $npes'
#export netcdf_name=netcdf/4.4.1.1

module purge
module load intel/17.0.4-fasrc01 impi/2017.2.174-fasrc01 netcdf/4.1.3-fasrc02
#module load nco/4.7.4-fasrc01
#module load ncl_ncarg/6.4.0-fasrc01

#module purge
#module load ${compiler_name}
#module load ${mpi_name}
#module load ${netcdf_name}
#module load mkl

#NETCDF_HOME=`which ncdump`;NETCDF_HOME=${NETCDF_HOME##* };NETCDF_HOME=${NETCDF_HOME%/*/*}
#echo NETCDF_HOME = $NETCDF_HOME
#export NETCDF_INC=${NETCDF_HOME}/include
#export NETCDF_LIB=${NETCDF_HOME}/lib

#MPI_HOME=`which mpirun`;MPI_HOME=${MPI_HOME##* };MPI_HOME=${MPI_HOME%/*/*}
#echo MPI_HOME = $MPI_HOME
#export MPI_INC=${MPI_HOME}/include
#export MPI_LIB=${MPI_HOME}/lib

#COMPILER_HOME=`which ${compiler_fsingle}`;COMPILER_HOME=${COMPILER_HOME##* };COMPILER_HOME=${COMPILER_HOME%/*/*/*}
#echo COMPILER_HOME = $COMPILER_HOME
#export COMPILER_INC=${COMPILER_HOME}/include
#export COMPILER_LIB=${COMPILER_HOME}/lib

#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${NETCDF_LIB}
#echo LD_LIBRARY_PATH = $LD_LIBRARY_PATH

#----------------------------------------------------------------
## customer code/namelist ##
#----------------------------------------------------------------
# define environment variable used in compiling
platform=impi
npes=${SLURM_NTASKS}
template=${MiMAROOT}/bin/mkmf.template.$platform	# path_name list template for compiling 
mkmf=${MiMAROOT}/bin/mkmf							# path to executable mkmf
sourcedir=${MiMAROOT}/src							# path to directory containing model source code
mppnccombine=${MiMAROOT}/bin/mppnccombine.$platform    # path to executable mppnccombine
time_stamp=${MiMAROOT}/bin/time_stamp.csh            # generates string date for file name labels

#pathnames_tpl=${MiMAROOT}/exp/path_names			# compile file list template
pathnames=$execdir/path_names						# compile file list
pathnames_tpl=${MiMAROOT}/input/cw_spectral_pathnames       # path to file containing list of source paths
namelist=${MiMAROOT}/input/ph_spectral_namelistT63L40hv2    # path to namelist file
diagtable=${MiMAROOT}/input/d2-diag_table-1xday             # path to diagnositics table
fieldtable=${MiMAROOT}/input/cw_spectral_field_table        # path to field table (specifies tracers)
#srcmoddir=${MiMAHOME}/srcmods						# where hold the modified source code
#inputdir=${MiMAHOME}/input							# where hold the modified input files

mkdir -p $execdir/
 #cp /glade/u/home/pchan/fms_caseroot/t8port//INPUT/* ./INPUT/
 #ln -s $MiMAROOT/input/INPUT .
#cp -f ${MiMAROOT}/input/* $execdir/ 2>/dev/null 
#cp -f ${inputdir}/* $execdir/ 2>/dev/null # overwrite, if any, input.nml, diag_table and field_table
#cp $pathnames_tpl $execdir/ 2>/dev/null
#mkdir -p $execdir/srcmods

# srcmod
#if [ "$(ls -A $execdir/srcmods)" ]; then
#	cp -rf $srcmoddir $execdir/srcmods 2>/dev/null # modified code cp and add to the head of path_names list
	find $execdir/srcmods -maxdepth 1 -iname "*.f90" -o -iname "*.inc" -o -iname "*.c" -o -iname "*.h" > ${pathnames}
	echo "Modifying the following SourceCode"
	cat ${pathnames}
#fi
cat $pathnames_tpl >> ${pathnames}

# prepare workdir
if [ -e $workdir ]; then
  echo "ERROR: Existing workdir already exist, running will overwrite workdir. Please Move or remove $workdir and try again."
  exit 1
fi
mkdir -p $workdir
mkdir -p $workdir/INPUT
mkdir -p $workdir/RESTART
#cp ~pedram/fms_TestsV1/scripts/AObaroefbV1HF/ForcingM2EOF.txt $workdir/
cp ForcingM2EOF.txt $workdir/

mkdir -p $output_dir
mkdir -p $output_dir/history
mkdir -p $output_dir/logfile
mkdir -p $output_dir/restart
cp -aL $diagtable $output_dir/
cp -aL $namelist $output_dir/

#--------------------------------------------------------------------------------------------------------
# set run length and time step, get input data and executable
#cd $workdir
#cat > input.nml <<EOF
# &main_nml
#     days   = 1,
#     dt_atmos = 1800 /
#EOF
#cat $namelist >> input.nml
#cp $diagtable diag_table
#cp $fieldtable field_table
#cp $execdir/fms.x fms.x

#cd $execdir
#cp input.nml field_table diag_table $workdir/

#----------------------------------------------------------------
## compile ##
#----------------------------------------------------------------

# C compile the nc file combiner
if [ ! -f $mppnccombine ]; then
  ${compiler_csingle} -O -o ${mppnccombine} -I${NETCDF_INCLUDE} -L${NETCDF_LIB} -lnetcdf -I${MPI_INCLUDE} -L${MPI_LIB} -lmpi ${MiMAROOT}/postprocessing/mppnccombine.c
fi

# compile the model code and create executable
cd $execdir
cppDefs="-Duse_libMPI -Duse_netCDF"
$mkmf -p fms.x -t $template -c "$cppDefs" -a $sourcedir $pathnames $sourcedir/shared/mpp/include $sourcedir/shared/include
#$mkmf -p fms.x -t $template -c "$cppDefs" -a $sourcedir $pathnames /usr/local/include $NETCDF_INCLUDE $sourcedir/shared/mpp/include $sourcedir/shared/include
time ( make -j${SLURM_NTASKS} -f Makefile >& $execdir/bldlog )
mkdir obj
mv *.o *.mod Makefile bldlog .cppdefs path_names  obj/

# set initial conditions and move to executable directory

cd $workdir
if [ -n "$init_cpio" ];  then
  cd INPUT
  echo $init_cpio
  echo $(basename $init_cpio) 
  cp $init_cpio $(basename $init_cpio)
  cpio -iv  < $(basename $init_cpio)
  rm -f $(basename $init_cpio)
fi

#----------------------------------------------------------------
## generate run script ##
#----------------------------------------------------------------
#mpirun -np $npes fms.x
#--------------------------------------------------------------------------------------------------------
# combine netcdf files
#if [ $npes -gt 1 ]; then
#  for ncfile in *.nc.0000; do
#    $mppnccombine $ncfile:r
#    if ($status == 0) rm -f $ncfile:r.????
#  end
#endif

cd $workdir    
#--------------------------------------------------------------------------------------------------------

#  --- begin loop over $runs_per_script ---

  while [ $daynow -lt $daystop -a $irun -le $runs_per_script ]; do

#--------------------------------------------------------------------------------------------------------

# set run length and time step, get input data and executable

if [ $irun -eq 1 -a $ireload -eq 1 ];  then
  cat > input.nml <<EOF
 &main_nml
     current_time = 0,
     days   = $days,
     dt_atmos = 900,
     override = .true. /
EOF

else  # irun -ne 1 .or. ireload -ne 1
  cat > input.nml <<EOF
 &main_nml
     days   = $days,
     dt_atmos = 900 /
EOF

fi

cat >> input.nml <<EOF
 &atmosphere_nml
     lconst_forcing  = .true.,
     moist_tropics   = .false.,
     turb            = .false.,
     ldry_convection = .false.,
     roughness       = 0.05 /

 &const_forcing_nml
     p_k    = 30000.0,     ! pressure of forcing maximum
     lat_j  = 60.0,         ! latitude of forcing maximum
     u_amp  = 1.0,     ! amplitude of u forcing (m/s/s) -- default = 0.05m/s/day
     t_amp  = 1.0, ! amplitude of t forcing (k/s) -- default = 0.05K/day
     forcing_t = .true.,    ! should t with constant forcing applied
     forcing_u = .true. / ! should u with constant forcing applied
        
 &spectral_init_cond_nml
  initial_temperature      = 280.0 /

 &hs_forcing_nml
   delh = 60.0,    
   ka   = -40.0,      
   t_strat = 200.0 /
   
EOF
 # still in $workdir
cat $namelist >> input.nml
cp $diagtable diag_table
cp $fieldtable field_table
cp $execdir/fms.x fms.x

#--------------------------------------------------------------------------------------------------------
#time srun -n $SLURM_NTASKS --mpi=pmi2 --cpu_bind=cores --hint=compute_bound ./fms.x
time srun -n $npes --mpi=pmi2 --cpu_bind=cores --hint=compute_bound ./fms.x
#time mpirun -np $npes --bind-to core ./fms.x
exitstatus=$?
if [ $exitstatus -ne 0 ]; then
  set -x
  mail -s "srun ended" ${USER} << EOF
Job_id=${SLURM_JOBID}
caseid=${CASENAME}
ExitCode $exitstatus

tail slurm.out-${SLURM_JOBID}
`tail ${SLURM_SUBMIT_DIR}/slurm.out-${SLURM_JOBID}`
EOF
  wait
  sleep 60s
  exit $exitstatus
fi
#--------------------------------------------------------------------------------------------------------

#   --- generate date for file names ---
 # read in time_stamp.out from exe
date_name=$($time_stamp -eh)
if [ ! -n "$date_name" ]; then
 date_name=tmp$(date '+%j%H%M%S')
fi
#echo $date_name
if [ -f time_stamp.out ]; then
 rm -f time_stamp.out
fi

#--------------------------

#   --- move output files to their own directories (don't combine) ---
echo $date_name " irun =" $irun
daynow=`sed 's/day0*//;s/h00//' <<< "$date_name"`

#mkdir $output_dir/history/$date_name
#mkdir $scratchdir

#echo `date +%FT%T` "start moving"

shopt -s nullglob
#for ncfile in *.nc *.nc.????; do
#  mv $ncfile $scratchdir/$date_name.$ncfile
#done

echo `date +%FT%T` "start combining"

# combine 1xday data
#time $mppnccombine -r $scratchdir/$date_name.1xday.nc
time $mppnccombine -r 1xday.nc

echo `date +%FT%T` "start moving"

#mv $scratchdir/$date_name.1xday.nc $output_dir/history/
mv 1xday.nc $output_dir/history/$date_name.1xday.nc

echo `date +%FT%T` "done moving"

#rm -r $scratchdir

#echo "Done with Deleting"
#date

# take the zonal average
#ncwa -a lon $output_dir/history/$date_name.1xday.nc $output_dir/history/$date_name.za.1xday.nc
#rm -f $output_dir/history/$date_name.1xday.nc


#remove data from first 100 days:
#if [ $ireload -eq 1 ]; then
#  rm -rf $output_dir/history/
#fi


#   --- save ascii output files to local disk ---
 # logfile.out
for out in  *.out; do
  mv $out $output_dir/logfile/$date_name.$out
done


#   --- move restart files to output directory --- 

cd RESTART  # $workdir/RESTART
resfiles=(*.res*)
if [ ${#resfiles[@]} -gt 0 ]; then
#     --- desired filename for cpio of output restart files ---	
  restart_cpio=$output_dir/restart/$date_name.cpio
   [ ! -d $(dirname $restart_cpio) ] && mkdir -p $(dirname $restart_cpio)
#     --- also save namelist and diag_table ---
  cp $workdir/{*.nml,diag_table} .
  files=(${resfiles[*]} input.nml diag_table)
  /bin/ls ${files[*]} | cpio -ov > $(basename $restart_cpio)
  mv $(basename $restart_cpio) $restart_cpio
#     --- up restart for next run ---
  if [ $irun -lt $runs_per_script ]; then
     mv -f *.res*  ../INPUT
  fi
fi

cd $workdir


#--------------------------------------------------------------------------------------------------------

#   --- write new reload information ---
let irun++

[ -f $reload_file ] && mv -f $reload_file $reload_file"_prev"

if [ $daynow -lt $daystop -a $irun -le $runs_per_script ];  then
  echo "irun=$irun"          >  $reload_file
else
  let ireload++
  echo "irun=1"              >  $reload_file
fi

echo   "ireload=$ireload"       >> $reload_file
echo   "init_cpio=$restart_cpio"  >> $reload_file


# post processing

#cd $scriptdir

#if($irun > 2) then # submit analysis and copy data to atmos

#  echo "delh = $delh"    > ./post_processing_info_7
#  echo "gamma = $gamma"    >> ./post_processing_info_7
#  echo "date_name = $date_name" >> ./post_processing_info_7

#  llsubmit post_processing_7

#endif

cd $workdir

done ################# loop over $runs_per_script ended ###################

#cd $scriptdir
cd
rm -rf $workdir


# resubmit 
if [ $daynow -lt $daystop -a $ireload -le $num_script_runs ]; then 
  echo `date +%FT%T` "re-submitting"
  cd $scriptdir
  sbatch $0
  exit 0
fi

#[ $ireload -gt $num_script_runs ] && exit 0
#if [ $ireload -gt $num_script_runs ]; then 
#  cd /n/home04/pedram/fms_lrf/scripts/NewRuns/P500L90n15
#  sbatch P500L90n15  

#  exit 0
#fi

: <<'EOF'
mod18
for jj in 1 12 14 15 16 2 13 9; do
  CASEprefix=`basename ~/fms-backup/b${jj}[a-z]*`
#  CASEprefix=$(basename $PWD)
  echo $CASEprefix
  for ii in `seq 1 20`; do
    CASENAME=${CASEprefix}ens$ii
    cd ~/zho/jetshift/fms-output/$CASENAME
    rm finish
    sleep 1s
    export CASENAME
    export iPerTask=5
    export iMax=250  #TODO

    JOB=
  done
done
EOF

#cd $output_dir/../..
cd $output_dir
module load nco/4.7.4-fasrc01
set -x
export CASENAME
export iPerTask=5
export iMax=$(( (`sed 's/day0*//;s/h00//' <<< "$date_name"` - 1000)/$days ))
[ $iMax -le 0 ] && exit 0
cp -a ../../nc.fms.vint.ncl ./
cp -a ../../nc.1cospectrum.ncl ./
mkdir -p analysis
JOB=$( sbatch --account=kuang_lab -p huce_intel -J 1vint_${CASENAME} --array=1-$iMax:$iPerTask -n 1 -x ~/sw/crontab/node-exclude -t 60 --mem=6000 -o "logfile/slurm.%x.%a" --mail-type=FAIL,TIME_LIMIT <<'EOF' |tee /dev/tty |egrep -o -e "\b[0-9]+$"
#!/bin/sh
set -x
export tid=$SLURM_ARRAY_TASK_ID
while [ $tid -lt $(( $SLURM_ARRAY_TASK_ID + $iPerTask )) -a $tid -le $iMax ]; do
#  ncwa -a lon history/day$(printf "%05d" $(($tid *100+1000)))h00.1xday.nc analysis/tzyorg_${CASENAME}_$(printf "%03d" $tid).nc
  ncwa -a time,lon history/day$(printf "%05d" $(($tid *100+1000)))h00.1xday.nc analysis/zyorg_${CASENAME}_$(printf "%03d" $tid).nc
  ncl -Q nc.fms.vint.ncl
  ncl -Q nc.1cospectrum.ncl
  let tid++
done
EOF
)
#JOB=`sbatch --account=kuang_lab -p huce_intel -J vint_${CASENAME} --array=1-$iMax:$iPerTask -n 1 -x ~/sw/crontab/node-exclude -t 60 --mem=6000 -o "logfile/slurm.%x.%a" --mail-type=FAIL,TIME_LIMIT --wrap='export tid=$SLURM_ARRAY_TASK_ID; while [ $tid -lt $(( $SLURM_ARRAY_TASK_ID + $iPerTask )) -a $tid -le $iMax ]; do  ncwa -a time,lon history/day$(printf "%05d" $(($tid *100+1000)))h00.1xday.nc analysis/zyorg_${CASENAME}_$(printf "%03d" $tid).nc; ncl -Q nc.fms.vint.ncl; let tid++; done' | egrep -o -e "\b[0-9]+$"`

 # if job fails
export JOB
sbatch --account=kuang_lab -p huce_intel,test -J notok_${CASENAME} -n 1 -x ~/sw/crontab/node-exclude -t 60 --mem=1000 --dependency=afternotok:${JOB} -o "slurm-%x" --mail-type=BEGIN <<'EOF'
#!/bin/sh
set -x
#\sacct -j $JOB -s F,CA,OOM,TO -XnP -o NodeList |uniq >> ~/sw/crontab/node-exclude
printf "$(date +%F) %s\n" $(\sacct -j $JOB -s F,CA,OOM,TO -XnP -o NodeList |uniq) >> ~/sw/crontab/node-exclude-log
awk 'NF<=3 {print $2}' ~/sw/crontab/node-exclude-log > ~/sw/crontab/node-exclude

export CASENAME=`\sacct -j $JOB -s F,CA,OOM,TO -XnP -o JobName |sed -n 's/vint_//p' |uniq`
export SBATCH_ARRAY_INX=`\sacct -j $JOB -s F,CA,OOM,TO -XnP -o JobID |awk -F_ -vORS=, '{print $2}'`
#TODO JOB=`sbatch --account=kuang_lab -p huce_intel -J vintt_${CASENAME} -n 1 -x ~/sw/crontab/node-exclude -t 60 --mem=6000 -o "logfile/slurm.%x.%A.%a" --mail-type=FAIL,TIME_LIMIT --wrap='export tid=$SLURM_ARRAY_TASK_ID; while [ $tid -lt $(( $SLURM_ARRAY_TASK_ID + $iPerTask )) -a $tid -le $iMax ]; do  ncwa -O -a time,lon history/day$(printf "%05d" $(($tid *100+1000)))h00.1xday.nc analysis/zyorg_${CASENAME}_$(printf "%03d" $tid).nc; ncl -Q nc.fms.vint.ncl; let tid++; done' | egrep -o -e "\b[0-9]+$"`
EOF
: <<'EOF'
\sacct -XP |cut -d_ -f 1| sort -u
\sacct -XP --name=vint_b2kidstonens1 |grep -v completed
#\sacct -j 20023728,20023526,20023832 -s F,CA,OOM,TO -Xn -o NodeList%-20 |uniq >> ~/sw/crontab/node-exclude
\sacct -j 20023728,20023526,20023832 -s F,CA,OOM,TO -XnP -o NodeList |uniq >> ~/sw/crontab/node-exclude
srun -p huce_intel -w ~/sw/crontab/node-exclude ~/sw/crontab/mount-check.sh

JOB=20023728
JOB=22762936
export iMax=100
#export CASENAME=`\sacct -j $JOB -s F,CA,OOM,TO -Xn -o JobName%25 |sed -n 's/vint_//p' |uniq |xargs`
export CASENAME=`\sacct -j $JOB -s F,CA,OOM,TO -XnP -o JobName |sed -n 's/vint_//p' |uniq`
export SBATCH_ARRAY_INX=`\sacct -j $JOB -s F,CA,OOM,TO -XnP -o JobID |awk -F_ -vORS=, '{print $2}'`
JOB=`sbatch --account=kuang_lab -p huce_intel -J vint_${CASENAME} -n 1 -x ~/sw/crontab/node-exclude -t 60 --mem=6000 -o "logfile/slurm.%x.%a" --mail-type=FAIL,TIME_LIMIT --wrap='  >& logfile/slurm-${SLURM_JOB_NAME}-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}' | egrep -o -e "\b[0-9]+$"`

EOF

unset SBATCH_ARRAY_INX
sbatch --account=kuang_lab -p huce_bigmem,huce_intel -J 2nces_${CASENAME} -n 1 -x ~/sw/crontab/node-exclude -t 60 --mem=$((iMax*250+5000)) -o "slurm-%x" --mail-type=FAIL --dependency=afterok:${JOB} <<'EOF'
#!/bin/sh -e
set -x
 #ncra -O `\ls -1v analysis/var_${CASENAME}_*.nc` tmp.ncra.${CASENAME}.nc
 #ncwa -O -a lon tmp.ncra.${CASENAME}.nc zyvar_${CASENAME}.nc
#nces -O -n $iMax,3,1 analysis/zyvar_${CASENAME}_001.nc zyvar_${CASENAME}.nc
nces -O -n $iMax,3,1 analysis/zyorg_${CASENAME}_001.nc zyorg_${CASENAME}.nc
nces -O -n $iMax,3,1 analysis/nc.1cospectrum_001.nc nc.2cospectrum.nc
#ncl -Q ../../nc.cospectrum.ncl
#ncl -Q ../../nc.structure.ncl
touch finish
#touch finish-${CASENAME}

#ncdiff -O zyvar_${CASENAME}.nc zyvar_b1ctrl.nc diff_zy_${CASENAME}.nc
#rm cy.cospectrum.data.nc; ncl -Q cy.cospectrum.ncl
#rm yl.eke240.data.nc; ncl -Q yl.eke240.ncl
EOF

