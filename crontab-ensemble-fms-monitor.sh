#!/bin/sh

# 

#crontab -e  # and input the following, boslogin04
: <<'EOF'
SHELL=/bin/bash
PATH=/sbin:/bin:/usr/sbin:/usr/bin
MAILTO=""

# For details see man 4 crontabs

# Example of job definition:
# .---------------- minute (0 - 59)
# |  .------------- hour (0 - 23)
# |  |  .---------- day of month (1 - 31)
# |  |  |  .------- month (1 - 12) OR jan,feb,mar,apr ...
# |  |  |  |  .---- day of week (0 - 6) (Sunday=0 or 7) OR sun,mon,tue,wed,thu,fri,sat
# |  |  |  |  |
# *  *  *  *  *   command to be executed
 */10 * *  *  * /n/holylfs/LABS/kuang_lab/pchan/jetshift/crontab-ensemble-fms-monitor.sh >& /n/holylfs/LABS/kuang_lab/pchan/jetshift/log/crontab-log-`date '+%F-%R'`
EOF


 # per ensemble postprocess
 #echo c6nplus3exp1 10 c5ctrl >> crontab-running-ensemble
set -x
find /n/holylfs/LABS/kuang_lab/pchan/jetshift/log/crontab-log-* -mmin +1440 -size -40000c -delete
cd `dirname $0`
CASEprefixArr=($(awk '{print $1}' crontab-running-ensemble))
EnsSize=($(awk '{print $2}' crontab-running-ensemble))
CtrlCase=($(awk '{print $3}' crontab-running-ensemble))
for ii in ${!CASEprefixArr[*]}; do
#  if [ $(\ls zyvar_${CASEprefixArr[$ii]}ens*.nc |wc -l) -eq ${EnsSize[$ii]} ]; then
#  if [ $(\ls cy.cospectrum.${CASEprefixArr[$ii]}ens*.data.nc |wc -l) -eq ${EnsSize[$ii]} ]; then
#  if [ $(\ls finish-${CASEprefixArr[$ii]}ens* |wc -l) -eq ${EnsSize[$ii]} ]; then
  if [ $(\ls fms-output/${CASEprefixArr[$ii]}ens*/finish |wc -l) -eq ${EnsSize[$ii]} ]; then
    sed -i "/^${CASEprefixArr[$ii]}/d" crontab-running-ensemble
    source /usr/local/bin/lmod.sh
    module load intel/17.0.4-fasrc01 nco/4.7.4-fasrc01
    nces -O fms-output/${CASEprefixArr[$ii]}ens*/zyvar_${CASEprefixArr[$ii]}ens*.nc ensemble-wise/zyvar_${CASEprefixArr[$ii]}.nc
    nces -O fms-output/${CASEprefixArr[$ii]}ens*/zyorg_${CASEprefixArr[$ii]}ens*.nc ensemble-wise/zyorg_${CASEprefixArr[$ii]}.nc
#    nces -O fms-output/${CASEprefixArr[$ii]}ens*/cy.cospectrum.${CASEprefixArr[$ii]}ens*.data.nc ensemble-wise/nc.3cospectrum.${CASEprefixArr[$ii]}.nc
    ncdiff -O ensemble-wise/zyvar_${CASEprefixArr[$ii]}.nc ensemble-wise/zyvar_${CtrlCase[$ii]}.nc ensemble-wise/diff_zy_${CASEprefixArr[$ii]}.nc
    export CASENAME=${CASEprefixArr[$ii]}
    export NCARG_ROOT="/n/home05/pchan/sw/ncl-6.6.2"
    export PATH="$NCARG_ROOT/bin:$PATH"
    ncl -Q cy.cospectrum.ncl pp=300.
    ncl -Q cy.cospectrum.ncl pp=800.
    mail -s "Ensemble ${CASEprefixArr[$ii]} ended" ${USER} << EOF
EOF
  fi
done

