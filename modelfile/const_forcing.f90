
module const_forcing_mod

!-----------------------------------------------------------------------

use     constants_mod, only: KAPPA, CP_AIR, GRAV, RDGAS, PI

use           fms_mod, only: error_mesg, FATAL, file_exist,       &
                             open_namelist_file, set_domain,      &
			     read_data, check_nml_error,          &
                             mpp_pe, mpp_root_pe, close_file,     &
                             write_version_number, stdlog,        &
                             uppercase,&  !pjk
                             mpp_clock_id,mpp_clock_begin,mpp_clock_end,CLOCK_COMPONENT!,mpp_chksum

use  time_manager_mod, only: time_type, get_time

use  diag_manager_mod, only: register_diag_field, send_data

use  field_manager_mod, only: MODEL_ATMOS, parse
use tracer_manager_mod, only: query_method, get_tracer_index

use transforms_mod, only:     grid_domain

implicit none
private

!-----------------------------------------------------------------------
!---------- interfaces ------------

   public :: const_forcing, const_forcing_init

!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

   real :: p_k   ! pressure of forcing maximum
   real :: lat_j ! latitude of forcing maximum
   real :: p_hw = 7500.0 ! half-width of forcing (pressure)
   real :: lat_hw = 10.0   ! half-width of forcing (latitude)
   real :: u_amp = 5.787e-7 ! amplitude of u forcing (m/s/s) -- default = 0.05m/s/day
   real :: t_amp = 5.787e-7 ! amplitude of t forcing (k/s) -- default = 0.05K/day
   logical :: forcing_t = .false. ! should t have constant forcing applied
   logical :: forcing_u = .false. ! should u have constant forcing applied
!-----------------------------------------------------------------------

   namelist /const_forcing_nml/ p_k, lat_j, u_amp, t_amp, forcing_t, forcing_u, p_hw, lat_hw

!-----------------------------------------------------------------------

   character(len=128) :: version='$Id: const_forcing.f90,v 10.0 2003/10/27 23:31:04 arl Exp $'
   character(len=128) :: tagname='$Name: lima $'

   integer :: id_udt_cf, id_tdt_cf

   real    :: missing_value = -1.e10
   character(len=14) :: mod_name = 'const_forcing'

contains

!#######################################################################

subroutine const_forcing (is, js, Time, lat, p_full, udt, tdt )

!-----------------------------------------------------------------------
  type(time_type), intent(in)                  :: Time
  integer, intent(in) :: is, js
  real, intent(in),    dimension(:)     :: lat
  real, intent(in),    dimension(:,:,:)   :: p_full
  real, intent(inout), dimension(:,:,:)   :: udt, tdt
!-----------------------------------------------------------------------

   integer :: i, j, k, k2, j2, np, nj 
   logical :: used
   real, dimension(size(p_full,2), size(p_full,3)) :: su, sT
   real, dimension(10, 10) :: Au, AT   
   real, dimension(size(udt,1), size(udt,2), size(udt,3)) :: &
        tudt, ttdt

!-----------------------------------------------------------------------
   open (unit=20, file='ForcingM2EOF.txt', status='old', action='read')

   do np=1,10
     do nj=1,10
       read(20, *), AT(nj,np), Au(nj,np)
     enddo     
   enddo
   
   close(20)
            
   tudt(:,:,:) = 0.0
   ttdt(:,:,:) = 0.0   
   do k2=1,10
      do j2=1,10
         su(:,:) = 0.0
         sT(:,:) = 0.0
         do j=1,size(p_full,2)
            do k=1,size(p_full,3)
              sT(j,k)=exp(-((p_full(1,j,k)-k2*10000.0)/p_hw)**2-((abs(lat(j))-(j2-1)*10.0)/lat_hw)**2)
              if (k2==1) then
                su(j,k)=exp(-((p_full(1,j,k)-k2*10000.0)/p_hw)**2-((abs(lat(j))-(j2-1)*10.0)/lat_hw)**2)
                su(j,1)= 0.0
              else 
                su(j,k)=exp(-((p_full(1,j,k)-k2*10000.0)/p_hw)**2-((abs(lat(j))-(j2-1)*10.0)/lat_hw)**2)
              end if
            enddo
         enddo
         
         do j=1,size(p_full,2)
            do k=1,size(p_full,3)
                tudt(:,j,k) = tudt(:,j,k) + u_amp * Au(j2,k2) * su(j,k) / (86400.0)
                ttdt(:,j,k) = ttdt(:,j,k) + T_amp * AT(j2,k2) * sT(j,k) / (86400.0)  
            enddo
         enddo
         
      enddo
   enddo

   udt = udt + tudt
   if (id_udt_cf > 0) used = send_data ( id_udt_cf, tudt, Time, is, js)

   tdt = tdt + ttdt
   if (id_tdt_cf > 0) used = send_data ( id_tdt_cf, ttdt, Time, is, js)
            
!-----------------------------------------------------------------------

    end subroutine const_forcing

!#######################################################################

 subroutine const_forcing_init ( axes, Time )

!-----------------------------------------------------------------------
!
!           routine for initializing the model with an
!              initial condition at rest (u & v = 0)
!
!-----------------------------------------------------------------------

           integer, intent(in) :: axes(4)
   type(time_type), intent(in) :: Time

!-----------------------------------------------------------------------
   integer  unit, io, ierr
! 

!     ----- read namelist -----

      if (file_exist('input.nml')) then
         unit = open_namelist_file ( )
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=const_forcing_nml, iostat=io, end=10)
            ierr = check_nml_error (io, 'const_forcing_nml')
         enddo
  10     call close_file (unit)
      endif

!     ----- write version info and namelist to log file -----

      call write_version_number (version,tagname)
      if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=const_forcing_nml)

!     ----- register diagnostic fields -----

      if(forcing_u) then
         id_udt_cf = register_diag_field ( mod_name, 'udt_forcing', axes(1:3), Time, &
              'zonal wind forcing', 'm/s/s'   , &
              missing_value=missing_value)
      end if

      if(forcing_t) then
         id_tdt_cf = register_diag_field ( mod_name, 'tdt_forcing', axes(1:3), Time, &
              'temperature forcing', 'k/s'   , &
              missing_value=missing_value)
      end if

!-----------------------------------------------------------------------

    end subroutine const_forcing_init

end module const_forcing_mod
