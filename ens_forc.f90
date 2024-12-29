 Program GenEnsForc
    ! PROGRAM FV3Tiles_To_Vector(LENS_OUT, vector_rand_ens)
    
    use stochastic_physics,  only : init_stochastic_physics_land,&
                            run_stochastic_physics_land, finalize_stochastic_physics
    use get_stochy_pattern_mod,  only : write_stoch_restart_atm

    use atmosphere_stub_mod, only: Atm,atmosphere_init_stub
    !use mpp_domains
    use mpp_mod,             only: mpp_set_current_pelist,mpp_get_current_pelist,mpp_init,mpp_pe,mpp_npes ,mpp_declare_pelist,mpp_root_pe
    use mpp_domains_mod,     only: mpp_broadcast_domain,MPP_DOMAIN_TIME,mpp_domains_init ,mpp_domains_set_stack_size
    use fms_mod,             only:  fms_init
    use xgrid_mod,           only: grid_box_type
    use netcdf
    use kinddef,             only : kind_dbl_prec,kind_phys
    use stochy_namelist_def, only : stochini, max_n_var_lndp
    !use M_DA, only: matsqrt     
    use netcdf
    use mpi

    implicit none

    integer, parameter :: dp = kind(1.0d+0)

    integer, parameter      :: nlevs=3
    ! integer, parameter :: max_n_var_lndp = 6
    integer                 :: ntasks,fid
    integer                 :: nthreads
    integer            :: ncid, xt_dim_id,yt_dim_id,time_dim_id,xt_var_id,&
                          yt_var_id,time_var_id,var_id_lat,var_id_lon,var_id_tile
    integer            :: varid1,varid2,varid3,varid4,varid_lon,varid_lat,varid_tile
    integer                 :: varidl(max_n_var_lndp)
    integer                 :: zt_dim_id,zt_var_id
    character*2             :: strid

    ! character(len=3), dimension(max_n_var_lndp)         ::  lndp_var_list
    ! real(kind=kind_dbl_prec), dimension(max_n_var_lndp) ::  lndp_prt_list
    
    real :: ak(nlevs+1),bk(nlevs+1)
    real(kind=4) :: ts,undef

    data ak(:) /0.0, 306.1489, 13687.72    , 0.99/
    data bk(:) /1.0,   0.9284,     0.013348, 0.0/
    integer     :: nb,blksz_1,nblks,ierr,my_id,i,j,k,l,nx,ny !,id
    integer     :: isc,iec,jsc,jec,isd,ied,jsd,jed
    integer     :: halo_update_type = 1
    real        :: dx,dy,pi,rd,cp
    ! logical   :: write_this_tile
    integer  :: nargs,ntile_out,nlunit,pe,npes,stackmax=4000000
    integer  :: i1,i2,j1,npts,istart,tpt
    ! character*80 :: fname
    ! character*1  :: ntile_out_str
    integer :: comm

!4.25.23 not sure below are being used
    real(kind=4), allocatable, dimension(:)                   :: grid_xt, grid_yt

    real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlat
    real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlon

    type(grid_box_type)                :: grid_box

    real (kind=kind_phys),allocatable  :: sfc_wts   (:,:,:)
    integer,allocatable                :: blksz(:)
    integer              :: me              !< MPI rank designator
    integer              :: root_pe         !< MPI rank of root atmosphere processor
    real(kind=kind_phys) :: dtp             !< physics timestep in seconds
    real(kind=kind_phys) :: fhour           !< previous forecast hour

    ! real(kind=kind_phys) :: sppt_amp        !< amplitude of sppt (to go to cld scheme)
    ! logical  :: do_sppt,do_shum,do_skeb,
    ! integer  ::  skeb_npass

    logical  :: use_zmtnblck    
    integer  ::  n_var_lndp, lndp_type
    character(len=65)              :: fn_nml                   !< namelist filename
    character(len=256),allocatable :: input_nml_file(:) !< character string containing full namelist
    
    Real, allocatable     :: vector_rand_ens(:)
    INTEGER :: IDIM, JDIM, NUM_TILES, IY, IM, ID, IH, &
               i_layout, j_layout, tot_subtiles, my_tile
    REAL    :: FH, DELTSFC
    ! INTEGER :: IERR
    INTEGER             :: NPROCS, MYRANK, NUM_THREADS, NUM_PARTHDS, MAX_TASKS
    REAL                :: horz_len_scale, ver_len_scale, temp_len_scale 
    Integer             :: ens_size, t_indx, t_len       !, n_surf_vars
    integer             :: n_forc_vars = 6
    integer             :: n_state_vars = 8
    ! LOGICAL             :: rcov_localize, ens_inflate
    CHARACTER(LEN=500)  :: static_filename, fv3_prefix, vector_prefix

    integer             :: PRINTRANK = 4
    LOGICAL             :: print_debg_info = .false.
  
    Integer                :: LENSFC, LENSFC_landm
    INTEGER, allocatable   :: tile_xy(:), Idim_xy(:), Jdim_xy(:), VETFCS(:)
    Real, allocatable      :: RLA_land(:), RLO_land(:), OROG_land(:)

    Real, allocatable      :: rand_Ta3D(:,:,:), sfc_wts_out(:,:,:)

    CHARACTER(LEN=500)     :: fv3_inp_file, vector_inp_file
    CHARACTER(LEN=2)       :: RANKCH 

    CHARACTER(LEN=2)       :: peg_str
    CHARACTER(LEN=4)       :: ensCH
    character(len=3)       :: mem_str
    CHARACTER(LEN=3)       :: proc_str
    CHARACTER(len=500)     :: forc_inp_file_ens
    CHARACTER(len=250)     :: forc_inp_path = "./"
    CHARACTER(len=250)     :: forc_inp_file = "C96_GDAS_forcing_2019-12-15.nc"    
    Integer                :: vector_size = 18360
    Real                   :: std_dev_f(6) = (/0.1, 0.05, 20.0, 2.0, 0.05, 0.05/)
    Real                   :: std_dev_s(8) = (/0.2, 0.2, 0.2, 0.2, 2.0, 2.0, 2.0, 2.0/)
    
    Integer                :: ixy, sxy, ipr, arr_indx, iproc, ie, iv, it
    INTEGER                :: mpiReal_size, rsize, isize, mpiInt_size
    integer, allocatable   :: tile_members(:,:), tile_group(:), comm_tile(:)
    integer                :: group_world  !, comm_tile !comm_world, 
    integer                :: Np_ext, Np_til, p_gN, p_gRank
    integer, allocatable   :: pg_bindx(:,:)

    LOGICAL                :: file_exists
    Integer                :: error, varid !ncid, 
 
    Real, allocatable      :: forcArray(:), stateArray(:, :)
    Real(dp), allocatable  :: C_cov(:,:), L_cov(:,:)

    character(len=500)      :: stoch_ini_file   ! input init pattern file

    logical                 :: perturb_forcing = .true., perturb_state = .false.    ! if true, state pert samples start at n_forc + 1 
    character(len=128)      :: state_file_name = "ufs_land_restart.2020-01-01_00-00-00.nc"
    CHARACTER(len=500)      :: state_file_ens
    Integer                 :: ncid_st

    Integer, dimension(6)  :: forc_ens_pert_type = (/1, 1, 0, 0, 1, 1/)
    Integer, dimension(8)  :: state_ens_pert_type = (/1, 1, 1, 1, 0, 0, 0, 0/)
! use snow and soil levels info 
    ! Integer, dimension(8)  :: state_layer_dims = (/3, 3, 4, 4, 1, 1, 1, 1/)
    Integer                :: st_layer_dim
    character(len=24)      :: st_layer_var

    character(len=128), dimension (6)   :: forc_var_list = &
        (/"precipitation", "solar_radiation",     &
        "longwave_radiation", "temperature",    &
        "wind_speed", "specific_humidity"/)
! ToDo Check smc or slc is the correct var

    character(len=128), dimension (8)   :: state_var_list =   &       
        (/'smc1', 'smc2', 'smc3', 'smc4', 'soilt1', 'soilt2', 'soilt3', 'soilt4'/)   !&
        ! (/"snow_level_ice", "temperature_snow",    &  !"snow_level_liquid",       &
        ! "soil_moisture_vol", "temperature_soil"/) 

! Forcing
! lrad and temp additive errors, rest multiplicative
! Threshold 
! !> Ensure downward longwave radition doesn’t have a negative (upward) values—due to its additive perturbation 

! ! States: Considered the following 
!     snow_level_ice(time, snow_levels, location)         mm
!     snow_level_liquid(time, snow_levels, location)      mm 
!     temperature_snow(time, snow_levels, location) ;
!     soil_liquid_vol(time, soil_levels, location) ;
!     soil_moisture_vol(time, soil_levels, location) volumetric moisture content m3/m3
!     temperature_soil(time, soil_levels, location)  soil level temp K
!
! Snow: forcing perturbation is enough to get sufficient variance. can in the future perturb snow_level_ice and temperature_snow
! Soil: perturb layer soil_moisture_vol, temperature_soil (smc and stc)
!       each layer var specified separately (to help handle correlation below)
        ! also correspond to the increment vars in GSI and land IAU 
        ! 'smc1', 'smc2', 'smc3', 'smc4', 'soilt1', 'soilt2', 'soilt3', 'soilt4'
! ! Correlations
!         1. Errors: multiplicative for snow mass and soil moisture, additive for temp
!         2. Threshold: ensure smc does not go beyond 1
!         3. Temp and Mass/moisture not correlated
!         4. snow and soil not correlated
!         5. Correlation between layers (for each variable)            
!             5.1 (first test) apply same random errors
!             5.2 generate same samples; scale by (inverse of ?) layer size. 
!                     To ensure proportionate error or to avoid unrealistically big errors on the lower layers?
!                     Ask Clara about this
!             5.3 specify correlation (high) and generate (correlated) individual samples
!             5.4 independent random errors
!  ! upper limits for smv   

    NAMELIST/NAMENS/ IDIM, JDIM, NUM_TILES, i_layout, j_layout, IY, IM, ID, IH, FH, DELTSFC, &
                    horz_len_scale, ver_len_scale, temp_len_scale, ens_size, &
                    t_len, t_indx, &
                    static_filename, fv3_prefix, vector_prefix, vector_size, &    !rand_var, &
                    PRINTRANK, print_debg_info, & 
                    perturb_forcing, forc_inp_path, forc_inp_file, n_forc_vars, forc_var_list, std_dev_f, forc_ens_pert_type, &
                    perturb_state, state_file_name, n_state_vars, state_var_list, std_dev_s, state_ens_pert_type               
    !
    DATA IDIM,JDIM,NUM_TILES, i_layout, j_layout/96,96,6, 1, 1/ 
    DATA IY,IM,ID,IH,FH/1997,8,2,0,0./
    DATA DELTSFC/0.0/, MAX_TASKS/99999/
    DATA horz_len_scale/55.0/
    DATA ver_len_scale/800./
    DATA temp_len_scale/24./
    DATA ens_size/20/
    Data t_len/24/
    DATA t_indx/2/
    
    DATA static_filename/""/
    DATA fv3_prefix/"./"/
    DATA vector_prefix/""/
    ! DATA rand_var/"smc"/
    ! Data vector_size/18360/
    ! DATA print_debg_info/.false./
    ! DATA PRINTRANK/4/
    ! Data n_forc_vars/6/
    ! Data forc_inp_path/"./"/
    ! Data forc_inp_file/"C96_GDAS_forcing_2019-12-15.nc"/   

    ! Data std_dev_f/0.1, 0.05, 20.0, 2.0, 0.05, 0.05/  
    ! Data perturb_state/.false./
    ! Data state_file_name/"ufs_land_restart.2020-01-01_00-00-00.nc"/
    ! Data n_state_vars/10/
    ! DATA num_assim_steps/1/  ! For multiple time steps of assimilation
    ! DATA dT_Asssim/24.0/     ! hrs. For multiple time steps of assimilation
    
    ! print*, "starting stanalone stochy"
    call fms_init()
    ! print*, "starting stanalone stochy 2"
    call mpp_init()
    ! print*, "starting stanalone stochy 3"
    call fms_init
    ! print*, "starting stanalone stochy 4"

    !CALL MPI_INIT(IERR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)
    call MPI_Comm_group(MPI_COMM_WORLD, group_world, IERR) 
    if (myrank==0) PRINT*," running GenEnsForc RANK WITH", NPROCS, "TASKS"
    ! PRINT*," running GenEnsForc  RANK ", MYRANK, " WITH ", NPROCS, "TASKS"

    ! define stuff
    ! pi=3.14159265359
    ! undef=9.99e+20
    ! p1000=100000.0

    ! !define mid-layer pressure
    ! rd=287.0
    ! cp=1004.0
    ! DO k=1,nlevs
    !     pressi(k)=ak(k)+p1000*bk(k)
    ! ENDDO
    ! ex3d=cp*(pressi/p1000)**(rd/cp)
    ! DO k=1,nlevs
    !     exn = (ex3d(k)*pressi(k)-ex3d(k+1)*pressi(k+1))/((cp+rd)*(pressi(k)-pressi(k+1)))
    !     pressl(k)=p1000*exn**(cp/rd)
    ! ENDDO 

    ! PRINT*,"READING NAMELIST."
    ! CALL BAOPENR(360, "tiles_to_vector.nml", IERR)             !"snowDA.nml"   !
    open(360, file="generate_ens_forc_state.nml", form="formatted")
    read(360, NAMENS)
    close(360)    
    IF (MYRANK==0) WRITE(6, NAMENS)
    ! LENSFC = IDIM*JDIM ! TOTAL NUMBER OF POINTS FOR THE CUBED-SPHERE TILE
    
    ! n_surf_vars = n_forc_vars + n_state_vars 
    ! n_var_lndp = n_surf_vars
    lndp_type = 2
    
    !   comm_world = comm_world
    my_id=mpp_pe()
    ntasks=mpp_npes()
    ! print*, "global scope"
    ! print*, "proc ", my_id, " of ", ntasks
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

    rsize = SIZEOF(horz_len_scale)
    Call MPI_TYPE_SIZE(MPI_REAL, mpiReal_size, IERR) 
    If (rsize == 4 ) then 
            mpiReal_size = MPI_REAL4
    elseif (rsize == 8 ) then 
            mpiReal_size = MPI_REAL8
    elseif (rsize == 16 ) then 
            mpiReal_size = MPI_REAL16
    else
            PRINT*," Possible mismatch between Fortran Real ", rsize," and Mpi Real ", mpiReal_size
            Stop
    endif
    isize = SIZEOF(vector_size) 
    Call MPI_TYPE_SIZE(MPI_INTEGER, mpiInt_size, IERR) 
    If (isize == 2 ) then 
        mpiInt_size = MPI_INTEGER2
    elseif (isize == 4 ) then 
        mpiInt_size = MPI_INTEGER4
    elseif (isize == 8 ) then 
        mpiInt_size = MPI_INTEGER8
    else
        PRINT*," Possible mismatch between Fortran Int ", isize," and Mpi Int ", mpiInt_size
        Stop
    endif
    
    tot_subtiles = NUM_TILES * i_layout * j_layout

    Np_ext = MOD(NPROCS, tot_subtiles)  ! extra/inactive procs
    if (MYRANK >  NPROCS - Np_ext - 1) goto 999

    Np_til = NPROCS / tot_subtiles  ! num tile groups    
    p_gN = MYRANK / tot_subtiles  ! group for proc.  
    p_gRank = MOD(MYRANK, tot_subtiles)  ! proc. rank within group
    
    Allocate(tile_members(Np_til, tot_subtiles))
    Allocate(tile_group(Np_til))
    Allocate(comm_tile(Np_til))  
      
    do ixy = 0, Np_til - 1   
        Do ipr = 0, tot_subtiles-1
            tile_members(ixy+1, ipr+1) = ixy*tot_subtiles + ipr   !(/ixy*NUM_TILES + ipr, ipr=0,NUM_TILES-1/)
        enddo
    enddo
    ! tile_members(1,:) = (/0, 1, 2, 3, 4, 5/)
    ! tile_members(2,:) = (/6, 7, 8, 9, 10, 11/)
    do ixy = 0, Np_til - 1   
        call MPI_Group_incl(group_world, tot_subtiles, tile_members(ixy+1,:), &
                           tile_group(ixy+1), IERR)
         call MPI_Comm_create(MPI_COMM_WORLD, tile_group(ixy+1), comm_tile(ixy+1), IERR)
    enddo    
    ! call MPI_Group_incl(group_world, NUM_TILES, tile_members(1,:), tile_group(1), IERR)
    ! call MPI_Group_incl(group_world, NUM_TILES, tile_members(2,:), tile_group(2), IERR)     
    ! call MPI_Comm_create(MPI_COMM_WORLD, tile_group(1), comm_tile(1), IERR)
    ! call MPI_Comm_create(MPI_COMM_WORLD, tile_group(2), comm_tile(2), IERR)
    if(myrank==0) print*, "Proc ",myrank, " done creating new tile group and comm"
    
    do ixy = 1, Np_til 
        WRITE(peg_str, '(I0)') ixy
        call mpp_declare_pelist( tile_members(ixy,:), "peList"//peg_str )
        ! WRITE(peg_str, '(I0)') 2
        ! call mpp_declare_pelist( tile_members(2,:), "peList"//peg_str )
    enddo
    
    do ixy = 0, Np_til - 1
        if (p_gN==ixy) then
            call mpp_set_current_pelist(tile_members(ixy+1,:))
            ! print*, "pelist ", ixy+1
            ntasks=mpp_npes()
            ! print*, "proc ", my_id, " of ", ntasks
            call atmosphere_init_stub (grid_box)
        endif
    Enddo 

    call mpp_set_current_pelist()

    ! CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    isd=Atm(1)%bd%isd
    ied=Atm(1)%bd%ied
    jsd=Atm(1)%bd%jsd
    jed=Atm(1)%bd%jed
    isc=Atm(1)%bd%isc
    iec=Atm(1)%bd%iec
    jsc=Atm(1)%bd%jsc
    jec=Atm(1)%bd%jec
    nx=iec-isc+1
    ny=jec-jsc+1

    if (p_gRank==0) then
        print*, "pe group", p_gN + 1
        print*,'nx, ny=',nx,ny
    endif
    
    my_tile = Atm(1)%tile_of_mosaic
    if (MYRANK ==0) print*, "size(Atm)", size(Atm), "npes_per_tile", Atm(1)%npes_per_tile, &
                    "tile_of_mosaic", my_tile  ! A Atm(1)%global_tile, 

    print*, "proc , isc,iec,jsc,jec,isd,ied,jsd,jed"
    print*, myrank, isc,iec,jsc,jec,isd,ied,jsd,jed

    Allocate(pg_bindx(tot_subtiles, 5))
    if (p_gRank /= 0) then  
        call MPI_SEND((/my_tile, isc,iec,jsc,jec/), 5, mpiInt_size, 0, &
                        10*(p_gRank+1), comm_tile(p_gN+1), IERR) 
    else
        pg_bindx(1,:) = (/my_tile, isc,iec,jsc,jec/)        
        Do iproc = 1, tot_subtiles-1 
            call MPI_RECV(pg_bindx(iproc+1,:), 5, mpiInt_size, &
            iproc, 10*(iproc+1), comm_tile(p_gN+1), MPI_STATUS_IGNORE, IERR)
        Enddo
    End if

    blksz_1=nx
    nblks=nx*ny/blksz_1
    allocate(blksz(nblks))
    do i=1,nblks
    blksz(i)=blksz_1
    enddo
    nthreads = 1
    me=my_id
    fhour=0
    dtp=3600            !900
    fn_nml='input.nml'
    nlunit=21
    
    !define model grid
    dx=360.0/nx
    dy=180.0/ny
    allocate(xlat(nblks,blksz_1))
    allocate(xlon(nblks,blksz_1))
    i1=isc
    j1=jsc
    do nb=1,nblks
        i2=i1+blksz_1-1
        if (i2 .le. iec) then 
        xlon(nb,1:blksz_1) = Atm(1)%gridstruct%agrid_64(i1:i2,j1,1)
        xlat(nb,1:blksz_1) = Atm(1)%gridstruct%agrid_64(i1:i2,j1,2)
        i1=i1+blksz_1
        else
        npts=iec-i1+1
        xlon(nb,1:npts) = Atm(1)%gridstruct%agrid_64(i1:iec,j1,1)
        xlat(nb,1:npts) = Atm(1)%gridstruct%agrid_64(i1:iec,j1,2)
        if (j1.LT. jec) then
            xlon(nb,npts+1:blksz_1) = Atm(1)%gridstruct%agrid_64(isc:isc+(blksz_1-npts+1),j1+1,1)
            xlat(nb,npts+1:blksz_1) = Atm(1)%gridstruct%agrid_64(isc:isc+(blksz_1-npts+1),j1+1,2)
        endif
        i1=npts+1
        j1=j1+1
        endif
        if (i2.EQ.iec) then
        i1=isc
        j1=j1+1
        endif
    end do
    
    if (p_gRank==0) then
        print*, "MYRANK",  MYRANK, "pe group", p_gN + 1, "min max xlon", &
        minval(xlon), maxval(xlon)
       
        print*, "MYRANK",  MYRANK, "pe group", p_gN + 1, "min max xlat", &
        minval(xlat), maxval(xlat)
    endif
    
    allocate(grid_xt(nx), grid_yt(ny))
    do i=1,nx
    grid_xt(i)=i
    enddo
    do j=1,ny
    grid_yt(j)=j
    enddo
    !     endif
    ! enddo  

    ! vector_size = LENS_OUT
    allocate(tile_xy(vector_size))
    allocate(Idim_xy(vector_size))
    allocate(Jdim_xy(vector_size))
    allocate(VETFCS(vector_size))
    allocate(RLA_land(vector_size))
    allocate(RLO_land(vector_size))
    allocate(OROG_land(vector_size))
    !   static_filename = "ufs-land_C"// trim(str(NDIM)) //"_static_fields.nc"
    call ReadTileInfo(trim(static_filename), vector_size, tile_xy, Idim_xy, Jdim_xy, &
                            RLA_land, RLO_land, OROG_land, VETFCS)              
    IF (MYRANK==0) then
        PRINT*, "finished reading tile info" !" proc ", myrank, 
        !  print*, "tile_xy", tile_xy
        !  print*, "Idim_xy", Idim_xy
        !  print*, "Jdim_xy", Jdim_xy
        ! print*, "vegetation_category", VETFCS_land   
        ! print*, "RLA", RLA_land
        ! print*, "RLO", RLO_land
        ! print*, "OROG", OROG_land
    Endif

    allocate(rand_Ta3D(NUM_TILES, JDIM, IDIM))
     
    !workg_T162_984x488.tile03.nc 
    allocate(forcArray(vector_size))  !, t_len))
    ! allocate(stateArray(vector_size, layer_dim_max))  !, t_len))
    allocate(vector_rand_ens(vector_size))  !, t_len))

    !apply inter-variable corr LL^T = C; Y = L*X
    allocate(C_cov(3,3)) !(n_surf_vars,n_surf_vars))
    allocate(L_cov(3,3)) !n_surf_vars,n_surf_vars))
    !P Qsi Qli Ta V RH Ps
    ! pert_corr_matrix = [[1.0, -0.8, 0.5, 0.0, 0.0, 0.0, 0.0],
    !                     [-0.8, 1.0, -0.5, 0.0, 0.0, 0.0, 0.0],
    !                     [0.5, -0.5, 1.0, 0.0, 0.0, 0.0, 0.0],
    !                     [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
    !                     [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    !                     [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
    !                     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]
    ! P Qsi Qli
    C_Cov(1,:) = (/1.0, -0.8, 0.5/)
    C_Cov(2,:) = (/-0.8, 1.0, -0.5/)
    C_Cov(3,:) = (/0.5, -0.5, 1.0/)
    !chol(C_Cov)
    L_Cov = matsqrt(C_Cov)         !calls lapack dpotrf(C_Cov)
    if (MYRANK == 0) then 
        print*, "C Cov matrix"
        print*, C_Cov
        print*, "Choleski factorization of corr matrix"
        print*, L_Cov
        print*, " set upper L to 0"
        L_Cov(1,2:3) = (/0.0, 0.0/)
        L_Cov(2,3) = 0.0
        print*, L_Cov
        print*, "LL**T"
        print*, matmul(L_Cov, transpose(L_Cov))
    endif
    
    ! print*,'calling init_stochastic_physics',nlevs
    allocate(input_nml_file(1))
    input_nml_file='input.nml'
    ! do_sppt=.false.
    ! do_shum=.false.
    ! do_skeb=.false.

    comm = comm_tile(p_gN+1)    !MPI_COMM_WORLD
    root_pe = 0  !if (p_gRank == 0)  root_pe = myrank  ! mpp_root_pe()   
    
    ! Np_til = NPROCS / NUM_TILES  ! num tile groups    
    ! p_gN = MYRANK / NUM_TILES  ! group for proc.  
    ! p_gRank = MOD(MYRANK, NUM_TILES)  ! proc. rank within group
    
    if (perturb_forcing) then 

        n_var_lndp = n_forc_vars
        allocate(sfc_wts(nblks, blksz_1, n_var_lndp))
        allocate(sfc_wts_out(JDIM / j_layout, IDIM / i_layout, n_var_lndp))

        ! WRITE(proc_str, '(I0)') myrank
        Do ie = p_gN*(ens_size/Np_til)+1, (p_gN+1)*(ens_size/Np_til)  !ens_size
            WRITE(ensCH, '(I0)') ie
            write(mem_str, '(I3.3)') ie
            stoch_ini_file = TRIM(forc_inp_path)//'/RESTART/stochy_final_ens_forc'//TRIM(ensCH)//'.nc'
            if (p_gRank == 0) then  
                forc_inp_file_ens=TRIM(forc_inp_path)//"/mem"//TRIM(mem_str)//"/"//TRIM(forc_inp_file)
                ! forcing file
                INQUIRE(FILE=trim(forc_inp_file_ens), EXIST=file_exists)
                if (.not. file_exists) then 
                    print *, 'error,file does not exist ', trim(forc_inp_file_ens) , ' exiting'
                    call MPI_ABORT(MPI_COMM_WORLD, 10, error) 
                endif
                error = nf90_open(trim(forc_inp_file_ens), NF90_WRITE, ncid)
                call netcdf_err(error, 'opening forcing file' )
            endif

            call init_stochastic_physics_land(nlevs, blksz, dtp,              &   !sppt_amp,                         &
                input_nml_file, stoch_ini_file, fn_nml, nlunit, xlon, xlat,   &  ! do_sppt, do_shum, do_skeb, lndp_type,                &
                n_var_lndp, use_zmtnblck,   &  !skeb_npass, lndp_var_list, lndp_prt_list,    &
                ak, bk, nthreads, root_pe, comm_tile(p_gN+1), ierr)
            if (ierr .ne. 0) then 
                print *, 'ERROR init_stochastic_physics call'
                call MPI_ABORT(MPI_COMM_WORLD, IERR, error)
            endif
            !if (p_gRank == 0) then
            !    print*, " pe group", p_gN + 1, " Done init_stochy, ens ", ie, " time ", it    
            !endif
            
            Do it=1, t_len
                call run_stochastic_physics_land(nlevs, it-1, fhour, blksz, nthreads=nthreads, sfc_wts=sfc_wts)
                            ! sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts, skebv_wts=skebv_wts,                         
                ! run_stochastic_physics_land(levs, kdt, fhour, blksz, nthreads, sfc_wts) 

                ! do l=1,n_var_lndp
                !     do j=1,ny
                !         sfc_wts_out(l,:,j)=sfc_wts(j,:,l)
                !     enddo
                ! enddo

                !apply inter-variable corr L on P Qsi Qli Y = L*X; LL^T = C
                sfc_wts_out =  sfc_wts
                sfc_wts_out(:,:,2) =  L_cov(2,1) * sfc_wts(:,:,1) + &
                                    L_cov(2,2) * sfc_wts(:,:,2) !+ &
                                    !   L_cov(2,3) * sfc_wts(:,:,3)  L_cov(2,3)=0.
                sfc_wts_out(:,:,3) =  L_cov(3,1) * sfc_wts(:,:,1) + &
                                    L_cov(3,2) * sfc_wts(:,:,2) + &
                                    L_cov(3,3) * sfc_wts(:,:,3)

                !if (p_gRank == 0) then
                !    print*, " pe group", p_gN + 1, " Done run_stochy, ens ", ie, " time ", it    
                !endif

                Do ixy = 1, n_var_lndp        !vector_length 
                    if (p_gRank /= 0) then  
                        call MPI_SEND(sfc_wts_out(:,:,ixy), JDIM*IDIM/(j_layout*i_layout), mpiReal_size, 0, &
                                        10*(p_gRank+1)+ixy, comm_tile(p_gN+1), IERR) 
                    else
    !rand_Ta3D(6, JDIM * j_layout, IDIM * i_layout))
                        ! isc,iec,jsc,jec             
                        rand_Ta3D(my_tile,jsc:jec,isc:iec) = sfc_wts_out(:,:,ixy)
                        Do iproc = 1, tot_subtiles-1 
                            ! call MPI_RECV(rand_Ta3D(iproc+1,:,:),JDIM*IDIM, mpiReal_size, &
                            call MPI_RECV(rand_Ta3D(pg_bindx(iproc+1,1), pg_bindx(iproc+1,2):pg_bindx(iproc+1,3), &
                                                    pg_bindx(iproc+1,4):pg_bindx(iproc+1,5)), &
                                                    JDIM*IDIM/(j_layout*i_layout), mpiReal_size, &                        
                            iproc, 10*(iproc+1)+ixy, comm_tile(p_gN+1), MPI_STATUS_IGNORE, IERR)
                        Enddo
                        Do iv=1, vector_size 
                            vector_rand_ens(iv) = rand_Ta3D(tile_xy(iv), Jdim_xy(iv), Idim_xy(iv))
                        end do

                        ! Start writing forcing file
                        error = nf90_inq_varid(ncid, trim(forc_var_list(ixy)), varid)
                        call netcdf_err(error, 'getting varid '//trim(forc_var_list(ixy)) )
                        !read
                        ERROR=NF90_GET_VAR(NCID, varid, forcArray, &
                                start = (/1,it/), count = (/vector_size, 1/))
                        CALL NETCDF_ERR(ERROR, 'reading values for '//trim(forc_var_list(ixy)) )                  
                        ! print*, trim(forc_var_list(ixy))
                        ! print*, forcArray
                        ! print*, "rand pattern"
                        ! print*, vector_rand_ens
                        if (forc_ens_pert_type(ixy) == 1) then 
                        ! pert_factors = 1 + std_dev * corr_rand
                            vector_rand_ens = 1.0 + (std_dev_f(ixy) * vector_rand_ens)
                            forcArray = vector_rand_ens * forcArray 
                        else
                            ! pert_factors = target_mean(0.0) + std_dev * corr_rand
                            vector_rand_ens = std_dev_f(ixy) * vector_rand_ens
                            forcArray = vector_rand_ens + forcArray
                        endif
                        ! print*, trim(forc_var_list(ixy)), " after mult"
                        ! print*, forcArray

! ensure downward longwave rad doesn't have negative values
! note this revision if the different name used for lwrad
                        if (trim(forc_var_list(ixy)) .eq. "longwave_radiation") then
                            Where(forcArray < 0) forcArray = 0
                        endif

                        error = nf90_put_var(ncid, varid , forcArray,       &
                            start = (/1,it/), count = (/vector_size, 1/))
                        call netcdf_err(error, 'writing '//trim(forc_var_list(ixy)) )
                    Endif                
                end do  
            Enddo
            if (p_gRank == 0) then
                error = nf90_close(ncid)
                call netcdf_err(error, 'closing '//trim(forc_inp_file_ens) ) 
            endif
            ! if (stochini) then
            !     call write_stoch_restart_atm('stochy_final_2_ens'//ensCH//'.nc')
            !  else
            call write_stoch_restart_atm(stoch_ini_file)
            !  endif
        
            ! if (p_gRank == 0) PRINT*,"proc ", p_gRank, " finished copying rand for ens ", ie

            call finalize_stochastic_physics()
        Enddo

        if (allocated(sfc_wts)) deallocate(sfc_wts)
        if (allocated(sfc_wts_out)) deallocate(sfc_wts_out)

    endif
    
    if(perturb_state) then      !.and. (it .eq. 1)

        n_var_lndp = n_state_vars
        allocate(sfc_wts(nblks, blksz_1, n_var_lndp))
        allocate(sfc_wts_out(JDIM / j_layout, IDIM / i_layout, n_var_lndp))

        Do ie = p_gN*(ens_size/Np_til)+1, (p_gN+1)*(ens_size/Np_til)  !ens_size

            WRITE(ensCH, '(I0)') ie
            write(mem_str, '(I3.3)') ie
            stoch_ini_file = TRIM(forc_inp_path)//'/RESTART/stochy_final_ens_state'//TRIM(ensCH)//'.nc'
            if (p_gRank == 0) then  
                ! state file            
                state_file_ens=TRIM(forc_inp_path)//"/mem"//TRIM(mem_str)//"/"//TRIM(state_file_name)
                INQUIRE(FILE=trim(state_file_ens), EXIST=file_exists)
                if (.not. file_exists) then
                    print *, 'error,file does not exist', trim(state_file_ens) , ' exiting'
                    call MPI_ABORT(MPI_COMM_WORLD, 10, error)
                endif
                error = nf90_open(trim(state_file_ens), NF90_WRITE, ncid_st)
                call netcdf_err(error, 'opening state file '//trim(state_file_ens) )
            endif

            call init_stochastic_physics_land(nlevs, blksz, dtp,              &       !sppt_amp,                         &
                input_nml_file, stoch_ini_file, fn_nml, nlunit, xlon, xlat,   &  ! do_sppt, do_shum, do_skeb, lndp_type,                &
                n_var_lndp, use_zmtnblck,   &  !skeb_npass, lndp_var_list, lndp_prt_list,    &
                ak, bk, nthreads, root_pe, comm_tile(p_gN+1), ierr)
            if (ierr .ne. 0) then 
                print *, 'ERROR init_stochastic_physics call'
                call MPI_ABORT(MPI_COMM_WORLD, IERR, error)
            endif
            !if (p_gRank == 0) then
            !    print*, " pe group", p_gN + 1, " Done init_stochy, ens ", ie, " time ", it    
            !endif
            
            it = 1   ! perturb only first time step state
            call run_stochastic_physics_land(nlevs, it-1, fhour, blksz, nthreads=nthreads, sfc_wts=sfc_wts)
                        ! sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts, skebv_wts=skebv_wts,                         
            ! run_stochastic_physics_land(levs, kdt, fhour, blksz, nthreads, sfc_wts) 

            ! do l=1,n_var_lndp
            !     do j=1,ny
            !         sfc_wts_out(l,:,j)=sfc_wts(j,:,l)
            !     enddo
            ! enddo

    !>>>> State correlations go here
    !>>>>
            ! !apply inter-variable corr L on P Qsi Qli Y = L*X; LL^T = C
            ! sfc_wts_out =  sfc_wts
            ! sfc_wts_out(:,:,2) =  L_cov(2,1) * sfc_wts(:,:,1) + &
            !                       L_cov(2,2) * sfc_wts(:,:,2) !+ &
            !                     !   L_cov(2,3) * sfc_wts(:,:,3)  L_cov(2,3)=0.
            ! sfc_wts_out(:,:,3) =  L_cov(3,1) * sfc_wts(:,:,1) + &
            !                       L_cov(3,2) * sfc_wts(:,:,2) + &
            !                       L_cov(3,3) * sfc_wts(:,:,3)

            !if (p_gRank == 0) then
            !    print*, " pe group", p_gN + 1, " Done run_stochy, ens ", ie, " time ", it    
            !endif  

            Do ixy = 1, n_var_lndp        !vector_length 
                if (p_gRank /= 0) then  
                    call MPI_SEND(sfc_wts_out(:,:,ixy), JDIM*IDIM/(j_layout*i_layout), mpiReal_size, 0, &
                                    10*(p_gRank+1)+ixy, comm_tile(p_gN+1), IERR) 
                else           
                    rand_Ta3D(my_tile,jsc:jec,isc:iec) = sfc_wts_out(:,:,ixy)
                    Do iproc = 1, tot_subtiles-1 
                        ! call MPI_RECV(rand_Ta3D(iproc+1,:,:),JDIM*IDIM, mpiReal_size, &
                        call MPI_RECV(rand_Ta3D(pg_bindx(iproc+1,1), pg_bindx(iproc+1,2):pg_bindx(iproc+1,3), &
                                                pg_bindx(iproc+1,4):pg_bindx(iproc+1,5)), &
                                                JDIM*IDIM/(j_layout*i_layout), mpiReal_size, &                        
                        iproc, 10*(iproc+1)+ixy, comm_tile(p_gN+1), MPI_STATUS_IGNORE, IERR)
                    Enddo
                    Do iv=1, vector_size 
                        vector_rand_ens(iv) = rand_Ta3D(tile_xy(iv), Jdim_xy(iv), Idim_xy(iv))
                    end do

                    ! get var name
                    call map_layer_var_names(state_var_list(ixy), st_layer_var, st_layer_dim)

                    error = nf90_inq_varid(ncid_st, trim(st_layer_var), varid)
                    call netcdf_err(error, 'getting varid '//trim(st_layer_var) )

                    !read
                    ERROR=NF90_GET_VAR(ncid_st, varid, forcArray, &
                            start = (/1, st_layer_dim, it/), count = (/vector_size, 1, 1/))
                    CALL NETCDF_ERR(ERROR, 'reading values for '//trim(state_var_list(ixy)) )                 

                    if (state_ens_pert_type(ixy) == 1) then 
                        ! pert_factors = 1 + std_dev * corr_rand
                        vector_rand_ens = 1.0 + (std_dev_f(ixy) * vector_rand_ens)
                        forcArray = vector_rand_ens * forcArray 
!Note: assuming only vol soil moist content with unit of mm3/mm3 etc
! revise if var like snow content is perturbed
                        Where(forcArray > 1.) forcArray = 1.
                    else
                        ! pert_factors = target_mean(0.0) + std_dev * corr_rand
                        vector_rand_ens = std_dev_f(ixy) * vector_rand_ens
                        forcArray = vector_rand_ens + forcArray
                    endif

                    error = nf90_put_var(ncid_st, varid , forcArray,       &
                    start = (/1, st_layer_dim, it/), count = (/vector_size, 1, 1/))
                    call netcdf_err(error, 'writing '//trim(state_var_list(ixy)) )
                Endif       

            end do         

            if (p_gRank == 0) then
                error = nf90_close(ncid_st)
                call netcdf_err(error, 'closing '//trim(state_file_ens) )
            endif
            ! if (stochini) then
            !     call write_stoch_restart_atm('stochy_final_2_ens'//ensCH//'.nc')
            !  else
            call write_stoch_restart_atm(stoch_ini_file)
            !  endif
        
            ! if (p_gRank == 0) PRINT*,"proc ", p_gRank, " finished copying rand for ens ", ie

            call finalize_stochastic_physics()

        Enddo

        if (allocated(sfc_wts)) deallocate(sfc_wts)
        if (allocated(sfc_wts_out)) deallocate(sfc_wts_out)

    End if


    
    Deallocate(tile_members)
    Deallocate(tile_group)
    Deallocate(comm_tile)

    DEAllocate(pg_bindx)

    deallocate(tile_xy)
    deallocate(Idim_xy)
    deallocate(Jdim_xy)
    deallocate(RLA_land)
    deallocate(RLO_land)
    deallocate(OROG_land)
    deallocate(VETFCS)
    ! deallocate(rand_Ta, rand_Prec, rand_sRad)
    ! deallocate(rand_lRad, rand_sHum, rand_wSpeed)    
    deallocate(rand_Ta3D)
    deallocate(forcArray)
    if (allocated(stateArray)) deallocate(stateArray)
    deallocate(vector_rand_ens)
    deallocate(C_cov, L_cov)    
   
999 CONTINUE
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    PRINT*
    if (p_gRank == 0) PRINT*,'ensForc COMPLETED ON rank group:', p_gN + 1

    CALL MPI_FINALIZE(IERR)

    ! STOP
    contains
        !http://fortranwiki.org/fortran/show/
        ! ! LAPACK.
        function matsqrt(A) result(Asqrt)
            
            IMPLICIT NONE
            
            integer, parameter :: dp = kind(1.0d+0)

            real(dp), dimension(:,:), intent(in) :: A
            real(dp), dimension(size(A,1),size(A,2)) :: Asqrt
            integer :: N, info
        
            ! External procedures defined in LAPACK
            external DPOTRF
        
            ! Backup
            Asqrt = A
            N = size(A,1)
        
            call DPOTRF ('L', N, Asqrt, N, info)
        
            if (info /= 0) then
                stop 'Matrix factorization failed!'
            end if

        end function matsqrt
        
        subroutine map_layer_var_names(in_var_name, out_var_name, layer_dim)  !, error)
            implicit none
            character (len=*),   intent(in) :: in_var_name
            character (len=24), intent(out) :: out_var_name
            integer,            intent(out) :: layer_dim
            ! integer,          intent(out) :: error
            integer   :: iostat
            character :: smcname(4), stcname(6)

        select case (trim(in_var_name))

            case('smc1', 'smc2', 'smc3', 'smc4')
                smcname = trim(in_var_name)
                out_var_name = "soil_moisture_vol"
                read(smcname(4:4), *, iostat=iostat)  layer_dim                
                if (iostat /= 0) then 
                    print*, "error getting layer information from "//trim(in_var_name)
                    call MPI_ABORT(MPI_COMM_WORLD, IERR, error)
                endif
                if (myrank == 0) print*,"doing state perturbation for "//trim(in_var_name)

            case('soilt1', 'soilt2', 'soilt3', 'soilt4') 
                stcname = trim(in_var_name)
                out_var_name = "temperature_soil"
                read(stcname(6:6), *, iostat=iostat)  layer_dim                
                if (iostat /= 0) then 
                    print*, "error getting layer information from "//trim(in_var_name)
                    call MPI_ABORT(MPI_COMM_WORLD, IERR, error)
                endif
                if (myrank == 0) print*,"doing state perturbation for "//trim(in_var_name)

            case default
                print*, 'ERROR: unknown land state variable. Variable must be one of :',   &
                        'smc1', 'smc2', 'smc3', 'smc4', 'soilt1', 'soilt2', 'soilt3', 'soilt4'
            
                call MPI_ABORT(MPI_COMM_WORLD, IERR, error)

        end select 

        end subroutine map_layer_var_names

 END Program GenEnsForc

 subroutine ReadTileInfo(filename, vector_length, tile_xy, Idim_xy, Jdim_xy, &
                        RLA, RLO, OROG, VETFCS) ! vegetation_category)

        use netcdf

        ! type(noahmp_type) :: noahmp
        CHARACTER(LEN=*)  :: filename
        integer           :: vector_length, vector_length_in
        integer           :: ncid, dimid, varid, status
        real              :: RLA(vector_length), RLO(vector_length), OROG(vector_length)
        Integer  :: tile_xy(vector_length), Idim_xy(vector_length), &
                    Jdim_xy(vector_length), VETFCS(vector_length) !, vegetation_category(:)
        LOGICAL                   :: file_exists

        INQUIRE(FILE=trim(filename), EXIST=file_exists)
        if (.not. file_exists) then 
            print *, 'error,file does not exist', &   
                    trim(filename) , ' exiting'
            ! call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
            STOP
        endif           

        status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
        CALL NETCDF_ERR(status, 'opening file '//trim(filename))

        status = nf90_inq_dimid(ncid, "location", dimid)
        status = nf90_inquire_dimension(ncid, dimid, len = vector_length_in)
        
        if (vector_length /= vector_length_in) then
            print*, "wrong vector size, stop"
            stop
        endif
        ! allocate(tile_xy(vector_length))
        ! allocate(Idim_xy(vector_length))
        ! allocate(Jdim_xy(vector_length))
        ! allocate(VETFCS(vector_length))
        ! allocate(RLA(vector_length))
        ! allocate(RLO(vector_length))
        ! allocate(OROG(vector_length))

        status = nf90_inq_varid(ncid, "cube_tile", varid)
        ! if(status /= nf90_noerr) call handle_err(status)
        CALL NETCDF_ERR(status, 'reading cube_tile var id' )
        status = nf90_get_var(ncid, varid, tile_xy, &
            start = (/1/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading cube tile value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "cube_i", varid)
        CALL NETCDF_ERR(status, 'reading cube i vari id' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, Idim_xy, &
            start = (/1/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading cube i value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "cube_j", varid)
        CALL NETCDF_ERR(status, 'reading cube j vari id' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, Jdim_xy, &
            start = (/1/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading cube j value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "vegetation_category", varid)
        CALL NETCDF_ERR(status, 'reading vegetation_category vari id' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, VETFCS, &
            start = (/1/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading vegetation_category value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "latitude", varid)
        ! if(status /= nf90_noerr) call NETCDF_ERR(status)
        status = nf90_get_var(ncid, varid, RLA, &
            start = (/1/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading latitude value' )
        ! if(status /= nf90_noerr) call NETCDF_ERR(status)

        status = nf90_inq_varid(ncid, "longitude", varid)
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, RLO, &
            start = (/1/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading longitude value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "elevation", varid)
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, OROG, &
            start = (/1/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading elevation value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        ! status = nf90_inq_varid(ncid, "land_mask", varid)
        ! if(status /= nf90_noerr) call handle_err(status)
        ! status = nf90_get_var(ncid, varid, this%land_mask, &
        !     start = (/namelist%begsub/), count = (/namelist%lensub/))
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_close(ncid)
        CALL NETCDF_ERR(status, 'closing file' )

    end subroutine ReadTileInfo

    subroutine Read_FV3_File(filename, var_name, IDIM, JDIM, t_indx, Var_Out)

        use netcdf

        ! type(noahmp_type) :: noahmp
        CHARACTER(LEN=*)  :: filename, var_name
        integer           :: IDIM, JDIM, t_indx
        Real              :: Var_Out(IDIM, JDIM)   !, DUMMY(IDIM, JDIM)
        ! real, allocatable :: RLA(:), RLO(:), OROG(:)
        LOGICAL                   :: file_exists
        integer           :: ncid, dimid, varid, status

        INQUIRE(FILE=trim(filename), EXIST=file_exists)
        if (.not. file_exists) then 
            print *, 'error,file does not exist', &   
                    trim(filename) , ' exiting'
            ! call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
            STOP
        endif           

        status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
        CALL NETCDF_ERR(status, 'opening file '//trim(filename))

        status = nf90_inq_varid(ncid, trim(var_name), varid)
        ! if(status /= nf90_noerr) call handle_err(status)
        CALL NETCDF_ERR(status, 'reading '//trim(var_name)//' var id' )
        status = nf90_get_var(ncid, varid, Var_Out, &
            start = (/1, 1, t_indx/), count = (/IDIM, JDIM, 1/))
        CALL NETCDF_ERR(status, 'reading cube tile value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        ! Var_Out = RESHAPE(DUMMY, (/IDIM * JDIM/))   

        status = nf90_close(ncid)
        CALL NETCDF_ERR(status, 'closing file '//trim(filename) )

    end subroutine Read_FV3_File

    SUBROUTINE NETCDF_ERR( ERR, STRING )

    !--------------------------------------------------------------
    ! IF AT NETCDF CALL RETURNS AN ERROR, PRINT OUT A MESSAGE
    ! AND STOP PROCESSING.
    !--------------------------------------------------------------
        
        use netcdf

        IMPLICIT NONE
        ! include 'mpif.h'

        INTEGER, INTENT(IN) :: ERR
        CHARACTER(LEN=*), INTENT(IN) :: STRING
        CHARACTER(LEN=80) :: ERRMSG

        IF( ERR == NF90_NOERR )RETURN
        ERRMSG = NF90_STRERROR(ERR)
        PRINT*,''
        PRINT*,'FATAL ERROR: ', TRIM(STRING), ': ', TRIM(ERRMSG)
        PRINT*,'STOP.'
        ! CALL MPI_ABORT(MPI_COMM_WORLD, 999)
        STOP

        RETURN
    END SUBROUTINE NETCDF_ERR
