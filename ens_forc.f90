 Program GenEnsForc
    ! PROGRAM FV3Tiles_To_Vector(LENS_OUT, vector_rand_ens)
    
    use stochastic_physics,  only : init_stochastic_physics,&
                            run_stochastic_physics, finalize_stochastic_physics
    use get_stochy_pattern_mod,  only : write_stoch_restart_atm

    use atmosphere_stub_mod, only: Atm,atmosphere_init_stub
    !use mpp_domains
    use mpp_mod,             only: mpp_set_current_pelist,mpp_get_current_pelist,mpp_init,mpp_pe,mpp_npes ,mpp_declare_pelist,mpp_root_pe
    use mpp_domains_mod,     only: mpp_broadcast_domain,MPP_DOMAIN_TIME,mpp_domains_init ,mpp_domains_set_stack_size
    use fms_mod,             only:  fms_init
    use xgrid_mod,           only: grid_box_type
    use netcdf
    use kinddef,             only : kind_dbl_prec,kind_phys
    use stochy_namelist_def, only : stochini
    !use M_DA, only: matsqrt     
    use netcdf
    use mpi

    implicit none
    ! include 'mpif.h'
    ! include 'netcdf.inc'
     !USE intrinsic::ieee_arithmetic
    integer, parameter :: dp = kind(1.0d+0)

    integer, parameter      :: nlevs=3
    integer, parameter :: max_n_var_lndp = 6
    integer                 :: ntasks,fid
    integer                 :: nthreads
    integer            :: ncid, xt_dim_id,yt_dim_id,time_dim_id,xt_var_id,&
                          yt_var_id,time_var_id,var_id_lat,var_id_lon,var_id_tile
    integer            :: varid1,varid2,varid3,varid4,varid_lon,varid_lat,varid_tile
    integer                 :: varidl(max_n_var_lndp)
    integer                 :: zt_dim_id,zt_var_id
    character*2             :: strid

    character(len=3), dimension(max_n_var_lndp)         ::  lndp_var_list
    real(kind=kind_dbl_prec), dimension(max_n_var_lndp) ::  lndp_prt_list
    
    real :: ak(nlevs+1),bk(nlevs+1)
    real(kind=4) :: ts,undef

    data ak(:) /0.0, 306.1489, 13687.72    , 0.99/
    data bk(:) /1.0,   0.9284,     0.013348, 0.0/
    integer     :: nb,blksz_1,nblks,ierr,my_id,i,j,k,l,nx,ny !,id
    integer     :: isc,iec,jsc,jec,isd,ied,jsd,jed
    integer :: halo_update_type = 1
    real        :: dx,dy,pi,rd,cp
    ! logical   :: write_this_tile
    integer  :: nargs,ntile_out,nlunit,pe,npes,stackmax=4000000
    integer  :: i1,i2,j1,npts,istart,tpt
    character*80 :: fname
    ! character*1  :: ntile_out_str
    integer :: comm

!4.25.23 these below are not being used
    ! real(kind=4),allocatable,dimension(:,:) :: workg,tile_number
    ! real(kind=4),allocatable,dimension(:,:,:) :: workg3d

    real(kind=4),allocatable,dimension(:) :: grid_xt,grid_yt
    real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlat
    real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlon
    real(kind=kind_dbl_prec)    :: ex3d(nlevs+1),pressi(nlevs+1),pressl(nlevs),&
                                   p1000,exn

    type(grid_box_type)           :: grid_box
    real (kind=kind_phys),allocatable :: shum_wts  (:,:,:)
    real (kind=kind_phys),allocatable :: sppt_wts  (:,:,:)
    real (kind=kind_phys),allocatable :: sppt_pattern(:,:)
    real (kind=kind_phys),allocatable :: skebu_wts (:,:,:)
    real (kind=kind_phys),allocatable :: skebv_wts (:,:,:)
    real (kind=kind_phys),allocatable :: sfc_wts   (:,:,:)
    integer,allocatable :: blksz(:)
    integer              :: me              !< MPI rank designator
    integer              :: root_pe         !< MPI rank of root atmosphere processor
    real(kind=kind_phys) :: dtp             !< physics timestep in seconds
    real(kind=kind_phys) :: fhour           !< previous forecast hour
    real(kind=kind_phys) :: sppt_amp        !< amplitude of sppt (to go to cld scheme)
    logical  :: do_sppt,do_shum,do_skeb,use_zmtnblck
    integer  ::  skeb_npass,n_var_lndp, lndp_type
    character(len=65) :: fn_nml                   !< namelist filename
    character(len=256),allocatable :: input_nml_file(:) !< character string containing full namelist
    
    ! Integer, intent(in)   :: LENS_OUT, num_ens
    ! Real, intent(out)     :: vector_rand_ens(LENS_OUT)
    Real, allocatable     :: vector_rand_ens(:)
    INTEGER :: IDIM, JDIM, NUM_TILES, IY, IM, ID, IH, &
               i_layout, j_layout, tot_subtiles, my_tile
    REAL    :: FH, DELTSFC
    ! INTEGER :: IERR
    INTEGER :: NPROCS, MYRANK, NUM_THREADS, NUM_PARTHDS, MAX_TASKS
    REAL                :: horz_len_scale, ver_len_scale, temp_len_scale 
    Integer             :: ens_size, n_surf_vars, t_indx, t_len
    ! LOGICAL             :: rcov_localize, ens_inflate
    CHARACTER(LEN=500)  :: static_filename, fv3_prefix, vector_prefix
    CHARACTER(LEN=50)   :: rand_var
    integer             :: PRINTRANK
    LOGICAL             :: print_debg_info
  
    Integer                :: LENSFC, LENSFC_landm, vector_size
    INTEGER, allocatable   :: tile_xy(:), Idim_xy(:), Jdim_xy(:), VETFCS(:)
    Real, allocatable      :: RLA_land(:), RLO_land(:), OROG_land(:)
    ! Real, allocatable      :: rand_Ta(:), rand_Prec(:), rand_sRad(:), rand_lRad(:)
    ! Real, allocatable      :: rand_sHum(:), rand_wSpeed 
    Real, allocatable      :: rand_Ta3D(:,:,:), sfc_wts_out(:,:,:)
    ! Real, allocatable      :: vector_rand_Ta(:)
    CHARACTER(LEN=500)     :: fv3_inp_file, vector_inp_file
    CHARACTER(LEN=2)       :: RANKCH 

    CHARACTER(len=250)     :: forc_inp_path, forc_inp_file
    CHARACTER(len=500)     :: forc_inp_file_ens
    CHARACTER(LEN=2)       :: peg_str
    CHARACTER(LEN=4)       :: ensCH
    character(len=3)       :: mem_str
    CHARACTER(LEN=3)       :: proc_str
    Real                   :: std_dev_f(8)
    
    Integer                :: ixy, ipr, arr_indx, iproc, ie, iv, it
    ! real,allocatable       :: DUMMY1(:,:)
    INTEGER                :: mpiReal_size, rsize, isize, mpiInt_size
    integer, allocatable   :: tile_members(:,:), tile_group(:), comm_tile(:)
    integer                :: group_world  !, comm_tile !comm_world, 
    integer                :: Np_ext, Np_til, p_gN, p_gRank
    integer, allocatable   :: pg_bindx(:,:)

    LOGICAL             :: file_exists
    Integer             :: error, varid !ncid, 
    ! character(len=50), dimension (8)   :: forc_var_list
    Real, allocatable      :: forcArray(:)
    Real(dp), allocatable  :: C_cov(:,:), L_cov(:,:)
    Integer, parameter, dimension(8)  :: ens_gen_type = (/1, 1, 0, 0, 1, 1, 1, 1/)
 !   character(len=*), parameter, dimension (8)   :: forc_var_list_par = &
 !         (/"precipitation", "solar_radiation",   &
 !           "longwave_radiation", "temperature",           &
 !           "wind_speed", "specific_humidity", "precipitation_conserve", &
 !           "surface_pressure"/) !not being perturbed atm

    character(len=128), dimension (8)   :: forc_var_list = &
          (/"precipitation", "solar_radiation",   &
            "longwave_radiation", "temperature",           &
            "wind_speed", "specific_humidity", "precipitation_conserve", &
            "surface_pressure"/) !not being perturbed atm
    
    character(len=500)      :: stoch_ini_file !11.8.21 TZG input init pattern file

    NAMELIST/NAMSNO/ IDIM, JDIM, NUM_TILES, i_layout, j_layout, IY, IM, ID, IH, FH, DELTSFC, &
                    horz_len_scale, ver_len_scale, temp_len_scale, ens_size, &
                    t_len, t_indx, &
                    static_filename, fv3_prefix, vector_prefix, rand_var, &
                    PRINTRANK, print_debg_info, n_surf_vars, &
                    vector_size, forc_inp_path, forc_inp_file, std_dev_f, forc_var_list             
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
    DATA print_debg_info/.false./
    DATA PRINTRANK/4/
    DATA static_filename/""/
    DATA fv3_prefix/"./"/
    DATA vector_prefix/""/
    DATA rand_var/"smc"/
    Data vector_size/18360/
    Data n_surf_vars/6/
    Data forc_inp_path/"./"/
    Data forc_inp_file/"C96_GDAS_forcing_2019-12-15.nc"/   
    Data std_dev_f/0.1, 0.05, 20.0, 2.0, 0.05, 0.05, 0.1, 0.01/ 
    Data forc_var_list/"precipitation", "solar_radiation",   &
            "longwave_radiation", "temperature",           &
            "wind_speed", "specific_humidity", "precipitation_conserve", &
            "surface_pressure"/ !not being perturbed atm   
    ! DATA obs_srch_rad/250.0/   
    ! DATA stdev_obsv_depth/40.0/
    ! DATA stdev_obsv_sncov/80.0/
    ! DATA stdev_back/30.0/
    ! DATA num_assim_steps/1/  ! For multiple time steps of assimilation
    ! DATA dT_Asssim/24.0/     ! hrs. For multiple time steps of assimilation
    
    namelist /gfs_physics_nml/do_sppt,do_skeb,do_shum,lndp_type,n_var_lndp
    
    nlunit=23
    open (unit=nlunit, file='input.nml', READONLY, status='OLD')
    n_var_lndp=0
    lndp_type=0
    do_sppt=.false.
    do_shum=.false.
    do_skeb=.false.
    read(nlunit,gfs_physics_nml)
    close(nlunit)
    ! define stuff
    pi=3.14159265359
    undef=9.99e+20
    p1000=100000.0
    !define mid-layer pressure
    rd=287.0
    cp=1004.0
    DO k=1,nlevs
    pressi(k)=ak(k)+p1000*bk(k)
    ENDDO
    ex3d=cp*(pressi/p1000)**(rd/cp)
    DO k=1,nlevs
    exn = (ex3d(k)*pressi(k)-ex3d(k+1)*pressi(k+1))/((cp+rd)*(pressi(k)-pressi(k+1)))
    pressl(k)=p1000*exn**(cp/rd)
    ENDDO    

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

    ! PRINT*,"READING NAMELIST."
    ! CALL BAOPENR(360, "tiles_to_vector.nml", IERR)             !"snowDA.nml"   !
    open(360, file="generate_ens_forc.nml", form="formatted")
    read(360, NAMSNO)
    ! READ(360, NML=NAMSNO)
    close(360)    
    IF (MYRANK==0) WRITE(6, NAMSNO)
    ! LENSFC = IDIM*JDIM ! TOTAL NUMBER OF POINTS FOR THE CUBED-SPHERE TILE
    
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
    ! if (myrank < 6 ) then
    !     call mpp_set_current_pelist(tile_members(1,:))
    !     print*, "pelist 1"
    !     my_id=mpp_pe()
    !     ntasks=mpp_npes()
    !     print*, "proc ", my_id, " of ", ntasks
    !     call atmosphere_init_stub (grid_box)
    ! else
    !     call mpp_set_current_pelist(tile_members(2,:))
    !     print*, "pelist 2"
    !     my_id=mpp_pe()
    !     ntasks=mpp_npes()
    !     print*, "proc ", my_id, " of ", ntasks
    !     call atmosphere_init_stub (grid_box)
    ! endif 

    call mpp_set_current_pelist()
    ! my_id=mpp_pe()
    ! ntasks=mpp_npes()
    ! print*, "global scope"
    ! print*, "proc ", my_id, " of ", ntasks
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

    ! allocate(workg(nx,ny))
    ! allocate(tile_number(nx,ny))
    ! allocate(workg3d(nx,ny,nlevs))

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
    
    allocate(grid_xt(nx),grid_yt(ny))
    do i=1,nx
    grid_xt(i)=i
    enddo
    do j=1,ny
    grid_yt(j)=j
    enddo
    !     endif
    ! enddo

    ! forc_var_list(1) = "temperature" 
    ! forc_var_list(2) = "precipitation_bilinear"
    ! forc_var_list(3) = "precipitation_conserve"
    ! forc_var_list(4) = "solar_radiation"
    ! forc_var_list(5) = "longwave_radiation"
    ! forc_var_list(6) = "specific_humidity"
    ! forc_var_list(7) = "wind_speed"
    ! forc_var_list(8) = "surface_pressure"     

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

    ! allocate(vector_rand_Ta(vector_size))

    ! allocate(rand_Ta(LENSFC), rand_Prec(LENSFC), rand_sRad(LENSFC))
    ! allocate(rand_lRad(LENSFC),rand_sHum(LENSFC), rand_wSpeed(LENSFC))
    allocate(rand_Ta3D(NUM_TILES, JDIM, IDIM))
    
    ! Read fv3    
    !workg_T162_984x488.tile03.nc 
    allocate(sfc_wts_out(JDIM / j_layout, IDIM / i_layout, n_surf_vars))
    allocate(forcArray(vector_size))  !, t_len))
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
    do_sppt=.false.
    do_shum=.false.
    do_skeb=.false.
    ! call init_stochastic_physics(nlevs, blksz, dtp, sppt_amp,                         &
    !         input_nml_file, fn_nml, nlunit, xlon, xlat, do_sppt, do_shum,                &
    !         do_skeb, lndp_type, n_var_lndp, use_zmtnblck, skeb_npass, &
    !         lndp_var_list, lndp_prt_list,    &
    !         ak, bk, nthreads, root_pe, comm, ierr)
    ! if (ierr .ne. 0) print *, 'ERROR init_stochastic_physics call' ! Draper - need proper error trapping here
    if (do_sppt) allocate(sppt_wts(nblks,blksz_1,nlevs))
    if (do_shum) allocate(shum_wts(nblks,blksz_1,nlevs))
    if (do_skeb) allocate(skebu_wts(nblks,blksz_1,nlevs))
    if (do_skeb) allocate(skebv_wts(nblks,blksz_1,nlevs))
    if (lndp_type > 0) allocate(sfc_wts(nblks,blksz_1,n_var_lndp))
    ! allocate(sppt_wts(nblks,blksz_1,nlevs))
    ! allocate(shum_wts(nblks,blksz_1,nlevs))
    ! allocate(skebu_wts(nblks,blksz_1,nlevs))
    ! allocate(skebv_wts(nblks,blksz_1,nlevs))
    ! allocate(sfc_wts(nblks,blksz_1,n_var_lndp))
    
    comm = comm_tile(p_gN+1)    !MPI_COMM_WORLD
    root_pe = 0  !if (p_gRank == 0)  root_pe = myrank  ! mpp_root_pe()   
    
    ! Np_til = NPROCS / NUM_TILES  ! num tile groups    
    ! p_gN = MYRANK / NUM_TILES  ! group for proc.  
    ! p_gRank = MOD(MYRANK, NUM_TILES)  ! proc. rank within group

    ! WRITE(proc_str, '(I0)') myrank
    Do ie = p_gN*(ens_size/Np_til)+1, (p_gN+1)*(ens_size/Np_til)  !ens_size
        WRITE(ensCH, '(I0)') ie
        write(mem_str, '(I3.3)') ie
        stoch_ini_file = TRIM(forc_inp_path)//'/RESTART/stochy_final_ens'//TRIM(ensCH)//'.nc'
        if (p_gRank == 0) then  
            forc_inp_file_ens=TRIM(forc_inp_path)//"/mem"//TRIM(mem_str)//"/"//TRIM(forc_inp_file)
            ! forcing file
            INQUIRE(FILE=trim(forc_inp_file_ens), EXIST=file_exists)
            if (.not. file_exists) then 
                print *, 'error,file does not exist', &   
                        trim(forc_inp_file_ens) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10, IERR) ! CSD - add proper error trapping?
            endif
            error = nf90_open(trim(forc_inp_file_ens), NF90_WRITE, ncid)
            call netcdf_err(error, 'opening forcing file' )
        endif
        ! IF (p_gRank==0) PRINT*," calling stand alone stochy, ens ", ie
        do_sppt=.false.
        do_shum=.false.
        do_skeb=.false.
        ! call standalone_stochy_sfc(IDIM, JDIM, n_surf_vars, sfc_wts_out)        
        call init_stochastic_physics(nlevs, blksz, dtp, sppt_amp,                         &
            input_nml_file, stoch_ini_file, fn_nml, nlunit, xlon, xlat, do_sppt, do_shum,                &
            do_skeb, lndp_type, n_var_lndp, use_zmtnblck, skeb_npass, &
            lndp_var_list, lndp_prt_list,    &
            ak, bk, nthreads, root_pe, comm_tile(p_gN+1), ierr)
        if (ierr .ne. 0) print *, 'ERROR init_stochastic_physics call' ! Draper - need proper error trapping here
        do_sppt=.false.
        do_shum=.false.
        do_skeb=.false.
        ! if (p_gRank == 0) then
        !     print*, " pe group", p_gN + 1
        !     PRINT*," done init_stochastic_physics, ens ", ie
        ! endif
        
        Do it=1, t_len
            call run_stochastic_physics(nlevs, it-1, fhour, blksz, &
                        sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts, &
                                        skebv_wts=skebv_wts, sfc_wts=sfc_wts, &
                                        nthreads=nthreads)
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
            ! if (p_gRank == 0) then
            !     print*, " pe group", p_gN + 1
            !     PRINT*," Done run_stochastic_physics, ens", ie, "time", it    
            ! endif

            Do ixy = 1, n_surf_vars        !vector_length 
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

                    ! Start writing restart file
                    error = nf90_inq_varid(ncid, trim(forc_var_list(ixy)), varid)
                    call netcdf_err(error, 'getting varid '//trim(forc_var_list(ixy)) )
                    !read
                    ERROR=NF90_GET_VAR(NCID, varid, forcArray, &
                            start = (/1,it/), count = (/vector_size, 1/))
                    CALL NETCDF_ERR(ERROR, 'ERROR READING '//trim(forc_var_list(ixy)) )                  
                    ! print*, trim(forc_var_list(ixy))
                    ! print*, forcArray
                    ! print*, "rand pattern"
                    ! print*, vector_rand_ens
                    if (ens_gen_type(ixy) == 1) then 
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
                    error = nf90_put_var(ncid, varid , forcArray,       &
                        start = (/1,it/), count = (/vector_size, 1/))
                    call netcdf_err(error, 'writing '//trim(forc_var_list(ixy)) )
                    ! if(ixy == 1) then
                    !     ! precip type 2
                    !     error = nf90_inq_varid(ncid, trim(forc_var_list(7)), varid)
                    !     call netcdf_err(error, 'getting varid '//trim(forc_var_list(7)) )
                    !     !read
                    !     ERROR=NF90_GET_VAR(NCID, varid, forcArray, &
                    !             start = (/1,it/), count = (/vector_size, 1/))
                    !     CALL NETCDF_ERR(ERROR, 'ERROR READING '//trim(forc_var_list(7)) )                  
                    !     ! print*, trim(forc_var_list(7))
                    !     ! print*, forcArray
                    !     forcArray = vector_rand_ens * forcArray 
                    !     ! print*, trim(forc_var_list(ixy)), " after mult"
                    !     ! print*, forcArray
                    !     error = nf90_put_var(ncid, varid , forcArray,       &
                    !         start = (/1,it/), count = (/vector_size, 1/))
                    !     call netcdf_err(error, 'writing '//trim(forc_var_list(7)) ) 
                    ! endif                                      
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
    if (do_sppt) deallocate(sppt_wts)
    if (do_shum) deallocate(shum_wts)
    if (do_skeb) deallocate(skebu_wts)
    if (do_skeb) deallocate(skebv_wts)
    if (lndp_type > 0) deallocate(sfc_wts)
    
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
    deallocate(sfc_wts_out)
    deallocate(forcArray)
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

 END Program GenEnsForc

!>@brief calls 'stochastic_physics' for the 
! initialization of the stochastic physics random pattern generators
! and running for one time step
! initiadapted from standalone_stochy by Phil Pigion 

!  subroutine  standalone_stochy_sfc(IDIM, JDIM, n_surf_vars, sfc_wts_out)
  
!     use stochastic_physics,  only : init_stochastic_physics,run_stochastic_physics
!     use get_stochy_pattern_mod,  only : write_stoch_restart_atm

!     use atmosphere_stub_mod, only: Atm,atmosphere_init_stub
!     !use mpp_domains
!     use mpp_mod,             only: mpp_set_current_pelist,mpp_get_current_pelist,mpp_init,mpp_pe,mpp_npes ,mpp_declare_pelist,mpp_root_pe
!     use mpp_domains_mod,     only: mpp_broadcast_domain,MPP_DOMAIN_TIME,mpp_domains_init ,mpp_domains_set_stack_size
!     use fms_mod,             only:  fms_init
!     use xgrid_mod,           only: grid_box_type
!     use netcdf
!     use kinddef,             only : kind_dbl_prec,kind_phys
!     use stochy_namelist_def, only : stochini

!     implicit none
!     integer, parameter      :: nlevs=3
!     integer, parameter :: max_n_var_lndp = 6
!     integer                 :: ntasks,fid
!     integer                 :: nthreads
!     integer                 :: ncid,xt_dim_id,yt_dim_id,time_dim_id,xt_var_id,yt_var_id,time_var_id,var_id_lat,var_id_lon,var_id_tile
!     integer                 :: varid1,varid2,varid3,varid4,varid_lon,varid_lat,varid_tile
!     integer                 :: varidl(max_n_var_lndp)
!     integer                 :: zt_dim_id,zt_var_id
!     character*2             :: strid

!     character(len=3), dimension(max_n_var_lndp)         ::  lndp_var_list
!     real(kind=kind_dbl_prec), dimension(max_n_var_lndp) ::  lndp_prt_list
!     include 'mpif.h'
!     include 'netcdf.inc'
!     real :: ak(nlevs+1),bk(nlevs+1)
!     real(kind=4) :: ts,undef

!     data ak(:) /0.0, 306.1489, 13687.72    , 0.99/
!     data bk(:) /1.0,   0.9284,     0.013348, 0.0/
!     integer     :: nb,blksz_1,nblks,ierr,my_id,i,j,k,l,nx,ny,id
!     integer     :: isc,iec,jsc,jec,isd,ied,jsd,jed
!     integer :: halo_update_type = 1
!     real        :: dx,dy,pi,rd,cp
!     logical   :: write_this_tile
!     integer  :: nargs,ntile_out,nlunit,pe,npes,stackmax=4000000
!     integer  :: i1,i2,j1,npts,istart,tpt
!     character*80 :: fname
!     character*1  :: ntile_out_str
!     integer :: comm

!     real(kind=4),allocatable,dimension(:,:) :: workg,tile_number
!     real(kind=4),allocatable,dimension(:,:,:) :: workg3d
!     real(kind=4),allocatable,dimension(:) :: grid_xt,grid_yt
!     real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlat
!     real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlon
!     real(kind=kind_dbl_prec)                    :: ex3d(nlevs+1),pressi(nlevs+1),pressl(nlevs),p1000,exn

!     type(grid_box_type)           :: grid_box
!     real (kind=kind_phys),allocatable :: shum_wts  (:,:,:)
!     real (kind=kind_phys),allocatable :: sppt_wts  (:,:,:)
!     real (kind=kind_phys),allocatable :: sppt_pattern(:,:)
!     real (kind=kind_phys),allocatable :: skebu_wts (:,:,:)
!     real (kind=kind_phys),allocatable :: skebv_wts (:,:,:)
!     real (kind=kind_phys),allocatable :: sfc_wts   (:,:,:)
!     integer,allocatable :: blksz(:)
!     integer              :: me              !< MPI rank designator
!     integer              :: root_pe         !< MPI rank of root atmosphere processor
!     real(kind=kind_phys) :: dtp             !< physics timestep in seconds
!     real(kind=kind_phys) :: fhour           !< previous forecast hour
!     real(kind=kind_phys) :: sppt_amp        !< amplitude of sppt (to go to cld scheme)
!     logical  :: do_sppt,do_shum,do_skeb,use_zmtnblck
!     integer  ::  skeb_npass,n_var_lndp, lndp_type
!     character(len=65) :: fn_nml                   !< namelist filename
!     character(len=256),allocatable :: input_nml_file(:) !< character string containing full namelist
    
!     !10.14.21
!     integer          :: IDIM, JDIM, n_surf_vars
!     real             :: sfc_wts_out(n_surf_vars, IDIM, JDIM)

!         namelist /gfs_physics_nml/do_sppt,do_skeb,do_shum,lndp_type,n_var_lndp
!     write_this_tile=.false.
!     ntile_out_str='0'
!     nlunit=23
!     nargs=iargc()
!     if (nargs.EQ.1) then
!     call getarg(1,ntile_out_str)
!     endif
!     read(ntile_out_str,'(I1.1)') ntile_out
!     open (unit=nlunit, file='input.nml', READONLY, status='OLD')
!     n_var_lndp=0
!     lndp_type=0
!     do_sppt=.false.
!     do_shum=.false.
!     do_skeb=.false.
!     read(nlunit,gfs_physics_nml)
!     close(nlunit)
!     ! define stuff
!     pi=3.14159265359
!     undef=9.99e+20
!     p1000=100000.0
!     !define mid-layer pressure
!     rd=287.0
!     cp=1004.0
!     DO k=1,nlevs
!     pressi(k)=ak(k)+p1000*bk(k)
!     ENDDO
!     ex3d=cp*(pressi/p1000)**(rd/cp)
!     DO k=1,nlevs
!     exn = (ex3d(k)*pressi(k)-ex3d(k+1)*pressi(k+1))/((cp+rd)*(pressi(k)-pressi(k+1)))
!     pressl(k)=p1000*exn**(cp/rd)
!     ENDDO
    
!     ! print*, "starting stanalone stochy"
!     ! call fms_init()
!     ! print*, "starting stanalone stochy 2"
!     ! call mpp_init()
!     ! print*, "starting stanalone stochy 3"
!     ! call fms_init
!     ! print*, "starting stanalone stochy 4"
!     my_id=mpp_pe()
!     ntasks=mpp_npes()
!     do_sppt=.false.
!     do_shum=.false.
!     do_skeb=.false.
!     print*, "proc ", my_id, " of ", ntasks, " doing only sfc rand"
    
!     call atmosphere_init_stub (grid_box)
!     isd=Atm(1)%bd%isd
!     ied=Atm(1)%bd%ied
!     jsd=Atm(1)%bd%jsd
!     jed=Atm(1)%bd%jed
!     isc=Atm(1)%bd%isc
!     iec=Atm(1)%bd%iec
!     jsc=Atm(1)%bd%jsc
!     jec=Atm(1)%bd%jec
!     nx=iec-isc+1
!     ny=jec-jsc+1
!     allocate(workg(nx,ny))
!     allocate(tile_number(nx,ny))
!     allocate(workg3d(nx,ny,nlevs))
!     print*,'nx,ny=',nx,ny
!     blksz_1=nx
!     nblks=nx*ny/blksz_1
!     allocate(blksz(nblks))
!     do i=1,nblks
!       blksz(i)=blksz_1
!     enddo
!     nthreads = 1
!     me=my_id
!     fhour=0
!     dtp=900
!     fn_nml='input.nml'
!     nlunit=21
    
!     !define model grid
!     dx=360.0/nx
!     dy=180.0/ny
!     allocate(xlat(nblks,blksz_1))
!     allocate(xlon(nblks,blksz_1))
!     i1=isc
!     j1=jsc
!     do nb=1,nblks
!         i2=i1+blksz_1-1
!         if (i2 .le. iec) then 
!            xlon(nb,1:blksz_1) = Atm(1)%gridstruct%agrid_64(i1:i2,j1,1)
!            xlat(nb,1:blksz_1) = Atm(1)%gridstruct%agrid_64(i1:i2,j1,2)
!            i1=i1+blksz_1
!         else
!            npts=iec-i1+1
!            xlon(nb,1:npts) = Atm(1)%gridstruct%agrid_64(i1:iec,j1,1)
!            xlat(nb,1:npts) = Atm(1)%gridstruct%agrid_64(i1:iec,j1,2)
!            if (j1.LT. jec) then
!               xlon(nb,npts+1:blksz_1) = Atm(1)%gridstruct%agrid_64(isc:isc+(blksz_1-npts+1),j1+1,1)
!               xlat(nb,npts+1:blksz_1) = Atm(1)%gridstruct%agrid_64(isc:isc+(blksz_1-npts+1),j1+1,2)
!            endif
!            i1=npts+1
!            j1=j1+1
!         endif
!         if (i2.EQ.iec) then
!            i1=isc
!            j1=j1+1
!         endif
!     end do
    
!     allocate(grid_xt(nx),grid_yt(ny))
!     do i=1,nx
!       grid_xt(i)=i
!     enddo
!     do j=1,ny
!       grid_yt(j)=j
!     enddo
!     print*,'calling init_stochastic_physics',nlevs
!     root_pe=mpp_root_pe()
!     allocate(input_nml_file(1))
!     input_nml_file='input.nml'
!     comm=MPI_COMM_WORLD
!     call init_stochastic_physics(nlevs, blksz, dtp, sppt_amp,                         &
!          input_nml_file, fn_nml, nlunit, xlon, xlat, do_sppt, do_shum,                &
!          do_skeb, lndp_type, n_var_lndp, use_zmtnblck, skeb_npass, &
!          lndp_var_list, lndp_prt_list,    &
!          ak, bk, nthreads, root_pe, comm, ierr)
!     if (ierr .ne. 0) print *, 'ERROR init_stochastic_physics call' ! Draper - need proper error trapping here
 
!     if (do_sppt)allocate(sppt_wts(nblks,blksz_1,nlevs))
!     if (do_shum)allocate(shum_wts(nblks,blksz_1,nlevs))
!     if (do_skeb)allocate(skebu_wts(nblks,blksz_1,nlevs))
!     if (do_skeb)allocate(skebv_wts(nblks,blksz_1,nlevs))
!     if (lndp_type > 0) then 
!         allocate(sfc_wts(nblks,blksz_1,n_var_lndp))
!         ! allocate(sfc_wts_out(nblks,blksz_1,n_var_lndp))
!     endif

!     ! if (stochini) then
!     !    istart=11
!     ! else
!     !    istart=1
!     ! endif
!     call run_stochastic_physics(nlevs, 0, fhour, blksz, &
!                                 sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts, skebv_wts=skebv_wts, sfc_wts=sfc_wts, &
!                                 nthreads=nthreads)
!     do l=1,n_var_lndp
!         do j=1,ny
!             sfc_wts_out(l,:,j)=sfc_wts(j,:,l)
!         enddo
!     enddo

! end subroutine standalone_stochy_sfc
!  subroutine write_rand_pattern(file_name)
   
!     ! call get_outfile(fname)
!     ! write(strid,'(I2.2)') my_id+1
!     ! if (ntile_out.EQ.0) write_this_tile=.true.
!     ! if ((my_id+1).EQ.ntile_out) write_this_tile=.true.
!     ! print*,trim(fname)//'.tile'//strid//'.nc',write_this_tile
!     CHARACTER(LEN=500)  :: rand_pattern_filename
!     if (write_this_tile) then
        
!         rand_pattern_filename = trim(fname)//'.tile'//strid//'.nc'
 
!         fid=30+my_id
!         ierr=nf90_create(trim(fname)//'.tile'//strid//'.nc',cmode=NF90_CLOBBER,ncid=ncid)
!         ierr=NF90_DEF_DIM(ncid,"grid_xt",nx,xt_dim_id)
!         ierr=NF90_DEF_DIM(ncid,"grid_yt",ny,yt_dim_id)
!         if (do_skeb)ierr=NF90_DEF_DIM(ncid,"p_ref",nlevs,zt_dim_id)
!         ierr=NF90_DEF_DIM(ncid,"time",NF90_UNLIMITED,time_dim_id)
!        !> - Define the dimension variables.
!         ierr=NF90_DEF_VAR(ncid,"grid_xt",NF90_FLOAT,(/ xt_dim_id /), xt_var_id)
!         ierr=NF90_PUT_ATT(ncid,xt_var_id,"long_name","T-cell longitude")
!         ierr=NF90_PUT_ATT(ncid,xt_var_id,"cartesian_axis","X")
!         ierr=NF90_PUT_ATT(ncid,xt_var_id,"units","degrees_E")
!         ierr=NF90_DEF_VAR(ncid,"grid_yt",NF90_FLOAT,(/ yt_dim_id /), yt_var_id)
!         ierr=NF90_PUT_ATT(ncid,yt_var_id,"long_name","T-cell latitude")
!         ierr=NF90_PUT_ATT(ncid,yt_var_id,"cartesian_axis","Y")
!         ierr=NF90_PUT_ATT(ncid,yt_var_id,"units","degrees_N")
!         ierr=NF90_DEF_VAR(ncid,"grid_lat",NF90_FLOAT,(/ xt_dim_id, yt_dim_id, time_dim_id /), var_id_lat)
!         ierr=NF90_PUT_ATT(ncid,var_id_lat,"long_name","T-cell latitudes")
!         ierr=NF90_PUT_ATT(ncid,var_id_lat,"units","degrees_N")
!         ierr=NF90_PUT_ATT(ncid,var_id_lat,"missing_value",undef)
!         ierr=NF90_PUT_ATT(ncid,var_id_lat,"_FillValue",undef)
!         ierr=NF90_DEF_VAR(ncid,"grid_lon",NF90_FLOAT,(/ xt_dim_id, yt_dim_id, time_dim_id /), var_id_lon)
!         ierr=NF90_PUT_ATT(ncid,var_id_lon,"long_name","T-cell longitudes")
!         ierr=NF90_PUT_ATT(ncid,var_id_lon,"units","degrees_N")
!         ierr=NF90_PUT_ATT(ncid,var_id_lon,"missing_value",undef)
!         ierr=NF90_PUT_ATT(ncid,var_id_lon,"_FillValue",undef)
!         ierr=NF90_DEF_VAR(ncid,"tile_num",NF90_FLOAT,(/ xt_dim_id, yt_dim_id, time_dim_id /), var_id_tile)
!         ierr=NF90_PUT_ATT(ncid,var_id_tile,"long_name","tile number")
!         ierr=NF90_PUT_ATT(ncid,var_id_tile,"missing_value",undef)
!         ierr=NF90_PUT_ATT(ncid,var_id_tile,"_FillValue",undef)
!         if (do_skeb)then
!            ierr=NF90_DEF_VAR(ncid,"p_ref",NF90_FLOAT,(/ zt_dim_id /), zt_var_id)
!            ierr=NF90_PUT_ATT(ncid,zt_var_id,"long_name","reference pressure")
!            ierr=NF90_PUT_ATT(ncid,zt_var_id,"cartesian_axis","Z")
!            ierr=NF90_PUT_ATT(ncid,zt_var_id,"units","Pa")
!         endif
!         ierr=NF90_DEF_VAR(ncid,"time",NF90_FLOAT,(/ time_dim_id /), time_var_id)
!         ierr=NF90_DEF_VAR(ncid,"time",NF90_FLOAT,(/ time_dim_id /), time_var_id)
!         ierr=NF90_PUT_ATT(ncid,time_var_id,"long_name","time")
!         ierr=NF90_PUT_ATT(ncid,time_var_id,"units","hours since 2014-08-01 00:00:00")
!         ierr=NF90_PUT_ATT(ncid,time_var_id,"cartesian_axis","T")
!         ierr=NF90_PUT_ATT(ncid,time_var_id,"calendar_type","JULIAN")
!         ierr=NF90_PUT_ATT(ncid,time_var_id,"calendar","JULIAN")
!         if (do_sppt)then
!            ierr=NF90_DEF_VAR(ncid,"sppt_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), varid1)
!            ierr=NF90_PUT_ATT(ncid,varid1,"long_name","sppt pattern")
!            ierr=NF90_PUT_ATT(ncid,varid1,"units","None")
!            ierr=NF90_PUT_ATT(ncid,varid1,"missing_value",undef)
!            ierr=NF90_PUT_ATT(ncid,varid1,"_FillValue",undef)
!            ierr=NF90_PUT_ATT(ncid,varid1,"cell_methods","time: point")
!         endif
!         if (do_shum)then
!            ierr=NF90_DEF_VAR(ncid,"shum_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), varid2)
!            ierr=NF90_PUT_ATT(ncid,varid2,"long_name","shum pattern")
!            ierr=NF90_PUT_ATT(ncid,varid2,"units","None")
!            ierr=NF90_PUT_ATT(ncid,varid2,"missing_value",undef)
!            ierr=NF90_PUT_ATT(ncid,varid2,"_FillValue",undef)
!            ierr=NF90_PUT_ATT(ncid,varid2,"cell_methods","time: point")
!         endif
!         if (do_skeb)then
!            ierr=NF90_DEF_VAR(ncid,"skebu_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,zt_dim_id,time_dim_id/), varid3)
!            ierr=NF90_DEF_VAR(ncid,"skebv_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,zt_dim_id,time_dim_id/), varid4)
!            ierr=NF90_PUT_ATT(ncid,varid3,"long_name","skeb u pattern")
!            ierr=NF90_PUT_ATT(ncid,varid3,"units","None")
!            ierr=NF90_PUT_ATT(ncid,varid3,"missing_value",undef)
!            ierr=NF90_PUT_ATT(ncid,varid3,"_FillValue",undef)
!            ierr=NF90_PUT_ATT(ncid,varid3,"cell_methods","time: point")
!            ierr=NF90_PUT_ATT(ncid,varid4,"long_name","skeb v pattern")
!            ierr=NF90_PUT_ATT(ncid,varid4,"units","None")
!            ierr=NF90_PUT_ATT(ncid,varid4,"missing_value",undef)
!            ierr=NF90_PUT_ATT(ncid,varid4,"_FillValue",undef)
!            ierr=NF90_PUT_ATT(ncid,varid4,"cell_methods","time: point")
!         endif
!         if (lndp_type > 0)then
!            do l=1,n_var_lndp
!               ierr=NF90_DEF_VAR(ncid,lndp_var_list(l),NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), varidl(l))
!               ierr=NF90_PUT_ATT(ncid,varidl(l),"long_name",lndp_var_list(l)//" pattern")
!               ierr=NF90_PUT_ATT(ncid,varidl(l),"units","None")
!               ierr=NF90_PUT_ATT(ncid,varidl(l),"missing_value",undef)
!               ierr=NF90_PUT_ATT(ncid,varidl(l),"_FillValue",undef)
!               ierr=NF90_PUT_ATT(ncid,varidl(l),"cell_methods","time: point")
!            enddo
!         endif
!         ierr=NF90_ENDDEF(ncid)
!         ierr=NF90_PUT_VAR(ncid,xt_var_id,grid_xt)
!         ierr=NF90_PUT_VAR(ncid,yt_var_id,grid_yt)

!          ! put lat lon and tile number
!     !ierr=NF90_PUT_VAR(ncid,var_id_lon,transpose(xlon(isc:iec,jsc:iec)),(/1,1,1/))
!     !ierr=NF90_PUT_VAR(ncid,var_id_lat,transpose(xlat(isc:iec,jsc:iec)),(/1,1,1/))
!     ierr=NF90_PUT_VAR(ncid,var_id_lon,transpose(xlon(:,:)),(/1,1,1/))
!     ierr=NF90_PUT_VAR(ncid,var_id_lat,transpose(xlat(:,:)),(/1,1,1/))
!     tile_number=my_id+1
       
!     ierr=NF90_PUT_VAR(ncid,var_id_tile,tile_number,(/1,1,1/))

!            do l=1,n_var_lndp
!               do j=1,ny
!                  workg(:,j)=sfc_wts(j,:,l)
!               enddo
!               ierr=NF90_PUT_VAR(ncid,varidl(l),workg,(/1,1,tpt/))
!            enddo
!         endif
!         ierr=NF90_PUT_VAR(ncid,time_var_id,ts,(/tpt/))
!         endif
!         tpt=tpt+1
!      endif
!   enddo
!   if (write_this_tile) ierr=NF90_CLOSE(ncid)


!  end subroutine write_rand_pattern

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
