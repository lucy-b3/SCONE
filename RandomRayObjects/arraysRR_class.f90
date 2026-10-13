module arraysRR_class

  use numPrecision
  use universalVariables
  use constantsRR
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile
  use rng_class,                      only : RNG
  use hashFunctions_func,             only : FNV_1

  ! Data
  use baseMgNeutronDatabase_class,    only : baseMgNeutronDatabase
  use dataRR_class,                   only : dataRR

  ! Geometry
  use coord_class,                    only : coordList
  use geometryStd_class,              only : geometryStd

  use field_inter,                    only : field
  use cartesianField_class,           only : cartesianField
  use pieceConstantField_inter,       only : pieceConstantField, pieceConstantField_CptrCast
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_geomIdx  => geomIdx, &
                                             gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr

  ! Visualisation
  use tallyMap_inter,                 only : tallyMap
  use visualiser_class,               only : visualiser
  use particle_class,                 only : particle, particleState
  use tallyAdmin_class,               only : tallyAdmin
  
  ! For locks
  use omp_lib

  implicit none
  private
  
  
  !!
  !! Object to store all arrays in random ray
  !! By default, will have the current flux, previous flux, accumulated flux,
  !! source arrays, and geometric arrays.
  !!
  !! Can be extended to having corresponding flux and source moments, as well
  !! as the fixed source and additional geometric info.
  !!
  !! Private Members
  !!   nG             -> Number of energy groups, kept for convenience.
  !!   nCells         -> Number of unique cells in the geometry, kept for convenience.
  !!   lengthPerIt    -> RR active length per iteration, kept for convenience
  !!   XSData         -> Pointer to nuclear data, for convenience.
  !!   geom           -> Pointer to geometry, for convenience.
  !!   rho            -> Stabilisation factor: 0 is no stabilisation, 1 is aggressive stabilisation
  !!   ani            -> Order of anisotropic flux moments to be stored
  !!   simulationType -> Identifies which simulation to perform: flat/linear, isotropic/anisotropic
  !!
  !!   scalarFlux     -> Array of scalar flux values of length [nG * nCells]
  !!   prevFlux       -> Array of previous scalar flux values of length [nG * nCells]
  !!   fluxScores     -> Array of scalar flux values and squared values to be reported 
  !!                     in results [nG * nCells, 2]
  !!
  !!   source      -> Array of sources [nG * nCells]
  !!   fixedSource -> Array of fixed sources [nG * nCell]
  !!   sourceIdx   -> Array of material indices containing fixed sources
  !!
  !!   volumeTracks  -> Array of sum of track lengths for computing volumes [nCells]
  !!   lengthSquared -> Array of sum of lengths^2 for use in cells with zero XS [nCells]
  !!   volume        -> Array of dimensionless cell volumes [nCells]
  !!
  !!   cellHit   -> Array of ints whether a cell was visited this iteration [nCells]
  !!   cellFound -> Array of logicals whether a cell was ever visited [nCells]
  !!   cellPos   -> Array of cell positions [3 * nCells]
  !!
  !!   locks -> Array of OpenMP locks for each geometric cell
  !!
  type, public :: arraysRR
    private
    ! Components
    class(geometryStd), pointer           :: geom        => null()
    type(dataRR)                          :: XSData      
    integer(shortInt)                     :: nG          = 0
    integer(shortInt)                     :: nT          = 0
    integer(shortInt)                     :: nCells      = 0
    real(defReal)                         :: lengthPerIt = ZERO
    real(defFlt)                          :: rho         = 0.0_defFlt
    integer(shortInt)                     :: ani         = 0
    integer(shortInt)                     :: simulationType = 0
    
    logical(defBool)                      :: feedback     = .false. 
    class(pieceConstantField), pointer        :: tempField    => null()

    integer(shortInt)           :: volPolicy   = hybrid !simAverage
    integer(shortInt)           :: missPolicy  = hybrid !srcPolicy
    logical(defBool)            :: set2D       = .false. ! Stabilises LS in 2D problems

    real(defReal)               :: totalVolume = ONE
    
    ! Flux arrays
    real(defFlt), dimension(:), allocatable    :: scalarFlux
    real(defFlt), dimension(:), allocatable    :: prevFlux
    real(defReal), dimension(:,:), allocatable :: fluxScores
   
    real(defFlt), dimension(:), allocatable    :: power
    real(defReal), dimension(:,:), allocatable :: powerScores

    ! Source arrays
    real(defFlt), dimension(:), allocatable      :: source
    real(defFlt), dimension(:), allocatable      :: fixedSource
    integer(shortInt), dimension(:), allocatable :: sourceIdx

    ! Geometry arrays
    real(defReal), dimension(:), allocatable     :: volumeTracks
    real(defReal), dimension(:), allocatable     :: allVolumeTracks
    real(defReal), dimension(:), allocatable     :: lengthSquared
    real(defReal), dimension(:), allocatable     :: volume
    integer(shortInt), dimension(:), allocatable :: cellHit
    integer(longInt), dimension(:), allocatable  :: cellTotalHit
    logical(defBool), dimension(:), allocatable  :: cellFound
    real(defReal), dimension(:,:), allocatable   :: cellPos

    ! Linear source arrays
    real(defFlt), dimension(:), allocatable   :: scalarX
    real(defFlt), dimension(:), allocatable   :: scalarY
    real(defFlt), dimension(:), allocatable   :: scalarZ
    real(defFlt), dimension(:), allocatable   :: prevX
    real(defFlt), dimension(:), allocatable   :: prevY
    real(defFlt), dimension(:), allocatable   :: prevZ
    real(defFlt), dimension(:), allocatable    :: sourceX
    real(defFlt), dimension(:), allocatable    :: sourceY
    real(defFlt), dimension(:), allocatable    :: sourceZ
    real(defFlt), dimension(:), allocatable    :: fixedX
    real(defFlt), dimension(:), allocatable    :: fixedY
    real(defFlt), dimension(:), allocatable    :: fixedZ
    real(defReal), dimension(:), allocatable   :: momMat
    real(defReal), dimension(:), allocatable   :: momTracks
    real(defReal), dimension(:), allocatable   :: centroid
    real(defReal), dimension(:), allocatable   :: centroidTracks

    real(defReal), dimension(:,:), allocatable :: xScores
    real(defReal), dimension(:,:), allocatable :: yScores
    real(defReal), dimension(:,:), allocatable :: zScores

    ! OMP locks
    integer(kind=omp_lock_kind), dimension(:), allocatable :: locks
    
    ! Other data
    real(defReal)     :: averageHit = ZERO
    integer(shortInt) :: iterations = 0

  contains
    
    ! Public procedures
    procedure :: init
    procedure :: kill

    ! Access procedures
    procedure :: getDataPointer
    procedure :: getGeomPointer
    procedure :: getTempPointer
    procedure :: getSourcePointer
    procedure :: getFluxPointer
    procedure :: getSource
    procedure :: getFixedSource
    procedure :: getPrevFlux
    procedure :: getFluxScore
    procedure :: getFluxSD
    procedure :: getPowerScore
    procedure :: getPowerSD
    procedure :: getNG
    procedure :: getVolume
    procedure :: getCellPos
    procedure :: hasHit
    procedure :: getCellHitRate
    procedure :: getSimulationType
    procedure :: found
    procedure :: wasFound
    procedure :: wasHit

    procedure :: hasFixedSource
    procedure :: getFluxAtAPoint

    procedure :: getFluxXYZPointers
    procedure :: getSourceXYZPointers
    procedure :: getCentroid
    procedure :: getMomentMatrix
    procedure :: getFluxMoments
    procedure :: getFluxMomentSDs

    ! Change individual elements of the type
    ! Predominantly for use in the transport sweep
    procedure :: incrementVolume
    procedure :: incrementLengthSquared
    procedure :: incrementCentroid
    procedure :: incrementMoments
    procedure :: hitCell
    procedure :: getAverageHitRate
    procedure :: wipeCellHits
    procedure :: newFound
    procedure :: setLock
    procedure :: unsetLock

    ! Basic RR procedures
    procedure :: resetFluxes
    procedure :: normaliseFluxAndVolume
    procedure :: updateSource
    procedure :: accumulateFluxScores
    procedure :: finaliseFluxScores
    procedure :: calculateKeff
    procedure :: calculatePower
    procedure :: zeroPrevFlux

    ! Tally results for use with MC tally machinery
    procedure :: tallyResults    

    ! Output procedures
    procedure :: outputMap
    procedure :: outputToVTK

    ! Private procedures
    procedure, private :: initialiseFixedSource
    
    procedure, private :: resetFluxesFlatIso
    procedure, private :: resetFluxesLinIso
    procedure, private :: resetFluxesLIFA
    procedure, private :: resetFluxesFlatAni
    
    procedure, private :: normaliseFluxAndVolumeFlatIso
    procedure, private :: normaliseFluxAndVolumeLinIso
    procedure, private :: normaliseFluxAndVolumeLIFA
    procedure, private :: normaliseFluxAndVolumeFlatAni
    
    procedure, private :: sourceUpdateKernelFlatIso
    procedure, private :: sourceUpdateKernelLinIso
    procedure, private :: sourceUpdateKernelLIFA
    procedure, private :: sourceUpdateKernelFlatAni

    procedure, private :: accumulateFluxScoresFlat
    procedure, private :: accumulateFluxScoresLinear

    procedure, private :: finaliseFluxScoresFlat
    procedure, private :: finaliseFluxScoresLinear
    
    procedure, private :: calculateKeffKernel
    procedure, private :: invertMatrix

  end type arraysRR

contains

  !!
  !! Initialise the arrays object
  !!
  !! The object is fed sizes and requirements by the physics package.
  !! This will allocate the necessary arrays
  !!
  subroutine init(self, db, geom, lengthPerIt, rho, lin, ani, doKinetics, loud, mapTemp, dictFS, volPolicy, missPolicy, set2D)
    class(arraysRR), intent(inout)                      :: self
    class(baseMgNeutronDatabase), pointer, intent(in)   :: db
    class(geometryStd), pointer, intent(in)             :: geom
    real(defReal), intent(in)                           :: lengthPerIt
    real(defReal), intent(in)                           :: rho
    logical(defBool), intent(in)                        :: lin
    integer(shortInt), intent(inout)                       :: ani
    logical(defBool), intent(in)                        :: doKinetics
    logical(defBool), intent(in)                        :: loud
    logical(defBool), intent(in)                        :: mapTemp
    class(dictionary), pointer, intent(inout), optional :: dictFS
    integer(shortInt), intent(in), optional             :: volPolicy, missPolicy
    logical(defBool), intent(in), optional              :: set2D
    integer(shortInt)                                   :: i
    real(defReal), dimension(6)                         :: bb
    class(field), pointer                               :: field
    character(100), parameter :: Here = 'init (arraysRR_class.f90)'

    call self % XSData % init(db, doKinetics, ani, loud, mapTemp)
    self % nG          = self % XSdata % getNG()
    self % nT          = self % XSdata % getNT()
    self % geom        => geom
    self % nCells      = self % geom % numberOfCells() 

    self % feedback = mapTemp
    if (self % feedback) then
      print *, "Making field pointer"
      field => gr_fieldPtr(gr_fieldIdx(nameTemperature)) ! Need to generalise
      self % tempField => pieceConstantField_CptrCast(field)
    end if
    
    self % lengthPerIt = lengthPerIt
    self % rho = real(rho, defFlt)

    if (present(volPolicy)) then
      self % volPolicy = volPolicy
    else
      self % volPolicy = simAverage
    end if
    if (present(missPolicy)) then
      self % missPolicy = missPolicy
    else
      self % missPolicy = srcPolicy
    end if

    if (present(set2D)) then
      self % set2D = set2D
    else
      self % set2D = .false.
    end if

    ! Assume bounding box of the geometry is filled (and a box)
    ! Can this be relaxed in future?
    bb = self % geom % bounds()
    self % totalVolume = (bb(4) - bb(1)) * (bb(5) - bb(2)) * (bb(6) - bb(3))

    ani = 0
    ! Set simulation type
    if (.not. lin .and. ani == 0) then
      self % simulationType = flatIso
    elseif (lin .and. ani == 0) then
      self % simulationType = linearIso
    elseif (.not. lin .and. ani > 0) then
      self % simulationType = flatAni
    else
      self % simulationType = linearAni
    end if
    if (ani >= 0) self % ani = ani

    ! Allocate and initialise arrays
    allocate(self % scalarFlux(self % nG * self % nCells))
    allocate(self % prevFlux(self % nG * self % nCells))
    allocate(self % fluxScores(2, self % nG * self % nCells))
    allocate(self % power(self % nCells))
    allocate(self % powerScores(2, self % nCells)) 
    allocate(self % source(self % nG * self % nCells))
    allocate(self % volumeTracks(self % nCells))
    allocate(self % lengthSquared(self % nCells))
    allocate(self % allVolumeTracks(self % nCells))
    allocate(self % volume(self % nCells))
    allocate(self % cellHit(self % nCells))
    allocate(self % cellTotalHit(self % nCells))
    allocate(self % cellFound(self % nCells))
    allocate(self % cellPos(3, self % nCells))
    
    self % scalarFlux    = 0.0_defFlt
    self % prevFlux      = 1.0_defFlt
    self % fluxScores    = ZERO
    self % power         = ZERO
    self % powerScores   = ZERO
    self % source        = 0.0_defFlt
    self % volumeTracks  = ZERO
    self % lengthSquared = ZERO
    self % volume        = ZERO
    self % cellHit       = 0
    self % cellTotalHit  = 0
    self % cellFound     = .false.
    self % cellPos       = -INFINITY

    ! Initialise the fixed source if present
    if (present(dictFS)) then
      allocate(self % fixedSource(self % nG * self % nCells))
      call self % initialiseFixedSource(dictFS)
    end if

    ! TODO: allocate linear and anisotropic components, if present
    if (lin) then
      allocate(self % scalarX(self % nCells * self % nG))
      allocate(self % scalarY(self % nCells * self % nG))
      allocate(self % scalarZ(self % nCells * self % nG))
      allocate(self % prevX(self % nCells * self % nG))
      allocate(self % prevY(self % nCells * self % nG))
      allocate(self % prevZ(self % nCells * self % nG))
      allocate(self % sourceX(self % nCells * self % nG))
      allocate(self % sourceY(self % nCells * self % nG))
      allocate(self % sourceZ(self % nCells * self % nG))
      allocate(self % momMat(self % nCells * matSize))
      allocate(self % momTracks(self % nCells * matSize))
      allocate(self % centroid(self % nCells * nDim))
      allocate(self % centroidTracks(self % nCells * nDim))
      allocate(self % xScores(2, self % nG * self % nCells))
      allocate(self % yScores(2, self % nG * self % nCells))
      allocate(self % zScores(2, self % nG * self % nCells))
      
      self % scalarX        = ZERO
      self % scalarY        = ZERO
      self % scalarZ        = ZERO
      self % prevX          = ZERO
      self % prevY          = ZERO
      self % prevZ          = ZERO
      self % sourceX        = 0.0_defFlt
      self % sourceY        = 0.0_defFlt
      self % sourceZ        = 0.0_defFlt
      self % momMat         = ZERO
      self % momTracks      = ZERO
      self % centroid       = ZERO
      self % centroidTracks = ZERO
      self % xScores        = ZERO
      self % yScores        = ZERO
      self % zScores        = ZERO
    end if

    if (ani > 0) then

    end if
    
    ! Initialise OMP locks
    allocate(self % locks(self % nCells))
    do i = 1, self % nCells
#ifdef _OPENMP
      call OMP_init_lock(self % locks(i))
#endif
    end do

  end subroutine init

  !!
  !! Overwrite the active length
  !!
  subroutine setActiveLength(self, lengthPerIt)
    class(arraysRR), intent(inout) :: self
    real(defReal), intent(in)      :: lengthPerIt

    self % lengthPerIt = lengthPerIt

  end subroutine setActiveLength

  !!
  !! Initialise the adjoint source and update nuclear data.
  !! For now, assumes the adjoint is for global variance reduction.
  !!
  subroutine initAdjoint(self)
    class(arraysRR), intent(inout)         :: self
    integer(shortInt)                      :: cIdx
    logical(defBool)                       :: doLinear
    integer(shortInt), save                :: i, g
    real(defFlt), dimension(matSize), save :: invM
    real(defFlt), save                     :: xMom, yMom, zMom
    !$omp threadprivate(i, g, invM, xMom, yMom, zMom)

    doLinear = .false.
    invM = 0.0_defFlt

    call self % xsData % setAdjointXS()

    if (.not. allocated(self % fixedSource)) then
      allocate(self % fixedSource(size(self % scalarFlux)))
    end if
    self % fixedSource = 0.0_defFlt

    if ((self % simulationType == linearIso) .or. &
            (self % simulationType == linearAni)) then
      doLinear = .true.
      if (.not. allocated(self % fixedX)) then
        allocate(self % fixedX(size(self % scalarFlux)))
        allocate(self % fixedY(size(self % scalarFlux)))
        allocate(self % fixedZ(size(self % scalarFlux)))
      end if
      self % fixedX = 0.0_defFlt
      self % fixedY = 0.0_defFlt
      self % fixedZ = 0.0_defFlt
    end if

    ! Create fixed source from the flux scores
    ! Presently assumes the response of interest is global flux
    !$omp parallel do
    do cIdx = 1, self % nCells

      if (.not. self % wasFound(cIdx)) cycle
      if (doLinear) invM = self % invertMatrix(cIdx)

      do g = 1, self % nG

        i = (cIdx - 1) * self % nG + g

        ! Check for inordinately small flux values.
        ! Note, these can have arbitrarily low magnitude.
        ! Maybe should be something more robust.
        if (self % fluxScores(1, i) == ZERO) cycle
        self % fixedSource(i) = real(ONE / self % fluxScores(1, i), defFlt)

        ! Linear source treatment relies on performing a Taylor expansion
        ! of q' = 1/phi = 1/(phi_0 + <gradPhi , (r - r0)>)
        ! = 1/phi_0 - 1/phi^2_0 * <gradPhi , (r - r0)>
        if (doLinear .and. (self % fixedSource(i) > 0)) then
          self % fixedX(i) = invM(xx) * xMom + invM(xy) * yMom + invM(xz) * zMom
          self % fixedY(i) = invM(xy) * xMom + invM(yy) * yMom + invM(yz) * zMom
          self % fixedZ(i) = invM(xz) * xMom + invM(yz) * yMom + invM(zz) * zMom

          self % fixedX(i) = -self % fixedX(i) * self % fixedSource(i) ** 2
          self % fixedY(i) = -self % fixedY(i) * self % fixedSource(i) ** 2
          self % fixedZ(i) = -self % fixedZ(i) * self % fixedSource(i) ** 2

        end if

      end do

    end do
    !$omp end parallel do

    ! Reinitialise arrays to be used during transport
    self % scalarFlux      = ZERO
    self % prevFlux        = ZERO
    self % fluxScores      = ZERO
    self % source          = 0.0_defFlt

    ! Ideally we would have a way of reusing the volume estimators
    ! No compact ideas at the moment, so these will simply be reinitialised
    self % volumeTracks    = ZERO
    self % allVolumeTracks = ZERO
    self % volume          = ZERO

    if (doLinear) then
      self % scalarX        = ZERO
      self % scalarY        = ZERO
      self % scalarZ        = ZERO
      self % prevX          = ZERO
      self % prevY          = ZERO
      self % prevZ          = ZERO
      self % sourceX        = 0.0_defFlt
      self % sourceY        = 0.0_defFlt
      self % sourceZ        = 0.0_defFlt
      self % momMat         = ZERO
      self % momTracks      = ZERO
      self % centroid       = ZERO
      self % centroidTracks = ZERO
      self % xScores        = ZERO
      self % yScores        = ZERO
      self % zScores        = ZERO
    end if

  end subroutine initAdjoint
  
  !!
  !! Initialises fixed sources to be used in the simulation.
  !! Takes a dictionary containing names of materials in the geometry and
  !! source strengths in each energy group and places these in the appropriate
  !! elements of the fixed source vector.
  !!
  !! Also sets source material identities for future use with uncollided calculations. 
  !!
  subroutine initialiseFixedSource(self, dict)
    class(arraysRR), intent(inout)               :: self
    class(dictionary), intent(inout)             :: dict
    character(nameLen),dimension(:), allocatable :: names
    real(defReal), dimension(:), allocatable     :: sourceStrength
    integer(shortInt)                            :: i, nSource, cIdx
    integer(shortInt), save                      :: g, matIdx, idx, id
    logical(defBool)                             :: found
    character(nameLen)                           :: sourceName
    character(nameLen), save                     :: localName
    character(100), parameter :: Here = 'initialiseFixedSource (arraysRR_class.f90)'
    !$omp threadprivate(matIdx, localName, idx, g, id)

    call dict % keys(names)

    nSource = size(names)

    ! Use for uncollided flux sampling
    allocate(self % sourceIdx(nSource))

    ! Cycle through entries of the dictionary
    do i = 1, nSource

      sourceName = names(i)
      call dict % get(sourceStrength, sourceName)

      ! Ensure correct number of energy groups
      if (size(sourceStrength) /= self % nG) call fatalError(Here,'Source '//sourceName//&
              ' has '//numToChar(size(sourceStrength))//' groups rather than '//numToChar(self % nG))

      ! Make sure that the source corresponds to a material present in the geometry
      found = .false.
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells

        id        = cIdx
        matIdx    = self % geom % geom % graph % getMatFromUID(id)
        localName = self % XSData % getName(matIdx)

        if (localName == sourceName) then

          if (.not. found) then
            !$omp critical
            self % sourceIdx(i) = matIdx         
            !$omp end critical
          end if

          found = .true.
          do g = 1, self % nG
            idx = (cIdx - 1) * self % nG + g
            self % fixedSource(idx) = real(sourceStrength(g),defFlt)
          end do

        end if

      end do
      !$omp end parallel do

      if (.not. found) call fatalError(Here,'The source '//trim(sourceName)//' does not correspond to '//&
              'any material found in the geometry.')

    end do

  end subroutine initialiseFixedSource


  !!
  !! Return a pointer to the flux vector for a given cell
  !!
  subroutine getFluxPointer(self, cIdx, fluxVec)
    class(arraysRR), intent(in), target              :: self
    integer(shortInt), intent(in)                    :: cIdx
    real(defFlt), dimension(:), pointer, intent(out) :: fluxVec
    integer(shortInt)                                :: baseIdx1, baseIdx2

    baseIdx1 = self % nG * (cIdx - 1) + 1
    baseIdx2 = self % nG * cIdx
    fluxVec => self % scalarFlux(baseIdx1:baseIdx2)

  end subroutine getFluxPointer

  !!
  !! Return a pointer to the flux spatial moment vectors for a given cell
  !!
  subroutine getFluxXYZPointers(self, cIdx, xVec, yVec, zVec)
    class(arraysRR), intent(in), target               :: self
    integer(shortInt), intent(in)                     :: cIdx
    real(defFlt), dimension(:), pointer, intent(out) :: xVec
    real(defFlt), dimension(:), pointer, intent(out) :: yVec
    real(defFlt), dimension(:), pointer, intent(out) :: zVec
    integer(shortInt)                                 :: baseIdx1, baseIdx2

    baseIdx1 = self % nG * (cIdx - 1) + 1
    baseIdx2 = self % nG * cIdx
    xVec => self % scalarX(baseIdx1:baseIdx2)
    yVec => self % scalarY(baseIdx1:baseIdx2)
    zVec => self % scalarZ(baseIdx1:baseIdx2)

  end subroutine getFluxXYZPointers
  
  !!
  !! Return a pointer to the nuclear data object
  !!
  function getDataPointer(self) result(dataPtr)
    class(arraysRR), intent(in), target :: self
    class(dataRR), pointer              :: dataPtr

    dataPtr => self % XSData

  end function getDataPointer
  
  !!
  !! Return a pointer to the geometry object
  !!
  function getGeomPointer(self) result(geomPtr)
    class(arraysRR), intent(in), target :: self
    class(geometryStd), pointer         :: geomPtr

    geomPtr => self % geom

  end function getGeomPointer

  !!
  !! Return a pointer to the geometry object
  !!
  function getTempPointer(self) result(tempPtr)
    class(arraysRR), intent(in), target :: self
    class(pieceConstantField), pointer         :: tempPtr

    tempPtr => self % tempField

  end function getTempPointer
  
  
  !!
  !! Return a pointer to the source vector for a given cell
  !!
  subroutine getSourcePointer(self, cIdx, sourceVec)
    class(arraysRR), intent(in), target              :: self
    integer(shortInt), intent(in)                    :: cIdx
    real(defFlt), dimension(:), pointer, intent(out) :: sourceVec
    integer(shortInt)                                :: baseIdx1, baseIdx2

    baseIdx1 = self % nG * (cIdx - 1) + 1
    baseIdx2 = self % nG * cIdx
    sourceVec => self % source(baseIdx1:baseIdx2)

  end subroutine getSourcePointer

  !!
  !! Return pointers to the source gradient vectors for a given cell
  !!
  subroutine getSourceXYZPointers(self, cIdx, xVec, yVec, zVec)
    class(arraysRR), intent(in), target              :: self
    integer(shortInt), intent(in)                    :: cIdx
    real(defFlt), dimension(:), pointer, intent(out) :: xVec
    real(defFlt), dimension(:), pointer, intent(out) :: yVec
    real(defFlt), dimension(:), pointer, intent(out) :: zVec
    integer(shortInt)                                :: baseIdx1, baseIdx2

    baseIdx1 = self % nG * (cIdx - 1) + 1
    baseIdx2 = self % nG * cIdx
    xVec => self % sourceX(baseIdx1:baseIdx2)
    yVec => self % sourceY(baseIdx1:baseIdx2)
    zVec => self % sourceZ(baseIdx1:baseIdx2)

  end subroutine getSourceXYZPointers
  
  !!
  !! Return source value given cell and group
  !!
  function getSource(self, cIdx, g) result(src)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defFlt)                  :: src

    src = self % source(self % nG * (cIdx - 1) + g)

  end function getSource

  !!
  !! Return fixed source value given cell and group
  !!
  elemental function getFixedSource(self, cIdx, g) result(src)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defFlt)                  :: src

    if (allocated(self % fixedSource)) then
      src = self % fixedSource(self % nG * (cIdx - 1) + g)
    else
      src = 0.0_defFlt
    end if

  end function getFixedSource
  
  !!
  !! Return previous flux value given cell and group
  !!
  function getPrevFlux(self, cIdx, g) result(flux)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defFlt)                  :: flux

    flux = self % prevFlux(self % nG * (cIdx - 1) + g)

  end function getPrevFlux
  
  !!
  !! Return final flux value given cell and group
  !!
  function getFluxScore(self, cIdx, g) result(flux)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defReal)                 :: flux

    flux = self % fluxScores(1, self % nG * (cIdx - 1) + g)

  end function getFluxScore
  
  !!
  !! Return final flux standard deviation given cell and group
  !! Will return square of flux scores if called before finaliseFluxScores
  !!
  function getFluxSD(self, cIdx, g) result(flux)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defReal)                 :: flux

    flux = self % fluxScores(2, self % nG * (cIdx - 1) + g)

  end function getFluxSD

  !!
  !! Return final flux moment values given cell and group
  !!
  pure function getFluxMoments(self, cIdx, g) result(flux)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defReal), dimension(3)   :: flux
    integer(shortInt)             :: idx

    idx = self % nG * (cIdx - 1) + g
    flux(1) = self % xScores(1, idx)
    flux(2) = self % yScores(1, idx)
    flux(3) = self % zScores(1, idx)

  end function getFluxMoments

  !!
  !! Return final flux moment standard deviations given cell and group
  !! Will return square of moment scores if called before finaliseFluxScores
  !!
  pure function getFluxMomentSDs(self, cIdx, g) result(fluxSD)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defReal), dimension(3)   :: fluxSD
    integer(shortInt)             :: idx

    idx = self % nG * (cIdx - 1) + g
    fluxSD(1) = self % xScores(2, idx)
    fluxSD(2) = self % yScores(2, idx)
    fluxSD(3) = self % zScores(2, idx)

  end function getFluxMomentSDs

  !!
  !! Return final power value given cell
  !!
  function getPowerScore(self, cIdx) result(power)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    real(defReal)                 :: power

    power = self % powerScores(1, cIdx)

    !print *, power
  end function getPowerScore

  !!
  !! Return final power standard deviation given cell
  !! Will return square of power scores if called before finaliseFluxScores
  !!
  function getPowerSD(self, cIdx) result(power)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    real(defReal)                 :: power

    power = self % powerScores(2, cIdx)

  end function getPowerSD

  !!
  !! Return cell position given cell ID
  !!
  !function getCellPos(self, cIdx, g) result(x)

  function getCellPos(self, cIdx) result(x)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    !integer(shortInt), intent(in) :: g
    real(defReal), dimension(3)   :: x

    x = self % cellPos(1:3, cIdx)

  end function getCellPos

  !!
  !! Return cell centroid given cell ID
  !!
  pure function getCentroid(self, cIdx) result(cent)
    class(arraysRR), intent(in)    :: self
    integer(shortInt), intent(in)  :: cIdx
    real(defReal), dimension(nDim) :: cent
    integer(shortInt)              :: idx0, idx1

    idx0 = nDim * (cIdx - 1) + 1
    idx1 = nDim * cIdx
    cent = self % centroid(idx0:idx1)

  end function getCentroid

  !!
  !! Return moment matrix given cell ID
  !!
  pure function getMomentMatrix(self, cIdx) result(mat)
    class(arraysRR), intent(in)       :: self
    integer(shortInt), intent(in)     :: cIdx
    real(defReal), dimension(matSize) :: mat
    integer(shortInt)              :: idx0, idx1

    idx0 = matSize * (cIdx - 1) + 1
    idx1 = matSize * cIdx
    mat = self % momMat(idx0:idx1)

  end function getMomentMatrix
  
  !!
  !! Return the simulation type
  !!
  function getSimulationType(self) result(simType)
    class(arraysRR), intent(in) :: self
    integer(shortInt)           :: simType

    simType = self % simulationType

  end function getSimulationType

  !!
  !! Increment the local volume estimate in cell cIdx.
  !! Assumes this is being called inside a lock for thread privacy.
  !!
  subroutine incrementVolume(self, cIdx, length)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: cIdx
    real(defReal), intent(in)      :: length     
    
    self % volumeTracks(cIdx) = self % volumeTracks(cIdx) + length
  
  end subroutine incrementVolume
  
  !!
  !! Increment the sum of length squared in cell cIdx.
  !! Assumes this is being called inside a lock for thread privacy.
  !!
  subroutine incrementLengthSquared(self, cIdx, length)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: cIdx
    real(defReal), intent(in)      :: length     
    
    self % lengthSquared(cIdx) = self % lengthSquared(cIdx) + length * length
  
  end subroutine incrementLengthSquared

  !!
  !! Increment the local centroid estimate in cell cIdx.
  !! rL is the tracklength-weighted centroid
  !! Assumes this is being called inside a lock for thread privacy.
  !!
  subroutine incrementCentroid(self, cIdx, rL)
    class(arraysRR), intent(inout)             :: self
    integer(shortInt), intent(in)              :: cIdx
    real(defReal), dimension(nDim), intent(in) :: rL
    integer(shortInt)                          :: idx0, idx1

    idx0 = nDim * (cIdx - 1) + 1
    idx1 = nDim * cIdx
    self % centroidTracks(idx0:idx1) = self % centroidTracks(idx0:idx1) + rL

  end subroutine incrementCentroid
  
  !!
  !! Increment the local moment matrix estimate in cell cIdx.
  !! mat is the tracklength-weighted matrix
  !! Assumes this is being called inside a lock for thread privacy.
  !!
  subroutine incrementMoments(self, cIdx, mat)
    class(arraysRR), intent(inout)                :: self
    integer(shortInt), intent(in)                 :: cIdx
    real(defReal), dimension(matSize), intent(in) :: mat
    integer(shortInt)                             :: idx0, idx1

    idx0 = matSize * (cIdx - 1) + 1
    idx1 = matSize * cIdx
    self % momTracks(idx0:idx1) = self % momTracks(idx0:idx1) + mat

  end subroutine incrementMoments
  
  !!
  !! Check if a cell has been hit
  !!
  elemental function hasHit(self, cIdx) result (hit)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt)             :: hit
    
    hit = self % cellHit(cIdx)
  
  end function hasHit

  !!
  !! Check if a cell has been hit
  !!
  elemental function wasHit(self, cIdx) result (hit)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    logical(defBool)              :: hit

    hit = (self % cellHit(cIdx) == 1)

  end function wasHit

  !!
  !! Has a cell ever been found?
  !!
  elemental function wasFound(self, cIdx) result(found)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    logical(defBool)              :: found

    found = (self % cellTotalHit(cIdx) > 0)

  end function wasFound
  
  !!
  !! Hit a cell 
  !!
  subroutine hitCell(self, cIdx)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: cIdx
    
    self % cellHit(cIdx) = 1
    self % cellTotalHit(cIdx) = self % cellTotalHit(cIdx) + 1
  
  end subroutine hitCell

  !!
  !! Return the cell hit rate for the given iteration
  !!
  function getCellHitRate(self, it) result(hitRate)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: it
    integer(shortInt)            :: totalHit, realCells
    real(defReal)                :: hitRate

    !totalHit = sum(self % cellHit)
    !hitRate  = real(totalHit, defReal) / self % nCells
    ! Reset averages after iteration 20
    if (it == 21) then
      self % averageHit = ZERO
      self % iterations = 0
    end if

    totalHit = sum(self % cellHit)
    if (it > 20) then
      realCells = count(self % cellTotalHit > 0)
    else
      realCells = self % nCells
    end if
    hitRate = real(totalHit,defReal) / realCells

    self % averageHit = self % averageHit + hitRate
    self % iterations = self % iterations + 1
    

  end function getCellHitRate

  !!
  !! Return the simulation average cell hit rate
  !!
  function getAverageHitRate(self) result(hitRate)
    class(arraysRR), intent(in)   :: self
    real(defReal)                 :: hitRate

    hitRate = self % averageHit / self % iterations

  end function getAverageHitRate

  !! 
  !! Wipe cell hits
  !!
  subroutine wipeCellHits(self)
    class(arraysRR), intent(inout) :: self

    self % cellHit = 0

  end subroutine wipeCellHits
  
  !!
  !! Has a cell ever been found?
  !!
  function found(self, cIdx) result(wasFound)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    logical(defBool)              :: wasFound
    
    wasFound = self % cellFound(cIdx)
  
  end function found

  !!
  !! Note that a new cell has been found
  !!
  subroutine newFound(self, cIdx, r)
    class(arraysRR), intent(inout)          :: self
    integer(shortInt), intent(in)           :: cIdx
    real(defReal), dimension(3), intent(in) :: r

    !$omp critical 
    self % cellFound(cIdx) = .true.
    self % cellPos(:,cIdx) = r
    !$omp end critical

  end subroutine newFound

  !!
  !! Return number of energy groups used
  !!
  function getNG(self) result(nG)
    class(arraysRR), intent(in) :: self
    integer(shortInt)           :: nG

    nG = self % nG

  end function getNG

  !!
  !! Check if a cell has an inhomogeneous source
  !!
  elemental function hasFixedSource(self, cIdx) result (hasSrc)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    logical(defBool)              :: hasSrc
    integer(shortInt)             :: idx1, idx2

    if (allocated(self % fixedSource)) then
      idx1 = self % nG * (cIdx - 1) + 1
      idx2 = self % nG * cIdx
      ! Take an absolute value in case of (possibly desirable?) negative sources
      hasSrc = any(abs(self % fixedSource(idx1:idx2)) > 0.0_defFlt)
    else
      hasSrc = .false.
    end if

  end function hasFixedSource

  !!
  !! Return volume given cell ID
  !!
  elemental function getVolume(self, cIdx) result(vol)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    real(defReal)                 :: vol

    vol = self % volume(cIdx)

  end function getVolume  

  !!
  !! Set the OMP lock in a given cell
  !!
  subroutine setLock(self, cIdx) 
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: cIdx

#ifdef _OPENMP
    call OMP_set_lock(self % locks(cIdx))
#endif

  end subroutine setLock
  
  !!
  !! Unset the OMP lock in a given cell
  !!
  subroutine unsetLock(self, cIdx) 
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: cIdx

#ifdef _OPENMP
    call OMP_unset_lock(self % locks(cIdx))
#endif

  end subroutine unsetLock

  !!
  !! Calls appropriate normalise flux and volume subroutines
  !!
  subroutine normaliseFluxAndVolume(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it
    character(100), parameter :: Here = 'normaliseFluxAndVolume (arraysRR_class.f90)'

    select case(self % simulationType)
      case(flatIso)
        call self % normaliseFluxAndVolumeFlatIso(it)
      case(linearIso)
        call self % normaliseFluxAndVolumeLinIso(it)
      case default
        call fatalError(Here,'Unsupported simulation type requested')
    end select

  end subroutine normaliseFluxAndVolume

  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolumeFlatIso(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it
    real(defReal)                             :: norm, normIt
    real(defReal), save                       :: vol, volAve, volNaive, D
    real(defFlt), save                        :: sigGG, tot, temp
    real(defFlt), dimension(self % nG)        :: total, scatterVec
    type(particle), save                      :: s
    integer(shortInt), save                   :: g, matIdx, idx, baseIdx
    integer(shortInt)                         :: cIdx
    logical(defBool), save                    :: hit, isSrc, smallCell
    real(defReal)                             :: pos(3), dir(3)
    character(100), parameter :: Here = 'normaliseFluxAndVolumeFlatIso (arraysRR_class.f90)'
    !$omp threadprivate(vol, idx, g, matIdx, baseIdx, sigGG, D, hit, isSrc, volAve, volNaive, tot, smallCell, temp, s)

    norm = ONE / self % lengthPerIt
    normIt = ONE / (self % lengthPerIt * it)
    
    !$omp parallel do private(total, scatterVec, pos, dir) 
    cellLoop: do cIdx = 1, self % nCells
      matIdx = self % geom % geom % graph % getMatFromUID(cIdx) 

      if (matIdx > self % XSData % getNMat()) then
         baseIdx = self % nG * (cIdx - 1)
         do g = 1, self % nG
           idx = baseIdx + g
           self % scalarFlux(idx) = 0.0_defFlt
         end do
         cycle
      end if
     
      if (matIdx < 1 .or. matIdx > self % XSData % getNMat()) then
        print *, "BAD MATIDX"
        print *, "cIdx      =", cIdx
        print *, "matIdx    =", matIdx
        print *, "nMat      =", self % XSData % getNMat()
        print *, "cellPos   =", self % cellPos(:,cIdx)
        print *, "cellHit   =", self % cellTotalHit(cIdx)
        !call fatalError(Here,'Invalid material index before getTotalXS')
      end if

      hit = self % wasHit(cIdx)
      isSrc = self % hasFixedSource(cIdx)

      ! Is the cell hit frequently?
      !smallCell = (real(self % cellTotalHit(cIdx) / it, defReal) < 1.5)
      smallCell = (real(self % cellTotalHit(cIdx), defReal) / real(it, defReal) < 1.5_defReal)
      
      !! Compute various volume types
      ! Actual integral volume
      self % allVolumeTracks(cIdx) = self % allVolumeTracks(cIdx) + self % volumeTracks(cIdx)
      self % volume(cIdx) = self % allVolumeTracks(cIdx) * normIt
      volAve = self % volume(cIdx)
      ! Cycle-wise volume
      volNaive = self % volumeTracks(cIdx) * norm

      ! Decide volume to use
      select case(self % volPolicy)
        case(simAverage)
          vol = volAve
        case(naive)
          vol = volNaive
        case(hybrid)
          if (isSrc) then
            vol = volNaive
          else
            vol = volAve
          end if
        case default
          call fatalError(Here,'Unsupported volume handling requested')
      end select

      if (smallCell) vol = volNaive
      
      ! Reset cycle-wise estimator
      self % volumeTracks(cIdx) = ZERO
            
      !call self % XSData % getTotalPointer(matIdx, total)
      ! Fudge particle state to search tally map
      pos = self % cellPos(:,cIdx)
      dir = [0.5,0.5,0.5]

      call s % coords % init(pos, dir)
      s % coords % matIdx = matIdx

      !print *, "FEEDBACK"
      !print *, self % feedback

      if (self % feedback) then

        print *, self % feedback
        print *, "Feedback present"
        temp = self % tempField % at(s % coords)

        total = self % XSData % getTotalXS(matIdx, temp)

      else

        total = self % XSData % getTotalXS(matIdx)

      end if
            
      groupLoop: do g = 1, self % nG

        idx = self % nG * (cIdx - 1) + g
        tot = total(g)

        ! Route for non-void materials
        if (matIdx <= self % XSData % getNMat() .and. tot > 0) then
          if (hit) then
     
            ! Can hit a cell but with a tiny volume, such that 
            ! things break a bit - would rather remove this arbitrary
            ! check in future
            if (vol < volume_tolerance) then
              self % scalarFlux(idx) = ZERO
              cycle groupLoop
            end if

            self % scalarFlux(idx) = self % scalarFlux(idx) * norm 
            self % scalarFlux(idx) = self % scalarFlux(idx) / (vol * tot)

            ! Presumes non-zero total XS
            !sigGG = self % XSData % getScatterXS(matIdx, g, g)
            if (self % feedback) then
              scatterVec = self % XSData % getScatterXS(matIdx, g, temp)
            else
              scatterVec = self % XSData % getScatterXS(matIdx, g)
            end if

            sigGG = scatterVec(g)

            if ((sigGG < 0) .and. (total(g) > 0)) then
              D = -real(self % rho * sigGG / tot, defReal)
            else
              D = ZERO
            end if

            self % scalarFlux(idx) =  (self % scalarFlux(idx) + self % source(idx) / tot &
                  + D * self % prevFlux(idx) ) / (1 + D)

          else
            
            ! Decide flux treatment to use on missing a cell
            select case(self % missPolicy)
              case(srcPolicy)
                self % scalarFlux(idx) = self % source(idx) / tot
              case(prevPolicy)
                self % scalarFlux(idx) = self % prevFlux(idx)
              case(hybrid)
                if (isSrc) then
                  self % scalarFlux(idx) = self % prevFlux(idx)
                else
                  self % scalarFlux(idx) = self % source(idx) / tot
                end if
              case default
                call fatalError(Here,'Unsupported miss handling requested')
            end select

          end if
        
        ! Alternatively, handle unidentified/void regions
        else
          
          if (vol < volume_tolerance) then
            self % scalarFlux(idx) = ZERO
            cycle groupLoop
          end if

          if (hit) then
            self % scalarFlux(idx) = self % scalarFlux(idx) * norm / vol
          else
            self % scalarFlux(idx) = self % prevFlux(idx)
          end if
        end if

      end do groupLoop

    end do cellLoop
    !$omp end parallel do

  end subroutine normaliseFluxAndVolumeFlatIso

  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source for linear isotropic sources
  !!
  subroutine normaliseFluxAndVolumeLinIso(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it
    real(defReal)                             :: norm, normVol
    real(defReal), save                       :: vol, invVol, volNaive, volAve, norm_V
    real(defFlt), save                        :: D, sigGG, temp
    type(particle), save                      :: s
    real(defFlt), allocatable, dimension(:), save  :: total
    real(defFlt), allocatable, dimension(:), save        :: scatterVec
    integer(shortInt)                         :: cIdx
    integer(shortInt), save                   :: g, matIdx, idx, dIdx, mIdx
    logical(defBool), save                    :: hit, isSrc, smallCell
    real(defReal)                             :: pos(3), dir(3)
    character(100), parameter :: Here = 'normaliseFluxAndVolumeLinearIso (arraysRR_class.f90)'
    !$omp threadprivate(s, total, temp, vol, idx, mIdx, dIdx, g, matIdx, invVol, norm_V, D, sigGG, hit, isSrc, &
    !$omp volNaive, volAve, smallCell, scatterVec)

    norm = ONE / self % lengthPerIt
    normVol = ONE / (self % lengthPerIt * it)
    
    !$omp parallel
    if (.not. allocated(total)) allocate(total(self % nG))
    if (.not. allocated(scatterVec)) allocate(scatterVec(self % nG))
    !$omp end parallel

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx)
      dIdx = (cIdx - 1) * nDim
      mIdx = (cIdx - 1) * matSize

      ! Fudge particle state to search tally map
      pos = self % cellPos(:,cIdx)
      dir = [0.5,0.5,0.5]

      call s % coords % init(pos, dir)
      s % coords % matIdx = matIdx

      if (self % feedback) then
        temp = self % tempField % at(s % coords)

        total = self % XSData % getTotalXS(matIdx, temp)

      else

        total = self % XSData % getTotalXS(matIdx)

      end if

      hit = self % wasHit(cIdx)
      isSrc = self % hasFixedSource(cIdx)
      ! Is the cell hit frequently?
      !smallCell = (real(self % cellTotalHit(cIdx) / it, defReal) < 1.5)
      smallCell = (real(self % cellTotalHit(cIdx), defReal) / real(it, defReal) < 1.5_defReal)

      ! Compute various volume types
      ! Actual integral volume
      self % allVolumeTracks(cIdx) = self % allVolumeTracks(cIdx) + &
              self % volumeTracks(cIdx)
      self % volume(cIdx) = self % allVolumeTracks(cIdx) * normVol
      volAve = self % volume(cIdx)
      ! Iteration-wise volume
      volNaive = self % volumeTracks(cIdx) * norm
      
      ! Decide volume to use
      select case(self % volPolicy)
        case(naive)
          vol = volNaive
          norm_V = ONE / self % volumeTracks(cIdx)
        case(simAverage)
          vol = volAve
          norm_V = it / self % allVolumeTracks(cIdx)
        case(hybrid)
          if (isSrc) then
            vol = volNaive
            norm_V = ONE / self % volumeTracks(cIdx)
          else
            vol = volAve
            norm_V = it / self % allVolumeTracks(cIdx)
          end if
        case default
          call fatalError(Here,'Unsupported volume handling requested')
      end select

      if (smallCell) then
        vol = volNaive
        norm_V = ONE / self % volumeTracks(cIdx)
      end if
      
      ! Reset cycle-wise estimator
      self % volumeTracks(cIdx) = ZERO
       
      ! Update geometric information provided volume has been visited
      if (self % allVolumeTracks(cIdx) > ZERO) then
      
        invVol = ONE / self % allVolumeTracks(cIdx)
        
        ! Update centroids
        self % centroid(dIdx + x) =  self % centroidTracks(dIdx + x) * invVol
        self % centroid(dIdx + y) =  self % centroidTracks(dIdx + y) * invVol
        self % centroid(dIdx + z) =  self % centroidTracks(dIdx + z) * invVol
      
        ! Update spatial moments
        self % momMat(mIdx + xx) = self % momTracks(mIdx + xx) * invVol
        self % momMat(mIdx + xy) = self % momTracks(mIdx + xy) * invVol
        self % momMat(mIdx + xz) = self % momTracks(mIdx + xz) * invVol
        self % momMat(mIdx + yy) = self % momTracks(mIdx + yy) * invVol
        self % momMat(mIdx + yz) = self % momTracks(mIdx + yz) * invVol
        self % momMat(mIdx + zz) = self % momTracks(mIdx + zz) * invVol

      end if

      !call self % XSData % getTotalPointer(matIdx, total)

      groupLoop: do g = 1, self % nG

        idx = self % nG * (cIdx - 1) + g
        
        if (matIdx <= self % XSData % getNMat() .and. total(g) > 0) then
          if (hit) then

            ! Can hit a cell but with a tiny volume, such that 
            ! things break a bit - would rather remove this arbitrary
            ! check in future
            if (vol < volume_tolerance) then
              self % scalarFlux(idx) = ZERO
              self % scalarX(idx) = ZERO
              self % scalarY(idx) = ZERO
              self % scalarZ(idx) = ZERO
              cycle groupLoop
            end if
      
            self % scalarFlux(idx) = self % scalarFlux(idx) * norm_V / total(g)
            self % scalarX(idx) = self % scalarX(idx) * norm_V / total(g)
            self % scalarY(idx) = self % scalarY(idx) * norm_V / total(g)
            self % scalarZ(idx) = self % scalarZ(idx) * norm_V / total(g)
          
            ! Apply the standard MoC post-sweep treatment and
            ! stabilisation for negative XSs
            ! Presumes non-zero total XS
            !sigGG = self % XSData % getScatterXS(matIdx, g, g)

            if (self % feedback) then
              scatterVec = self % XSData % getScatterXS(matIdx, g, temp)
            else
              scatterVec = self % XSData % getScatterXS(matIdx, g)
            end if

            sigGG = scatterVec(g)

            if (sigGG < 0) then
              D = -self % rho * sigGG / total(g)
            else
              D = 0.0_defFlt
            end if
           
            self % scalarFlux(idx) =  (self % scalarFlux(idx) + self % source(idx)/total(g) &
                  + D * self % prevFlux(idx) ) / (1 + D)

          else
            ! Decide flux treatment to use
            associate(mat => self % momMat((mIdx + 1):(mIdx + matSize)))
            select case(self % missPolicy)
              ! Note: this is policy to use the source, not policy for hitting a fixed source
              case(srcPolicy)
                self % scalarFlux(idx) = self % source(idx) / total(g)
                ! OPENMC SETS MOMENTS TO ZERO
                !self % scalarX(idx) = 0.0_defFlt
                !self % scalarY(idx) = 0.0_defFlt
                !self % scalarZ(idx) = 0.0_defFlt
                ! Need to multiply source gradients by moment matrix
                self % scalarX(idx) = real(mat(xx) *self % sourceX(idx) + &
                        mat(xy) * self % sourceY(idx) + mat(xz) * self % sourceZ(idx),defFlt)/ total(g)
                self % scalarY(idx) = real(mat(xy) *self % sourceX(idx) + &
                        mat(yy) * self % sourceY(idx) + mat(yz) * self % sourceZ(idx),defFlt)/ total(g)
                self % scalarZ(idx) = real(mat(xz) *self % sourceX(idx) + &
                        mat(yz) * self % sourceY(idx) + mat(zz) * self % sourceZ(idx),defFlt)/ total(g)
              case(prevPolicy)
                self % scalarFlux(idx) = self % prevFlux(idx)
                self % scalarX(idx) = self % prevX(idx)
                self % scalarY(idx) = self % prevY(idx)
                self % scalarZ(idx) = self % prevZ(idx)
              case(hybrid)
                if (isSrc) then
                  self % scalarFlux(idx) = self % prevFlux(idx)
                  self % scalarX(idx) = self % prevX(idx)
                  self % scalarY(idx) = self % prevY(idx)
                  self % scalarZ(idx) = self % prevZ(idx)
                else
                  self % scalarFlux(idx) = self % source(idx) / total(g)
                  ! OPENMC SETS MOMENTS TO ZERO
                  !self % scalarX(idx) = 0.0_defFlt
                  !self % scalarY(idx) = 0.0_defFlt
                  !self % scalarZ(idx) = 0.0_defFlt
                  ! Need to multiply source gradients by moment matrix
                  self % scalarX(idx) = real(mat(xx) *self % sourceX(idx) + &
                        mat(xy) * self % sourceY(idx) + mat(xz) * self % sourceZ(idx),defFlt)/ total(g)
                  self % scalarY(idx) = real(mat(xy) *self % sourceX(idx) + &
                        mat(yy) * self % sourceY(idx) + mat(yz) * self % sourceZ(idx),defFlt)/ total(g)
                  self % scalarZ(idx) = real(mat(xz) *self % sourceX(idx) + &
                        mat(yz) * self % sourceY(idx) + mat(zz) * self % sourceZ(idx),defFlt)/ total(g)
                end if
              case default
                call fatalError(Here,'Unsupported miss handling requested')
            end select
            end associate
          end if

        else

          ! Apply void treatment
          if (vol < volume_tolerance) then
            self % scalarFlux(idx) = ZERO
            self % scalarX(idx) = ZERO
            self % scalarY(idx) = ZERO
            self % scalarZ(idx) = ZERO
            cycle groupLoop
          end if

          if (hit) then
            self % scalarFlux(idx) = self % scalarFlux(idx) * norm / vol
            self % scalarX(idx) = ZERO
            self % scalarY(idx) = ZERO
            self % scalarZ(idx) = ZERO
          else
            self % scalarFlux(idx) = self % prevFlux(idx)
            self % scalarX(idx) = ZERO
            self % scalarY(idx) = ZERO
            self % scalarZ(idx) = ZERO
          end if

        end if
        
        ! For stability while still accumulating geometric info
        if (it < 10) then
          self % scalarX(idx) = ZERO
          self % scalarY(idx) = ZERO
          self % scalarZ(idx) = ZERO
        end if

      end do groupLoop

    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolumeLinIso

  
  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source for flat anisotropic sources
  !!
  subroutine normaliseFluxAndVolumeFlatAni(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it

  end subroutine normaliseFluxAndVolumeFlatAni
  
  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source for Linear sources with flat
  !! anisotropic sources
  !!
  subroutine normaliseFluxAndVolumeLIFA(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it

  end subroutine normaliseFluxAndVolumeLIFA

  !!
  !! Calculates thermal power from flux distribution
  !! 
  subroutine calculatePower(self, targetPow)

    class(arraysRR), intent(inout) :: self
    real(defReal), intent(in)      :: targetPow

    real(defReal), save            :: vol
    real(defFlt), save             :: pow

    real(defFlt), dimension(self % nG):: fissXS

    type(particle), save           :: s

    integer(shortInt), save        :: g, matIdx, idx
    real(defFlt), save             :: temp

    integer(shortInt)              :: cIdx

    real(defReal)                  :: pos(3), dir(3)

    character(100), parameter :: Here = 'calculatePower (arraysRR_class.f90)'

    !$omp threadprivate(vol, pow, idx, g, matIdx, s, temp)

    !!$omp parallel do
    !$omp parallel do private(fissXS, pos, dir)
    do cIdx = 1, self % nCells

      matIdx = self % geom % geom % graph % getMatFromUID(cIdx)

      vol = self % volume(cIdx)

      ! Save effort by skipping if volume is too small
      if (vol < volume_tolerance) then
        self % power(cIdx) = 0.0_defFlt
        cycle
      end if

      ! Fudge particle state to search tally map
      pos = self % cellPos(:,cIdx)
      dir = [0.5,0.5,0.5]

      call s % coords % init(pos, dir)
      s % coords % matIdx = matIdx

      if (self % feedback) then

        temp = self % tempField % at(s % coords)

        fissXS = self % XSData % getFissionXS(matIdx, temp)

      else

        fissXS = self % XSData % getFissionXS(matIdx)

      end if

      self % power(cIdx) = 0.0_defFlt

      do g = 1, self % nG

        idx = self % nG * (cIdx - 1) + g

        pow = self % scalarFlux(idx) &
             * energyPerFission &
             * self % volume(cIdx) &
             * joulesPerMeV &
             * fissXS(g)

        self % power(cIdx) = self % power(cIdx) + pow

      end do

      if (self % power(cIdx) < 0.0_defFlt) then
        print *, "Negative power found"
      end if

    end do
    !$omp end parallel do

    self % power = self % power * (targetPow / sum(self % power))

  end subroutine calculatePower  
  
  !!
  !! Update all sources given a prevFlux
  !! This nesting allows using combined OMP + SIMD
  !!
  subroutine updateSource(self, ONE_KEFF, it)
    class(arraysRR), intent(inout) :: self
    real(defReal), intent(in)      :: ONE_KEFF
    integer(shortInt), intent(in)  :: it
    real(defFlt)                   :: ONE_K
    integer(shortInt)              :: cIdx
    character(100), parameter      :: Here = 'updateSource (arraysRR_class.f90)'

    ONE_K = real(ONE_KEFF, defFlt)

    !print *, self % cellPos
    select case(self % simulationType)
      case(flatIso)
        !$omp parallel do 
        do cIdx = 1, self % nCells
          call self % sourceUpdateKernelFlatIso(cIdx, ONE_K, it)
        end do
        !$omp end parallel do
      case(linearIso)
        !$omp parallel do 
        do cIdx = 1, self % nCells
          call self % sourceUpdateKernelLinIso(cIdx, ONE_K, it)
        end do
        !$omp end parallel do
      case default
        call fatalError(Here,'Unsupported simulation type requested')
    end select

  end subroutine updateSource

  !!
  !! Kernel to update sources given a cell index
  !!
  subroutine sourceUpdateKernelFlatIso(self, cIdx, ONE_KEFF, it)
    class(arraysRR), target, intent(inout)   :: self
    integer(shortInt), intent(in)            :: cIdx
    real(defFlt), intent(in)                 :: ONE_KEFF
    integer(shortInt), intent(in)            :: it
    real(defFlt)                             :: scatter, fission, chi, scatterXS, temp, vol
    !real(defFlt), dimension(:), pointer      :: nuFission chi, scatterXS, fluxVec 
    real(defFlt), dimension(:), pointer      :: fluxVec
    real(defFlt), dimension(self % nG)       :: nuFission, scatterVec
    type(particle)                           :: s
    integer(shortInt)                        :: matIdx, g, gIn, baseIdx, idx, sIdx1, sIdx2
    real(defReal)                            :: pos(3), dir(3)
    
    if (it>1) then
      vol = self % volume(cIdx)
      if (vol < volume_tolerance) return
    end if

    ! Identify material
    matIdx = self % geom % geom % graph % getMatFromUID(cIdx) 
    !print *, 'source material has matidx ' // numToChar(matIdx)

    ! Guard against void cells
    if (matIdx > self % XSData % getNMat()) then
      baseIdx = self % nG * (cIdx - 1)
      do g = 1, self % nG
        idx = baseIdx + g
        self % source(idx) = 0.0_defFlt
      end do
      return
    end if

    if ((matIdx < 1) .or. (matIdx >  self % XSData % getNMat())) then
      print *, "matIdx = " // numToChar(matIdx) // " in sourceUpdate"
    end if

    ! Fudge a particle state to search tally map
    pos = self % cellPos(:,cIdx)
    dir = [0.5,0.5,0.5]
    call s % coords % init(pos, dir)
    s % coords % matIdx = matIdx
    if (self % feedback) then
      if (it == 1) then
        temp = self % XSData % getTemperature(matIdx,2)
      else
        temp = self % tempField % at(s % coords)
      end if
    end if

    baseIdx = self % nG * (cIdx - 1)
    fluxVec => self % prevFlux((baseIdx + 1):(baseIdx + self % nG))

    ! Calculate fission source
    fission = 0.0_defFlt
    if (self % feedback) then
      nuFission = self % XSData % getNuFissXS(matIdx, temp)
    else 
      nuFission = self % XSData % getNuFissXS(matIdx)
    end if
    
    !$omp simd reduction(+:fission)
    do gIn = 1, self % nG
      fission = fission + fluxVec(gIn) * nuFission(gIn)
    end do
    fission = fission * ONE_KEFF

    do g = 1, self % nG

      sIdx1 = self % nG * (g - 1) + 1
      sIdx2 = self % nG * g

      ! Calculate scattering source
      scatter = 0.0_defFlt

      if (self % feedback) then
        scatterVec = self % XSData % getScatterXS(matIdx, g, temp)
      else
        scatterVec = self % XSData % getScatterXS(matIdx, g)
      end if

      !$omp simd reduction(+:scatter)
      do gIn = 1, self % nG
        scatter = scatter + fluxVec(gIn) * scatterVec(gIn)
      end do

      ! Output index
      idx = baseIdx + g

      chi = self % XSData % getChi(matIdx, g)

      self % source(idx) = chi * fission + scatter
      if (allocated(self % fixedSource)) then
        self % source(idx) = self % source(idx) + self % fixedSource(idx)
      end if

    end do

  end subroutine sourceUpdateKernelFlatIso
  
  !!
  !! Kernel to update sources given a cell index for linear sources
  !! with isotropic scattering
  !!
  subroutine sourceUpdateKernelLinIso(self, cIdx, ONE_KEFF, it)
    class(arraysRR), target, intent(inout)  :: self
    integer(shortInt), intent(in)           :: cIdx
    real(defFlt), intent(in)                :: ONE_KEFF
    integer(shortInt), intent(in)           :: it
    real(defFlt)                            :: scatter, xScatter, yScatter, zScatter, &
                                               fission, xFission, yFission, zFission, &
                                               xSource, ySource, zSource, temp
    type(particle)                          :: s
    real(defFlt), dimension(matSize)        :: invM
    real(defFlt), dimension(self % nG)      :: fluxFlt, xFlt, yFlt, zFlt
    real(defFlt), dimension(self % nG)              :: nuFission, chi, scatterXS, scatterVec
    integer(shortInt)                       :: matIdx, g, gIn, baseIdx, idx, sIdx1, sIdx2
    real(defReal)                           :: pos(3), dir(3)

    ! Invert moment matrix
    invM = self % invertMatrix(cIdx)
    
    ! Identify material
    matIdx = self % geom % geom % graph % getMatFromUID(cIdx)
    
    ! Guard against void cells
    if (matIdx > self % XSData % getNMat()) then
      baseIdx = self % nG * (cIdx - 1)
      do g = 1, self % nG
        idx = baseIdx + g
        self % source(idx) = 0.0_defFlt
        if (allocated(self % fixedSource)) then
          self % source(idx) = self % source(idx) + self % fixedSource(idx)
        end if
        self % sourceX(idx) = 0.0_defFlt
        self % sourceY(idx) = 0.0_defFlt
        self % sourceZ(idx) = 0.0_defFlt
      end do
      return
    end if

    ! Fudge a particle state to search tally map
    pos = self % cellPos(:,cIdx)
    dir = [0.5,0.5,0.5]
    call s % coords % init(pos, dir)
    s % coords % matIdx = matIdx
    if (self % feedback) then
      if (it == 1) then
        temp = self % XSData % getTemperature(matIdx,2)
      else
        temp = self % tempField % at(s % coords)
      end if
    end if
     
    ! Calculate fission source
    fission = 0.0_defFlt
    if (self % feedback) then
      nuFission = self % XSData % getNuFissXS(matIdx, temp)
    else
      nuFission = self % XSData % getNuFissXS(matIdx)
    end if

    ! Obtain XSs
    !call self % XSData % getProdPointers(matIdx, nuFission, scatterXS, chi)

    baseIdx = self % nG * (cIdx - 1)
    fluxFlt  = real(self % prevFlux((baseIdx + 1):(baseIdx + self % nG)), defFlt)
    xFlt = real(self % prevX((baseIdx + 1):(baseIdx + self % nG)), defFlt)
    yFlt = real(self % prevY((baseIdx + 1):(baseIdx + self % nG)), defFlt)
    zFlt = real(self % prevZ((baseIdx + 1):(baseIdx + self % nG)), defFlt)
    
    ! Calculate fission source
    fission = 0.0_defFlt
    xFission = 0.0_defFlt
    yFission = 0.0_defFlt
    zFission = 0.0_defFlt

    !$omp simd reduction(+:fission, xFission, yFission, zFission)
    do gIn = 1, self % nG
      fission = fission + fluxFlt(gIn) * nuFission(gIn)
      xFission = xFission + xFlt(gIn) * nuFission(gIn)
      yFission = yFission + yFlt(gIn) * nuFission(gIn)
      zFission = zFission + zFlt(gIn) * nuFission(gIn)
    end do
    fission = fission * ONE_KEFF
    xFission = xFission * ONE_KEFF
    yFission = yFission * ONE_KEFF
    zFission = zFission * ONE_KEFF

    !call self % XSData % getTotalPointer(matIdx, total)
    do g = 1, self % nG

      sIdx1 = self % nG * (g - 1) + 1
      sIdx2 = self % nG * g
      !scatterVec => scatterXS(sIdx1:sIdx2)

      if (self % feedback) then
        scatterVec = self % XSData % getScatterXS(matIdx, g, temp)
      else
        scatterVec = self % XSData % getScatterXS(matIdx, g)
      end if

      chi = self % XSData % getChi(matIdx, g)

      ! Calculate scattering source
      scatter = 0.0_defFlt
      xScatter = 0.0_defFlt
      yScatter = 0.0_defFlt
      zScatter = 0.0_defFlt
      !$omp simd reduction(+:scatter, xScatter, yScatter, zScatter)
      do gIn = 1, self % nG
        scatter = scatter + fluxFlt(gIn) * scatterVec(gIn)
        xScatter = xScatter + xFlt(gIn) * scatterVec(gIn)
        yScatter = yScatter + yFlt(gIn) * scatterVec(gIn)
        zScatter = zScatter + zFlt(gIn) * scatterVec(gIn)
      end do

      ! Output index
      idx = baseIdx + g

      self % source(idx) = chi(g) * fission + scatter
      if (allocated(self % fixedSource)) then
        self % source(idx) = self % source(idx) + self % fixedSource(idx)
      end if
      
      xSource = chi(g) * xFission + xScatter
      ySource = chi(g) * yFission + yScatter
      zSource = chi(g) * zFission + zScatter
      
      if (allocated(self % fixedX)) then
        xSource = xSource + self % fixedX(idx)
        ySource = ySource + self % fixedY(idx)
        zSource = zSource + self % fixedZ(idx)
      end if
      
      ! Calculate source gradients by inverting the moment matrix
      self % sourceX(idx) = invM(xx) * xSource + &
              invM(xy) * ySource + invM(xz) * zSource
      self % sourceY(idx) = invM(xy) * xSource + &
              invM(yy) * ySource + invM(yz) * zSource
      self % sourceZ(idx) = invM(xz) * xSource + &
           invM(yz) * ySource + invM(zz) * zSource
      
    end do

  end subroutine sourceUpdateKernelLinIso
  
  !!
  !! Kernel to update sources given a cell index for flat sources
  !! with anisotropic scattering
  !!
  subroutine sourceUpdateKernelFlatAni(self, cIdx, ONE_KEFF)
    class(arraysRR), target, intent(inout)   :: self
    integer(shortInt), intent(in)            :: cIdx
    real(defFlt), intent(in)                 :: ONE_KEFF

  end subroutine sourceUpdateKernelFlatAni
  
  !!
  !! Kernel to update sources given a cell index for linear sources
  !! with flat anisotropic scattering
  !!
  subroutine sourceUpdateKernelLIFA(self, cIdx, ONE_KEFF)
    class(arraysRR), target, intent(inout)   :: self
    integer(shortInt), intent(in)            :: cIdx
    real(defFlt), intent(in)                 :: ONE_KEFF

  end subroutine sourceUpdateKernelLIFA
  
  !!
  !! Inverts the spatial moment matrix for use in linear source calculations.
  !!
  function invertMatrix(self, cIdx) result(invM)
    class(arraysRR), target, intent(in) :: self
    integer(shortInt), intent(in)       :: cIdx
    real(defFlt), dimension(matSize)    :: invM
    integer(shortInt)                   :: condX, condY, condZ, inversionTest
    real(defReal)                       :: det
    real(defFlt)                        :: one_det

    associate(momVec => self % momMat(((cIdx - 1) * matSize + 1):(cIdx * matSize)))

    ! Pre-invert the moment matrix
    ! Need to check for poor conditioning by evaluating the
    ! diagonal elements of the matrix
    condX = 0
    condY = 0
    condZ = 0

    ! Trying out simpler matrix test
    if (momVec(xx) > condition_tolerance) condX = 1
    if (momVec(yy) > condition_tolerance) condY = 1
    if (momVec(zz) > condition_tolerance) condZ = 1
    ! Significantly stabilises 2D linear source problems. Z moments can vary wildly.
    if (self % set2D) condZ = 0

    ! Map conditions to test variable
    inversionTest = condX * 4 + condY * 2 + condZ
    invM = 0.0_defFlt

    select case(inversionTest)
    case(invertXYZ)
      det = momVec(xx) * (momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)) &
            - momVec(yy) * momVec(xz) * momVec(xz) &
            - momVec(zz) * momVec(xy) * momVec(xy) &
            + 2 * momVec(xy) * momVec(xz) * momVec(yz)
      invM(xx) = real(momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz),defFlt)
      invM(xy) = real(momVec(xz) * momVec(yz) - momVec(xy) * momVec(zz),defFlt)
      invM(xz) = real(momVec(xy) * momVec(yz) - momVec(yy) * momVec(xz),defFlt)
      invM(yy) = real(momVec(xx) * momVec(zz) - momVec(xz) * momVec(xz),defFlt)
      invM(yz) = real(momVec(xy) * momVec(xz) - momVec(xx) * momVec(yz),defFlt)
      invM(zz) = real(momVec(xx) * momVec(yy) - momVec(xy) * momVec(xy),defFlt)

    case(invertYZ)
      det = momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)
      invM(yy) = real(momVec(zz),defFlt)
      invM(yz) = real(-momVec(yz),defFlt)
      invM(zz) = real(momVec(yy),defFlt)

    case(invertXY)
      det = momVec(xx) * momVec(yy) - momVec(xy) * momVec(xy)
      invM(xx) = real(momVec(yy),defFlt)
      invM(xy) = real(-momVec(xy),defFlt)
      invM(yy) = real(momVec(xx),defFlt)

    case(invertXZ)
      det = momVec(xx) * momVec(zz) - momVec(xz) * momVec(xz)
      invM(xx) = real(momVec(zz),defFlt)
      invM(xz) = real(-momVec(xz),defFlt)
      invM(zz) = real(momVec(xx),defFlt)

    case(invertX)
      det = momVec(xx)
      invM(xx) = 1.0_defFlt

    case(invertY)
      det = momVec(yy)
      invM(yy) = 1.0_defFlt

    case(invertZ)
      det = momVec(zz)
      invM(zz) = 1.0_defFlt

    case default
      det = ONE
    end select

    one_det = real(ONE/det, defFlt)
    invM = invM * one_det

    ! Check for zero determinant
    if (abs(det) < det_tolerance) invM = 0.0_defFlt

    end associate

  end function invertMatrix
  
  !! 
  !! Calculate keff
  !! Wraps the main kernel call to allow for OMP + SIMD (thanks Fortran)
  !!
  function calculateKeff(self, k0) result(k1)
    class(arraysRR), intent(in)           :: self
    real(defReal), intent(in)             :: k0
    real(defReal)                         :: k1
    integer(shortInt)                     :: cIdx
    real(defReal)                         :: fissTotal, prevFissTotal
    real(defReal), save                   :: fissLocal, prevFissLocal
    character(100), parameter             :: Here = 'calculateKeff (arraysRR_class.f90)'
    !$omp threadprivate (fissLocal, prevFissLocal)

    fissTotal     = ZERO
    prevFissTotal = ZERO
    !$omp parallel do reduction(+:fissTotal, prevFissTotal)
    do cIdx = 1, self % nCells
      call self % calculateKeffKernel(cIdx, fissLocal, prevFissLocal)
      fissTotal     = fissTotal + fissLocal
      prevFissTotal = prevFissTotal + prevFissLocal
    end do 
    !$omp end parallel do

    !print *, fissTotal
    !print *, prevFissTotal
    
    k1 = k0 * fissTotal / prevFissTotal 
    if ((k1 <= 0) .or. (k1 > 5)) call fatalError(Here, 'Unphysical keff: '//numToChar(k1))
    if (k1 /= k1) call fatalError(Here, 'NaN keff')

  end function calculateKeff
  
  !!
  !! Calculate keff for a single cell
  !!
  subroutine calculateKeffKernel(self, cIdx, fissionRate, prevFissionRate)
    class(arraysRR), target, intent(in)  :: self
    integer(shortInt), intent(in)        :: cIdx
    real(defReal), intent(out)           :: fissionRate
    real(defReal), intent(out)           :: prevFissionRate
    real(defReal)                        :: vol
    type(particle)                       :: s
    integer(shortInt)                    :: g, matIdx
    real(defFlt)                         :: temp
    real(defFlt), dimension(self % nG)   :: nuSigmaF
    real(defFlt), dimension(:), pointer  :: flux
    real(defFlt), dimension(:), pointer  :: prevFlux
    real(defReal)                        :: pos(3), dir(3)

    vol = self % volume(cIdx)
    fissionRate     = ZERO
    prevFissionRate = ZERO

    if (vol < volume_tolerance) return

    ! Identify material
    matIdx = self % geom % geom % graph % getMatFromUID(cIdx)

    ! Check whether to continue in this cell
    if (matIdx > self % XSData % getNMat()) return
    if (.not. self % XSData % isFissile(matIdx)) return

    vol = self % volume(cIdx)

    if (vol < volume_tolerance) return

    if ((matIdx < 1) .or. (matIdx > self % XSData % getNMat())) then
      print *, "matIdx = " // numToChar(matIdx) // " in calculateKeff"
    end if

    ! Fudge a particle state to search tally map
    pos = self % cellPos(:,cIdx)
    dir = [0.5,0.5,0.5]

    call s % coords % init(pos, dir)
    s % coords % matIdx = matIdx

    if (self % feedback) then

      temp = self % tempField % at(s % coords)

      nuSigmaF = self % XSData % getNuFissXS(matIdx, temp)

    else

      nuSigmaF = self % XSData % getNuFissXS(matIdx)

    end if

    flux => self % scalarFlux((self % nG * (cIdx - 1) + 1):(self % nG * cIdx))

    prevFlux => self % prevFlux((self % nG * (cIdx - 1) + 1):(self % nG * cIdx))

    !$omp simd reduction(+:fissionRate, prevFissionRate)
    do g = 1, self % nG

      fissionRate = fissionRate &
        + real(flux(g) * nuSigmaF(g), defReal)

      prevFissionRate = prevFissionRate &
        + real(prevFlux(g) * nuSigmaF(g), defReal)

    end do

    fissionRate     = fissionRate * vol
    prevFissionRate = prevFissionRate * vol

  end subroutine calculateKeffKernel 

  !!
  !! Zero the previous-step flux
  !!
  subroutine zeroPrevFlux(self)
    class(arraysRR), intent(inout) :: self
    character(100), parameter      :: Here = 'zeroPrevFlux (arraysRR_class.f90)'

    if (allocated(self % prevFlux)) then
      self % prevFlux = ZERO
    else
      call fatalError(Here,'prevFlux has not been initialised')
    end if

  end subroutine zeroPrevFlux

  !!
  !! Reset fluxes
  !!
  subroutine resetFluxes(self)
    class(arraysRR), intent(inout) :: self
    character(100), parameter      :: Here = 'resetFluxes (arraysRR_class.f90)'

    select case(self % simulationType)
      case(flatIso)
        call self % resetFluxesFlatIso()
      case(linearIso)
        call self % resetFluxesLinIso()
      case default
        call fatalError(Here,'Unsupported simulation type requested')
    end select

  end subroutine resetFluxes
  
  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !!
  subroutine resetFluxesFlatIso(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: idx

    !$omp parallel do 
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = 0.0_defFlt
    end do
    !$omp end parallel do

  end subroutine resetFluxesFlatIso
  
  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !! for linear sources with isotropic scattering
  !!
  subroutine resetFluxesLinIso(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: idx

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = 0.0_defFlt
    end do
    !$omp end parallel do

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevX(idx) = self % scalarX(idx)
      self % scalarX(idx) = ZERO
    end do
    !$omp end parallel do

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevY(idx) = self % scalarY(idx)
      self % scalarY(idx) = ZERO
    end do
    !$omp end parallel do

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevZ(idx) = self % scalarZ(idx)
      self % scalarZ(idx) = ZERO
    end do
    !$omp end parallel do

  end subroutine resetFluxesLinIso

  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !! for flat sources with anisotropic scattering
  !!
  subroutine resetFluxesFlatAni(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: idx

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = 0.0_defFlt
    end do
    !$omp end parallel do

  end subroutine resetFluxesFlatAni
  
  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !! for linear sources with flat anisotropic scattering
  !!
  subroutine resetFluxesLIFA(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: idx

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = 0.0_defFlt
    end do
    !$omp end parallel do

  end subroutine resetFluxesLIFA


  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxScores(self)
    class(arraysRR), intent(inout) :: self
    real(defReal), save            :: flux, pow
    integer(shortInt)              :: idx
    character(100), parameter      :: Here = 'accumulateFluxScores (arraysRR_class.f90)'
    !$omp threadprivate(flux, pow)

    select case(self % simulationType)
      case(flatIso, flatAni)
        call self % accumulateFluxScoresFlat()
      case(linearIso, linearAni)
        call self % accumulateFluxScoresLinear()
      case default
        call fatalError(Here,'Unsupported simulation type requested')
    end select

    !$omp parallel do schedule(static)
    do idx = 1, size(self % power)
      pow = real(self % power(idx),defReal)
      !print *, pow
      self % powerScores(1, idx) = self % powerScores(1, idx) + pow
      self % powerScores(2, idx) = self % powerScores(2, idx) + pow * pow
    end do
    !$omp end parallel do

    !print *, self % powerScores(1,:)

  end subroutine accumulateFluxScores

  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxScoresFlat(self)
    class(arraysRR), intent(inout) :: self
    real(defReal), save            :: flux
    integer(shortInt)              :: idx
    !$omp threadprivate(flux)

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      flux = self % scalarFlux(idx)
      self % fluxScores(1, idx) = self % fluxScores(1, idx) + flux
      self % fluxScores(2, idx) = self % fluxScores(2, idx) + flux * flux
    end do
    !$omp end parallel do

  end subroutine accumulateFluxScoresFlat
 
  !!
  !! Accumulate flux scores for stats
  !! Includes linear components
  !!
  subroutine accumulateFluxScoresLinear(self)
    class(arraysRR), intent(inout) :: self
    real(defReal), save            :: flux
    integer(shortInt)              :: idx
    !$omp threadprivate(flux)

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      flux = self % scalarFlux(idx)
      self % fluxScores(1, idx) = self % fluxScores(1, idx) + flux
      self % fluxScores(2, idx) = self % fluxScores(2, idx) + flux * flux
    end do
    !$omp end parallel do
    
    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      flux = self % scalarX(idx)
      self % xScores(1, idx) = self % xScores(1, idx) + flux
      self % xScores(2, idx) = self % xScores(2, idx) + flux * flux
    end do
    !$omp end parallel do
    
    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      flux = self % scalarY(idx)
      self % yScores(1, idx) = self % yScores(1, idx) + flux
      self % yScores(2, idx) = self % yScores(2, idx) + flux * flux
    end do
    !$omp end parallel do
    
    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      flux = self % scalarZ(idx)
      self % zScores(1, idx) = self % zScores(1, idx) + flux
      self % zScores(2, idx) = self % zScores(2, idx) + flux * flux
    end do
    !$omp end parallel do

  end subroutine accumulateFluxScoresLinear

  !!
  !! Tallies flux-related results, using an input tallyAdmin.
  !! Allows for use of MC tally machinery with RR, although limited
  !! to relatively simple estimators.
  !!
  !! Loops over all phase space points, contributing each to tallies.
  !!
  !! Assumes the cycle will be ended afterwards.
  !!
  subroutine tallyResults(self, tally)
    class(arraysRR), intent(in)                :: self
    type(tallyAdmin), pointer, intent(inout)   :: tally
    type(particle), save                       :: p
    real(defReal), save                        :: vol
    real(defReal), dimension(3), save          :: pos
    integer(shortInt), save                    :: g, matIdx
    real(defFlt), dimension(:), pointer, save :: fluxVec
    integer(longInt)                              :: seed
    integer(shortInt)                             :: seed_temp
    character(10)                                 :: time
    character(8)                                  :: date
    character(:),allocatable                      :: string
    type(RNG)                                     :: rand
    type(RNG), target, save                       :: pRNG
    integer(shortInt)                          :: i
    !$omp threadprivate(p, pRNG, vol, pos, g, matIdx, fluxVec)

    ! Initialise RNG
    call date_and_time(date, time)
    string = date // time
    call FNV_1(string,seed_temp)
    seed = seed_temp
    call rand % init(seed)

    !$omp parallel
    call p % build([-INFINITY, -INFINITY, -INFINITY], [ONE, ZERO, ZERO], 1, ZERO, ZERO)
    !$omp end parallel

    !$omp parallel do
    do i = 1, self % nCells

      vol = self % getVolume(i)
      pos = self % getCellPos(i)
      call p % teleport(pos)

      call self % getFluxPointer(i, fluxVec)
      matIdx = self % geom % geom % graph % getMatFromUID(i) 
      p % coords % matIdx = matIdx 

      pRNG = rand
      call pRNG % stride(i)
      p % pRNG => pRNG

      do g = 1, self % nG

        ! The weight should be flux * V, which, for a single score in a tally
        ! will produce SigmaX * flux * V or volume-integrated reaction rate.
        p % w = fluxVec(g) * vol

        p % G = g

        call tally % reportInColl(p, .false.)

      end do

    end do
    !$omp end parallel do
    
  end subroutine tallyResults
  
  !!
  !! Finalise flux scores for stats
  !!
  subroutine finaliseFluxScores(self,it)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: it
    integer(shortInt)              :: idx
    real(defReal)                  :: N1, Nm1
    character(100), parameter      :: Here = 'finaliseFluxScores (arraysRR_class.f90)'

    select case(self % simulationType)
      case(flatIso, flatAni)
        call self % finaliseFluxScoresFlat(it)
      case(linearIso, linearAni)
        call self % finaliseFluxScoresLinear(it)
      case default
        call fatalError(Here,'Unsupported simulation type requested')
    end select

    if (it /= 1) then
      Nm1 = 1.0_defReal/(it - 1)
    else
      Nm1 = 1.0_defReal
    end if
    N1 = 1.0_defReal/it

    !Finalise powers
    !$omp parallel do schedule(static)
    do idx = 1, size(self % power)
      self % powerScores(1, idx) = self % powerScores(1, idx) * N1
      self % powerScores(2, idx) = self % powerScores(2, idx) * N1
      self % powerScores(2, idx) = Nm1 * (self % powerScores(2, idx) - &
            self % powerScores(1, idx) * self % powerScores(1, idx))
      if (self % powerScores(2, idx) <= ZERO) then
        self % powerScores(2, idx) = ZERO
      else
        self % powerScores(2, idx) = sqrt(self % powerScores(2, idx))
      end if
    end do
    !$omp end parallel do

  end subroutine finaliseFluxScores

  !!
  !! Finalise flux scores for stats
  !!
  subroutine finaliseFluxScoresFlat(self,it)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: it
    integer(shortInt)              :: idx
    real(defReal)                  :: N1, Nm1

    if (it /= 1) then
      Nm1 = 1.0_defReal/(it - 1)
    else
      Nm1 = 1.0_defReal
    end if
    N1 = 1.0_defReal/it

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % fluxScores(1, idx) = self % fluxScores(1, idx) * N1
      self % fluxScores(2, idx) = self % fluxScores(2, idx) * N1
      self % fluxScores(2, idx) = Nm1 * (self % fluxScores(2, idx) - &
            self % fluxScores(1, idx) * self % fluxScores(1, idx))
      if (self % fluxScores(2, idx) <= ZERO) then
        self % fluxScores(2, idx) = ZERO
      else
        self % fluxScores(2, idx) = sqrt(self % fluxScores(2, idx))
      end if
    end do
    !$omp end parallel do

  end subroutine finaliseFluxScoresFlat

  !!
  !! Finalise flux scores for stats
  !! Includes linear components
  !!
  subroutine finaliseFluxScoresLinear(self,it)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: it
    integer(shortInt)              :: idx
    real(defReal)                  :: N1, Nm1

    if (it /= 1) then
      Nm1 = 1.0_defReal/(it - 1)
    else
      Nm1 = 1.0_defReal
    end if
    N1 = 1.0_defReal/it

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % fluxScores(1, idx) = self % fluxScores(1, idx) * N1
      self % fluxScores(2, idx) = self % fluxScores(2, idx) * N1
      self % fluxScores(2, idx) = Nm1 * (self % fluxScores(2, idx) - &
            self % fluxScores(1, idx) * self % fluxScores(1, idx))
      if (self % fluxScores(2, idx) <= ZERO) then
        self % fluxScores(2, idx) = ZERO
      else
        self % fluxScores(2, idx) = sqrt(self % fluxScores(2, idx))
      end if

      self % xScores(1, idx) = self % xScores(1, idx) * N1
      self % xScores(2, idx) = self % xScores(2, idx) * N1
      self % xScores(2, idx) = Nm1 * (self % xScores(2, idx) - &
            self % xScores(1, idx) * self % xScores(1, idx))
      if (self % xScores(2, idx) <= ZERO) then
        self % xScores(2, idx) = ZERO
      else
        self % xScores(2, idx) = sqrt(self % xScores(2, idx))
      end if

      self % yScores(1, idx) = self % yScores(1, idx) * N1
      self % yScores(2, idx) = self % yScores(2, idx) * N1
      self % yScores(2, idx) = Nm1 * (self % yScores(2, idx) - &
            self % yScores(1, idx) * self % yScores(1, idx))
      if (self % yScores(2, idx) <= ZERO) then
        self % yScores(2, idx) = ZERO
      else
        self % yScores(2, idx) = sqrt(self % yScores(2, idx))
      end if

      self % zScores(1, idx) = self % zScores(1, idx) * N1
      self % zScores(2, idx) = self % zScores(2, idx) * N1
      self % zScores(2, idx) = Nm1 * (self % zScores(2, idx) - &
            self % zScores(1, idx) * self % zScores(1, idx))
      if (self % zScores(2, idx) <= ZERO) then
        self % zScores(2, idx) = ZERO
      else
        self % zScores(2, idx) = sqrt(self % zScores(2, idx))
      end if
    end do
    !$omp end parallel do

  end subroutine finaliseFluxScoresLinear
  
  !!
  !! Outputs integral flux or fission rate when
  !! given a tally map
  !!
  !subroutine outputMap(self, out, map, doFission)
  subroutine outputMap(self, out, map, option)

    class(arraysRR), intent(in)                :: self
    class(outputFile), intent(inout)           :: out
    class(tallyMap), pointer, intent(in)       :: map

    integer(shortInt), intent(in)              :: option

    character(nameLen)                         :: name

    integer(shortInt)                          :: cIdx

    integer(shortInt), dimension(:), allocatable :: resArrayShape

    type(particleState), save                  :: p
    type(particle), save                       :: s

    real(defReal), save                        :: vol

    real(defFlt), save                         :: sig, temp

    real(defFlt), dimension(self % nG)         :: fissXS

    integer(shortInt), save                    :: i, matIdx, g

    real(defReal), dimension(:), allocatable   :: res, resSD

    real(defReal), dimension(3), save          :: pos, dir

    !$omp threadprivate(temp, s, vol, sig, i, matIdx, g, p)

    resArrayShape = [map % binArrayShape()]

    allocate(res(map % bins(0)))
    allocate(resSD(map % bins(0)))

    res   = ZERO
    resSD = ZERO

    ! Find whether cells are in map and sum contributions
    !$omp parallel do reduction(+:res,resSD)
    do cIdx = 1, self % nCells

      vol = self % volume(cIdx)

      if (vol < volume_tolerance) cycle

      matIdx = self % geom % geom % graph % getMatFromUID(cIdx)

      ! Fudge particle state to search tally map
      p % r = self % cellPos(:,cIdx)

      i = map % map(p)

      ! Fudge particle to search temperature map
      pos = self % cellPos(:,cIdx)
      dir = [0.5,0.5,0.5]

      call s % coords % init(pos, dir)

      s % coords % matIdx = matIdx

      if (self % feedback) then

        temp = self % tempField % at(s % coords)

      end if

      if (i > 0) then

        matIdx = self % geom % geom % graph % getMatFromUID(cIdx)

        if (option == 3) then

          res(i) = res(i) + self % getPowerScore(cIdx)

          resSD(i) = resSD(i) + self % getPowerSD(cIdx)**2

        else

          ! Fetch fission XS vector once per cell
          if (option == 1) then

            if (self % feedback) then
              fissXS = self % XSData % getFissionXS(matIdx, temp)
            else
              fissXS = self % XSData % getFissionXS(matIdx)
            end if

          end if

          do g = 1, self % nG

            if (option == 1) then
              sig = fissXS(g)
            else
              sig = 1.0_defFlt
            end if

            res(i) = res(i) + vol * self % getFluxScore(cIdx,g) * sig

            ! Neglect uncertainty in volume
            resSD(i) = resSD(i) + &
              self % getFluxSD(cIdx,g)**2 * vol * vol * sig * sig

          end do

        end if

      end if

    end do
    !$omp end parallel do

    do i = 1, size(resSD)

      resSD(i) = sqrt(resSD(i))

      if (res(i) > 0) then
        resSD(i) = resSD(i) / res(i)
      end if

    end do

    if (option == 1) then
      name = 'fiss1G'
    else if (option == 2) then
      name = 'flux1G'
    else
      name = 'power'
    end if

    call out % startBlock(name)

    call out % startArray(name, resArrayShape)

    do i = 1, map % bins(0)
      call out % addResult(res(i), resSD(i))
    end do

    call out % endArray()

    call map % print(out)

    call out % endBlock()

    deallocate(res)
    deallocate(resSD)

  end subroutine outputMap
  
  !!
  !! Send all arrays of interest to VTK output
  !!
  subroutine outputToVTK(self, viz)
    class(arraysRR), intent(in)               :: self
    class(visualiser), intent(inout)          :: viz
    real(defReal), dimension(:), allocatable  :: resVec
    character(nameLen)                        :: name
    integer(shortInt)                         :: cIdx, g

    allocate(resVec(self % nCells))

    ! Output all fluxes (assuming finalisation of scores happened)
    do g = 1, self % nG
      name = 'flux_g'//numToChar(g)
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells
        resVec(cIdx) = self % getFluxScore(cIdx,g)
      end do
      !$omp end parallel do
      call viz % addVTKData(resVec,name)
    end do

    ! Output all flux uncertainties
    do g = 1, self % nG
      name = 'std_g'//numToChar(g)
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells
        resVec(cIdx) = self % getFluxSD(cIdx,g) /self % getFluxScore(cIdx,g)
      end do
      !$omp end parallel do
      call viz % addVTKData(resVec,name)
    end do

    ! Output final iteration sources
    do g = 1, self % nG
      name = 'source_'//numToChar(g)
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells
        resVec(cIdx) = real(self % getSource(cIdx,g),defReal)
      end do
      !$omp end parallel do
      call viz % addVTKData(resVec,name)
    end do

    ! Output final volume estimates
    ! TODO: scale to be absolute, not relative
    name = 'volume'
    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      resVec(cIdx) = self % volume(cIdx)
    end do
    !$omp end parallel do
    call viz % addVTKData(resVec,name)

    call viz % finaliseVTK

  end subroutine outputToVTK

  !!
  !! Returns the flux vector at a point in space
  !!
  subroutine getFluxAtAPoint(self, r, flux, fluxSD)
    class(arraysRR), intent(in)                      :: self
    real(defReal), dimension(3), intent(in)          :: r
    real(defReal), dimension(self % nG), intent(out) :: flux
    real(defReal), dimension(self % nG), intent(out) :: fluxSD
    integer(shortInt)                                :: g, matIdx, cIdx, i
    real(defReal), dimension(3)                      :: mom, momSD, centroid, fluxGrad
    real(defFlt), dimension(matSize)                 :: invM

    ! Identify cell at the given point
    call self % geom % whatIsAt(matIdx, cIdx, r)

    if (cIdx > 0) then
      do g = 1, self % nG

        flux(g) = self % getFluxScore(cIdx, g)
        fluxSD(g) = self % getFluxSD(cIdx, g)

        ! Include linear moments if available
        if ((self % simulationType == linearIso) .or. &
                (self % simulationType == linearAni)) then

          fluxSD(g) = fluxSD(g) * fluxSD(g)

          mom = self % getFluxMoments(cIdx, g)
          momSD = self % getFluxMomentSDs(cIdx, g)
          centroid = self % getCentroid(cIdx)
          invM = self % invertMatrix(cIdx)

          ! Get flux gradients
          fluxGrad(x) = real(invM(xx) * mom(x) + invM(xy) * mom(y) + invM(xz) * mom(z), defReal)
          fluxGrad(y) = real(invM(xy) * mom(x) + invM(yy) * mom(y) + invM(yz) * mom(z), defReal)
          fluxGrad(z) = real(invM(xz) * mom(x) + invM(yz) * mom(y) + invM(zz) * mom(z), defReal)

          ! Note this will not correctly estimate uncertainty as moments are covariant
          do i = 1, 3
            flux(g) = flux(g) + fluxGrad(i) * (r(i) - centroid(i))
            ! Not sure exactly how to propagate uncertainties - need to do some maths
            fluxSD(g) = fluxSD(g) + momSD(i)**2 * (r(i) - centroid(i))**2
          end do

          if (fluxSD(g) > ZERO) fluxSD(g) = sqrt(fluxSD(g))

        end if

      end do

    else
      print *,'WARNING: No cell found at position '//numToChar(r)
      flux = -ONE
      fluxSD = -ONE
    end if

  end subroutine getFluxAtAPoint
  
  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt) :: i

    ! Clean contents
    if(allocated(self % scalarFlux)) deallocate(self % scalarFlux)
    if(allocated(self % prevFlux)) deallocate(self % prevFlux)
    if(allocated(self % fluxScores)) deallocate(self % fluxScores)
    if(allocated(self % power)) deallocate(self % power)
    if(allocated(self % powerScores)) deallocate(self % powerScores)
    if(allocated(self % source)) deallocate(self % source)
    if(allocated(self % fixedSource)) deallocate(self % fixedSource)
    if(allocated(self % sourceIdx)) deallocate(self % sourceIdx)
    if(allocated(self % volumeTracks)) deallocate(self % volumeTracks)
    if(allocated(self % lengthSquared)) deallocate(self % lengthSquared)
    if(allocated(self % volume)) deallocate(self % volume)
    if(allocated(self % cellHit)) deallocate(self % cellHit)
    if(allocated(self % cellTotalHit)) deallocate(self % cellTotalHit)
    if(allocated(self % cellFound)) deallocate(self % cellFound)
    if(allocated(self % cellPos)) deallocate(self % cellPos)
    
    if(allocated(self % locks)) then
      do i = 1, self % nCells
#ifdef _OPENMP
        call OMP_destroy_lock(self % locks(i))
#endif
      end do
      deallocate(self % locks)
    end if
    
    self % geom   => null()
    call self % XSData % kill()
    self % nG     = 0
    self % nCells = 0
    self % lengthPerIt = ZERO
    self % rho         = 0.0_defFlt
    self % simulationType = 0
    self % ani = 0

  end subroutine kill

end module arraysRR_class
