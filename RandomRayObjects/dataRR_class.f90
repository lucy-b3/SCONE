module dataRR_class

  use numPrecision
  use universalVariables
  use rng_class,                      only : RNG
  use genericProcedures,              only : fatalError

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat => nMat, mm_matName => matName
  use materialHandle_inter,           only : materialHandle
  use baseMgNeutronDatabase_class,    only : baseMgNeutronDatabase
  use baseMgNeutronMaterial_class,    only : baseMgNeutronMaterial, baseMgNeutronMaterial_CptrCast
  
  implicit none
  private
  
  !! XS storage for each material (with or without TH feedback)
  type :: matXS
    real(defFlt), allocatable :: sigmaT(:)    ! (g, t)
    real(defFlt), allocatable :: nuSigmaF(:)  ! (g, t)
    real(defFlt), allocatable :: sigmaF(:)    ! (g, t)
    real(defFlt), allocatable :: chi(:)         ! (g)
    real(defFlt), allocatable :: sigmaS(:)  ! (g1, g, t)
    !logical :: fissile
  end type matXS

  !!
  !! Nuclear data in a random ray-friendly format.
  !!
  !! Stores data and provides access in a manner which is more performant
  !! than is done for MC MG data at present.
  !!
  !! TODO: Add kinetic data and higher-order scattering matrices
  !!
  type, public :: dataRR
    private
    ! Components
    integer(shortInt)                     :: nG     = 0
    integer(shortInt)                     :: nG2    = 0
    integer(shortInt)                     :: nMat   = 0
    integer(shortInt)                     :: nT     = 0
    ! real(defReal), dimension(:), allocatable :: T

    ! Data space - absorb all nuclear data for speed
    !real(defFlt), dimension(:), allocatable       :: sigmaT
    !real(defFlt), dimension(:,:), allocatable       :: sigmaT
    !real(defFlt), dimension(:,:), allocatable       :: nuSigmaF
    !real(defFlt), dimension(:,:), allocatable       :: sigmaF
    !real(defFlt), dimension(:,:), allocatable       :: sigmaS
    !real(defFlt), dimension(:,:), allocatable       :: chi
    logical(defBool), dimension(:), allocatable   :: fissile
    type(matXS), dimension(:), allocatable         :: matRR ! 
    character(nameLen), dimension(:), allocatable :: names

    ! Optional kinetic data
    logical(defBool)                        :: doKinetics  = .false.
    integer(shortInt)                       :: nP = 0
    real(defFlt), dimension(:), allocatable :: chiD
    real(defFlt), dimension(:), allocatable :: chiP
    real(defFlt), dimension(:), allocatable :: beta
    real(defFlt), dimension(:), allocatable :: invSpeed

    ! Optional higher-order scattering matrices up to P3
    real(defFlt), dimension(:), allocatable :: sigmaS1
    real(defFlt), dimension(:), allocatable :: sigmaS2
    real(defFlt), dimension(:), allocatable :: sigmaS3

    ! Stores the temperatures of the XS files
    real(defFlt), dimension(:,:), allocatable :: temperatures

  contains
    
    procedure :: init
    procedure :: setAdjointXS
    procedure :: kill

    ! TODO: add an XS update procedure, e.g., given multiphysics
    ! TODO: add full handling of kinetic data
    ! TODO: add full handling of higher-order anisotropy

    ! Access procedures
    procedure :: getProdPointers
    !procedure :: getAllPointers
    procedure :: getTotalPointer
    procedure :: getNuFissPointer
    procedure :: getChiPointer
    procedure :: getScatterPointer
    procedure :: getScatterVecPointer
    procedure :: getTotalXS
    procedure :: getNuFissXS
    procedure :: getFissionXS
    procedure :: getScatterXS
    procedure :: getChi
    procedure :: getNG
    procedure :: getNT
    procedure :: getNMat
    procedure :: getNPrec
    procedure :: getName
    procedure :: isFissile
    procedure :: getTemperature

    ! Private procedures
    procedure, private :: getIdxs
    procedure, private :: getScatterIdxs
    !procedure, private :: getKineticIdxs


  end type dataRR

contains

  !!
  !! Initialise necessary nuclear data.
  !! Can optionally include kinetic parameters.
  !!
  subroutine init(self, db, doKinetics, loud, temp)
    class(dataRR), intent(inout)                     :: self
    class(baseMgNeutronDatabase),pointer, intent(in) :: db
    logical(defBool), intent(in)                     :: doKinetics
    !integer(shortInt), intent(in)                    :: aniOrder
    logical(defBool), intent(in)                     :: loud
    logical(defBool), optional, intent(in)           :: temp
    integer(shortInt)                                :: g, g1, m, matP1, t, aniOrder
    type(RNG)                                        :: rand
    logical(defBool)                                 :: fiss
    class(baseMgNeutronMaterial), pointer            :: mat
    class(materialHandle), pointer                   :: matPtr
    character(100), parameter :: Here = 'init (dataRR_class.f90)'

    self % doKinetics = doKinetics

    !print *, "Made it to dataRR"
    ! Store number of energy groups for convenience
    self % nG = db % nGroups()
    self % nG2 = self % nG * self % nG
    !self % nT = db % nTemps()

    ! Initialise local nuclear data
    ! Allocate nMat + 1 materials to catch any undefined materials
    ! TODO: clean nuclear database afterwards! It is no longer used
    !       and takes up memory.
    self % nMat = mm_nMat()
    !print *, "Constants found"
        
    matP1 = self % nMat + 1
    
    allocate(self % matRR(self % nMat))
    allocate(self % fissile(matP1))
    self % fissile = .false.
    allocate(self % names(matP1))
    self % names = 'unnamed'

    !print *, "Allocated matRR and fissle"

    if ((temp)) then

      self % nT = db % nTemps()
      allocate(self % temperatures(self % nMat, self % nT))
      ! Create a dummy RNG togetMacroXSs_byG satisfy the mgDatabase access interface
      call rand % init(1_longInt)
      if (loud) print *,'Initialising random ray nuclear data'
      !print *, 'starting first material'
      do m = 1, self % nMat
        !print *, 'BEGINNING ' // mm_matName(m)
        matPtr  => db % getMaterial(m)
        !print *, "material pointer found"
        mat     => baseMgNeutronMaterial_CptrCast(matPtr)
        !print *, "mat found"
        fiss = .false.

        allocate(self % matRR(m) % sigmaT(self % nG * self % nT))
        !print *, "SigmaT allocated"
        allocate(self % matRR(m) % nuSigmaF(self % nG * self % nT))
        allocate(self % matRR(m) % sigmaF(self % nG * self % nT))
        allocate(self % matRR(m) % chi(self % nG))
        allocate(self % matRR(m) % sigmaS(self % nG * self % nG * self % nT))

        !print *, "Allocated XSs with temps"

        do t = 1, self % nT
          self % temperatures(m,t) = real(mat % getTemp(t),defFlt)
        end do
        !print  *, "mat temps stored"

        !print *, self % temperatures(m,:)

        do g = 1, self % nG
          do t = 1, self % nT
            self % matRR(m) % sigmaT((g - 1) * self % nT + t) = real(mat % getTotalXS(g, rand, t),defFlt)
            !print *, self % matRR(m) % sigmaT((g - 1) * self % nT + t)
            self % matRR(m) % nuSigmaF((g - 1) * self % nT + t) = real(mat % getNuFissionXS(g, rand, t),defFlt)
            !print *, "nuSigmaF stored"
            self % matRR(m) % sigmaF((g - 1) * self % nT + t) = real(mat % getFissionXS(g, rand, t),defFlt)
            !print *, self % matRR(m) % sigmaF((g - 1) * self % nT + t)
            if (self % matRR(m) % nuSigmaF((g - 1) * self % nT + t) > 0) fiss = .true.
            self % matRR(m) % chi(g) = real(mat % getChi(g, rand),defFlt)
            ! Include scattering multiplicity
            do g1 = 1, self % nG
              !self % matRR(m) % sigmaS(((g - 1) * self % nG + (g1 - 1)) * self % nT + t)  = &
                      !real(mat % getScatterXS(g1, g, rand, t) * mat % getProd(g1, g, rand, t) , defFlt)
              self % matRR(m) % sigmaS(((g1 - 1) * self % nG + (g - 1)) * self % nT + t) = &
                      real(mat % getScatterXS(g1, g, rand, t) * mat % getProd(g1, g, rand, t), defFlt)
              !print *, "scattering stored"
              !print *, self % matRR(m) % sigmaS(((g - 1) * self % nG + (g1 - 1)) * self % nT + t)
            end do
          end do
        end do
        self % fissile(m) = fiss
        self % names(m) = mm_matName(m)
        !print *, "MATERIAL FINISHED"
      end do
      

    else

      self % nT = 1
      ! Create a dummy RNG togetMacroXSs_byG satisfy the mgDatabase access interface
      call rand % init(1_longInt)
      if (loud) print *,'Initialising random ray nuclear data'
      !print *, 'no temp found'
      do m = 1, self % nMat
        matPtr  => db % getMaterial(m)
        mat     => baseMgNeutronMaterial_CptrCast(matPtr)
        fiss = .false.

        !print *, "Mat pointers found"
        !print *, mm_matName(m)

        allocate(self % matRR(m) % sigmaT(self % nG))
        allocate(self % matRR(m) % nuSigmaF(self % nG))
        allocate(self % matRR(m) % sigmaF(self % nG))
        allocate(self % matRR(m) % chi(self % nG))
        allocate(self % matRR(m) % sigmaS(self % nG * self % nG))

        !print *, "XS allocations made"

        do g = 1, self % nG
          self % matRR(m) % sigmaT(g) = real(mat % getTotalXS(g, rand),defFlt)
          !print *, "SigmaT stored"
          !print *, self % matRR(m) % sigmaT(g)
          self % matRR(m) % nuSigmaF(g) = real(mat % getNuFissionXS(g, rand),defFlt)
          !print *, self % matRR(m) % nuSigmaF(g) 
          self % matRR(m) % sigmaF(g) = real(mat % getFissionXS(g, rand),defFlt)
          !print *, self % matRR(m) % sigmaF(g)
          if (self % matRR(m) % nuSigmaF(g) > 0) fiss = .true.
          self % matRR(m) % chi(g) = real(mat % getChi(g, rand),defFlt)
          !print *, self % matRR(m) % chi(g)
          !print *, "Chi stored"
          ! Include scattering multiplicity
          do g1 = 1, self % nG
            self % matRR(m) % sigmaS(self % nG * (g - 1) +  g1)  = &
                    real(mat % getScatterXS(g1, g, rand) * mat % scatter % prod(g1, g) , defFlt)
            !print *, self % matRR(m) % sigmaS(self % nG * (g - 1) +  g1)
            !print *, "Storing scattering"
          end do
        end do
        self % fissile(m) = fiss
        !print *, fiss
        self % names(m) = mm_matName(m)
      end do
      !print *, "XSs stored"

    end if 

    !Initialise data necessary for kinetic/noise calculations
    if (self % doKinetics) then
      print *,'Including kinetic data'
      call fatalError(Here,'Kinetic data not yet supported')

    end if

    aniOrder = 0
    ! Initialise higher-order scattering matrices
    if (aniOrder > 0) then
      print *,'Including anisotropic scattering data'
      call fatalError(Here,'Anisotropy not yet supported')

    end if

  end subroutine init

  !!
  !! Change cross sections for adjoints simulations:
  !! - swap chi and nuSigmaF
  !! - transpose scattering matrices
  !!
  subroutine setAdjointXS(self)
    class(dataRR), intent(inout)            :: self
    integer(shortInt)                       :: g, g1, m
    real(defFlt), dimension(:), allocatable :: buffer
    character(100), parameter :: Here = 'setAdjointXS (dataRR_class.f90)'

    if (self % nT > 1) then
      call fatalError(Here,'Adjoint not supported for coupled problems')
    end if

    allocate(buffer(self % nG))
    do m = 1, self % nMat
      allocate(buffer(self % nG))
      buffer = self % matRR(m) % chi
      self % matRR(m) % chi = self % matRR(m) % nuSigmaF
      self % matRR(m) % nuSigmaF = self % matRR(m) % chi
    end do
    deallocate(buffer)

    ! Transpose scattering matrix
    allocate(buffer(self % nG * self % nG))
    do m = 1, self % nMat
      buffer = self % matRR(m) % sigmaS
        do g = 1, self % nG
          do g1 = 1, self % nG
            self % matRR(m) % sigmaS(self % nG * (g1 - 1) + g)  =  buffer(self % nG * (g - 1) + g1)
        end do
      end do
    end do


    ! TODO: apply the same treatment to anisotropic scattering matrices

  end subroutine setAdjointXS

  !!
  !! Calculate the lower and upper indices for accessing the XS array
  !! (excluding scattering and kinetic data)
  !!
  pure subroutine getIdxs(self, matIdx, idx1, idx2)
    class(dataRR), intent(in)      :: self
    integer(shortInt), intent(in)  :: matIdx
    integer(shortInt), intent(out) :: idx1, idx2

    idx1 = (matIdx - 1) * self % nG + 1
    idx2 = matIdx * self % nG 

  end subroutine getIdxs

  !!
  !! Calculate the lower and upper indices for accessing the scattering XS array
  !!
  pure subroutine getScatterIdxs(self, matIdx, idx1, idx2)
    class(dataRR), intent(in)      :: self
    integer(shortInt), intent(in)  :: matIdx
    integer(shortInt), intent(out) :: idx1, idx2

    idx1 = (matIdx - 1) * self % nG2 + 1
    idx2 = matIdx * self % nG2 

  end subroutine getScatterIdxs

  !!
  !! Return if a material is fissile
  !!
  elemental function isFissile(self, matIdx) result(isFiss)
    class(dataRR), intent(in)     :: self
    integer(shortInt), intent(in) :: matIdx
    logical(defBool)              :: isFiss
    integer(shortInt)             :: mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    isFiss = self % fissile(mIdx)

  end function isFissile
  
  !!
  !! Return the number of groups
  !!
  elemental function getNG(self) result(nG)
    class(dataRR), intent(in) :: self
    integer(shortInt)         :: nG

    nG = self % nG

  end function getNG

  !!
  !! Return the number of temperature files
  !!
  elemental function getNT(self) result(nT)
    class(dataRR), intent(in) :: self
    integer(shortInt)         :: nT

    nT = self % nT

  end function getNT

  !!
  !! Return the number of materials
  !!
  elemental function getNMat(self) result(nM)
    class(dataRR), intent(in) :: self
    integer(shortInt)         :: nM

    nM = self % nMat

  end function getNMat

  !!
  !! Return the number of precursors
  !!
  elemental function getNPrec(self) result(nP)
    class(dataRR), intent(in) :: self
    integer(shortInt)         :: nP

    nP = self % nP

  end function getNPrec
  
  !!
  !! Return the name of a given material
  !!
  elemental function getName(self, matIdx) result(matName)
    class(dataRR), intent(in)     :: self
    integer(shortInt), intent(in) :: matIdx
    character(nameLen)            :: matName

    matName = self % names(matIdx)

  end function getName
  
  
  !!
  !! Get scatter pointer
  !!
  subroutine getScatterPointer(self, matIdx, sigS)
    class(dataRR), target, intent(in)                :: self
    integer(shortInt), intent(in)                    :: matIdx
    real(defFlt), dimension(:), pointer, intent(out) :: sigS
    integer(shortInt)                                :: idx1, idx2, mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    call self % getScatterIdxs(mIdx, idx1, idx2)
    sigS => self % matRR(mIdx) % sigmaS(idx1:idx2)

  end subroutine getScatterPointer
  
  !!
  !! Get scatter vector pointer
  !!
  subroutine getScatterVecPointer(self, matIdx, gOut, sigS)
    class(dataRR), target, intent(in)                :: self
    integer(shortInt), intent(in)                    :: matIdx
    integer(shortInt), intent(in)                    :: gOut
    real(defFlt), dimension(:), pointer, intent(out) :: sigS
    integer(shortInt)                                :: idx1, idx2, mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    idx1 = (matIdx - 1) * self % nG2 + (gOut - 1) * self % nG + 1
    idx2 = (matIdx - 1) * self % nG2 + gOut * self % nG 
    sigS => self % matRR(mIdx) % sigmaS(idx1:idx2)

  end subroutine getScatterVecPointer


  !!
  !! Get chi pointer
  !!
  subroutine getChiPointer(self, matIdx, chi)
    class(dataRR), target, intent(in)                :: self
    integer(shortInt), intent(in)                    :: matIdx
    real(defFlt), dimension(:), pointer, intent(out) :: chi
    integer(shortInt)                                :: idx1, idx2, mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    call self % getIdxs(mIdx, idx1, idx2)
    chi => self % matRR(mIdx) % chi(idx1:idx2)

  end subroutine getChiPointer

  !!
  !! Return pointers to all commonly used XSs for neutron production
  !! This is done for a given material, across all energies
  !!
  subroutine getProdPointers(self, matIdx, nuSigF, sigS, chi)
    class(dataRR), target, intent(in)                :: self
    integer(shortInt), intent(in)                    :: matIdx
    real(defFlt), dimension(:), pointer, intent(out) :: nuSigF, sigS, chi
    integer(shortInt)                                :: idx1, idx2, idx1s, idx2s, mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    call self % getIdxs(mIdx, idx1, idx2)
    call self % getScatterIdxs(mIdx, idx1s, idx2s)
    nuSigF => self % matRR(mIdx) % nuSigmaF(idx1:idx2)
    chi    => self % matRR(mIdx) % chi(idx1:idx2)
    sigS   => self % matRR(mIdx) % sigmaS(idx1s:idx2s)

  end subroutine getProdPointers
  
  !!
  !! Return pointers to only the total XS
  !! This is done for a given material, across all energies
  !!
  subroutine getTotalPointer(self, matIdx, sigT)
    class(dataRR), target, intent(in)                :: self
    integer(shortInt), intent(in)                    :: matIdx
    real(defFlt), dimension(:), pointer, intent(out) :: sigT
    integer(shortInt)                                :: idx1, idx2, mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    call self % getIdxs(mIdx, idx1, idx2)
    sigT => self % matRR(mIdx) % sigmaT(idx1:idx2)

  end subroutine getTotalPointer
  
  !!
  !! Return pointers to only the nuFission XS
  !! This is done for a given material, across all energies
  !!
  subroutine getNuFissPointer(self, matIdx, nuFiss)
    class(dataRR), target, intent(in)                :: self
    integer(shortInt), intent(in)                    :: matIdx
    real(defFlt), dimension(:), pointer, intent(out) :: nuFiss
    integer(shortInt)                                :: idx1, idx2, mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    call self % getIdxs(mIdx, idx1, idx2)
    nuFiss => self % matRR(mIdx) % nuSigmaF(idx1:idx2)

  end subroutine getNuFissPointer

  !! Return nuFission XS in a given material and group 
  !! 
  function getNuFissXS(self, matIdx, temp) result(nuFiss) 
    class(dataRR), intent(in)           :: self 
    integer(shortInt), intent(in)       :: matIdx 
    real(defFlt), optional, intent(in)  :: temp 
    real(defFlt), dimension(self % nG)  :: nuFiss 
    integer(shortInt)                   :: mIdx, g, tIdx 
 
    if (matIdx > self % nMat) then 
      mIdx = self % nMat + 1 
    else 
      mIdx = matIdx 
    end if 
 
    if (matIdx < 1 .or. matIdx > self % nMat) then
      nuFiss = ZERO
      return
    end if

    associate (m => self % matRR(mIdx)) 
      if (present(temp)) then 
 
        if (temp < self % temperatures(mIdx,1)) then 
 
          !$omp simd
          do g = 1, self % nG
            nuFiss(g) = m % nuSigmaF((g - 1) * self % nT + 1)
          end do
 
        else if (temp > self % temperatures(mIdx,3)) then 
 
          !$omp simd
          do g = 1, self % nG
            nuFiss(g) = m % nuSigmaF((g - 1) * self % nT + 3)
          end do
 
        else
 
          tIdx = 1 + merge(1, 0, temp > self % temperatures(mIdx,2))
 
          !$omp simd
          do g = 1, self % nG
            nuFiss(g) = m % nuSigmaF((g - 1) * self % nT + tIdx) + &
                        ((temp - self % temperatures(mIdx,tIdx)) / &
                         (self % temperatures(mIdx,tIdx+1) - &
                          self % temperatures(mIdx,tIdx))) * &
                        (m % nuSigmaF((g - 1) * self % nT + tIdx+1) - &
                         m % nuSigmaF((g - 1) * self % nT + tIdx))
          end do
 
        end if
 
      else
 
        !$omp simd
        do g = 1, self % nG
          nuFiss(g) = m % nuSigmaF(g)
        end do
 
      end if
    end associate
 
  end function getNuFissXS

  !! 
  !! Return fission XS in a given material and group 
  !! 
  function getFissionXS(self, matIdx, temp) result(sigF) 
    class(dataRR), intent(in)           :: self 
    integer(shortInt), intent(in)       :: matIdx 
    real(defFlt), optional, intent(in)  :: temp 
    real(defFlt), dimension(self % nG)  :: sigF 
    integer(shortInt)                   :: mIdx, g, tIdx 
 
    if (matIdx > self % nMat) then 
      mIdx = self % nMat + 1 
    else 
      mIdx = matIdx 
    end if 

    if (matIdx < 1 .or. matIdx > self % nMat) then
      sigF = ZERO
      return
    end if
 
    associate (m => self % matRR(mIdx)) 
      if (present(temp)) then 
 
        if (temp < self % temperatures(mIdx,1)) then 
 
          !$omp simd
          do g = 1, self % nG
            sigF(g) = m % sigmaF((g - 1) * self % nT + 1)
          end do
 
        else if (temp > self % temperatures(mIdx,3)) then 
 
          !$omp simd
          do g = 1, self % nG
            sigF(g) = m % sigmaF((g - 1) * self % nT + 3)
          end do
 
        else
 
          tIdx = 1 + merge(1, 0, temp > self % temperatures(mIdx,2))
 
          !$omp simd
          do g = 1, self % nG
            sigF(g) = m % sigmaF((g - 1) * self % nT + tIdx) + &
                      ((temp - self % temperatures(mIdx,tIdx)) / &
                       (self % temperatures(mIdx,tIdx+1) - &
                        self % temperatures(mIdx,tIdx))) * &
                      (m % sigmaF((g - 1) * self % nT + tIdx+1) - &
                       m % sigmaF((g - 1) * self % nT + tIdx))
          end do
 
        end if
 
      else
 
        !$omp simd
        do g = 1, self % nG
          sigF(g) = m % sigmaF(g)
        end do
 
      end if
    end associate
 
  end function getFissionXS


  !! Return total XS in a given material and group 
  !! 
  function getTotalXS(self, matIdx, temp) result(sigT) 
    class(dataRR), intent(in)          :: self 
    integer(shortInt), intent(in)     :: matIdx 
    real(defFlt), optional, intent(in) :: temp 
    real(defFlt), dimension(self % nG) :: sigT 
    integer(shortInt)                  :: mIdx, g, tIdx 
 
    if (matIdx > self % nMat) then 
      !print *, matIdx
      mIdx = self % nMat + 1 
      !print *, mm_matname(mIdx)
    else 
      mIdx = matIdx 
    end if

    if (matIdx < 1 .or. matIdx > self % nMat) then
      sigT = ZERO
      return
    end if
 
    associate (m => self % matRR(mIdx)) 
      if (present(temp)) then 

        if (temp < self % temperatures(mIdx,1)) then 

          !$omp simd
          do g = 1, self % nG
            sigT(g) = m % sigmaT((g - 1) * self % nT + 1)
          end do

        else if (temp > self % temperatures(mIdx,3)) then 

          !$omp simd
          do g = 1, self % nG
            sigT(g) = m % sigmaT((g - 1) * self % nT + 3)
          end do

        else

          tIdx = 1 + merge(1, 0, temp > self % temperatures(mIdx,2))

          !$omp simd
          do g = 1, self % nG
            sigT(g) = m % sigmaT((g - 1) * self % nT + tIdx) + &
                      ((temp - self % temperatures(mIdx,tIdx)) / &
                       (self % temperatures(mIdx,tIdx+1) - &
                        self % temperatures(mIdx,tIdx))) * &
                      (m % sigmaT((g - 1) * self % nT + tIdx+1) - &
                       m % sigmaT((g - 1) * self % nT + tIdx))
          end do

        end if

      else

        !$omp simd
        do g = 1, self % nG
          sigT(g) = m % sigmaT(g)
        end do

      end if
    end associate

  end function getTotalXS

  !!
  !! Return scatter XS in a given material, ingoing group, and outgoing group
  !!
  function getScatterXS(self, matIdx, gOut, temp) result(sigS)

    class(dataRR), intent(in)     :: self
    integer(shortInt), intent(in) :: matIdx, gOut
    real(defFlt), optional, intent(in) :: temp

    real(defFlt), dimension(self % nG) :: sigS

    integer(shortInt) :: mIdx, gIn, tIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if

    if (matIdx < 1 .or. matIdx > self % nMat) then
      sigS = ZERO
      return
    end if

    associate (m => self % matRR(mIdx))

      if (present(temp)) then

        if (temp < self % temperatures(mIdx,1)) then

          !$omp simd
          do gIn = 1, self % nG
            sigS(gIn) = m % sigmaS(((gIn - 1) * self % nG + (gOut - 1)) * &
                                   self % nT + 1)
          end do

        else if (temp > self % temperatures(mIdx,3)) then

          !$omp simd
          do gIn = 1, self % nG
            sigS(gIn) = m % sigmaS(((gIn - 1) * self % nG + (gOut - 1)) * &
                                   self % nT + 3)
          end do

        else

          tIdx = 1 + merge(1, 0, temp > self % temperatures(mIdx,2))

          !$omp simd
          do gIn = 1, self % nG

            sigS(gIn) = &
              m % sigmaS(((gIn - 1) * self % nG + (gOut - 1)) * &
                         self % nT + tIdx) &
              + ((temp - self % temperatures(mIdx,tIdx)) / &
                 (self % temperatures(mIdx,tIdx+1) - &
                  self % temperatures(mIdx,tIdx))) &
              * (m % sigmaS(((gIn - 1) * self % nG + (gOut - 1)) * &
                             self % nT + tIdx + 1) &
                 - m % sigmaS(((gIn - 1) * self % nG + (gOut - 1)) * &
                               self % nT + tIdx))

          end do

        end if

      else

        !$omp simd
        do gIn = 1, self % nG
          sigS(gIn) = m % sigmaS(self % nG * (gOut - 1) + gIn)
        end do

      end if

    end associate

  end function getScatterXS
  

  !!
  !! Return chi in a given material
  !!
  function getChi(self, matIdx, g) result(chi)
    class(dataRR), intent(in)     :: self
    integer(shortInt), intent(in) :: matIdx, g
    real(defFlt)                  :: chi
    integer(shortInt)             :: mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if

    chi = self % matRR(mIdx) % chi(g)

  end function getChi

  !!
  !! Get evaluated XS temperature for a given
  !! material and temperature index
  !!
  function getTemperature(self, matIdx, t) result(temp)
    class(dataRR), intent(in)     :: self
    integer(shortInt), intent(in) :: matIdx, t
    real(defFlt)                  :: temp
    integer(shortInt)             :: mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if

    temp = self % temperatures(mIdx,t)
  
  end function getTemperature

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(dataRR), intent(inout) :: self

    ! Clean contents
    self % nG         = 0
    self % nG2        = 0
    self % nMat       = 0
    self % nP         = 0
    self % nT         = 0
    self % doKinetics = .false.
    if (allocated(self % matRR)) deallocate(self % matRR)

  end subroutine kill

end module dataRR_class
