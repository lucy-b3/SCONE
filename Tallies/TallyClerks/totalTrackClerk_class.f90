module totalTrackClerk_class

  use numPrecision
  use tallyCodes
  use genericProcedures,          only : fatalError
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use outputFile_class,           only : outputFile
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk, kill_super => kill

  ! Nuclear Data interface
  use nuclearDatabase_inter,      only : nuclearDatabase

  ! Tally Filters
  use tallyFilter_inter,          only : tallyFilter
  use tallyFilterFactory_func,    only : new_tallyFilter

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  ! Tally Responses
  use tallyResponseSlot_class,    only : tallyResponseSlot

  implicit none
  private

  !!
  !! Stores total track length travelled by all neutrons
  !!
  !! Interface
  !!   tallyClerk Interface
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myTotalTrackClerk {
  !!   type totalTrackClerk;
  !! }
  !!
  type, public, extends(tallyClerk) :: totalTrackClerk
    private

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: kill
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportTrans

    ! Output procedures
    procedure  :: display
    procedure  :: print

  end type totalTrackClerk

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(totalTrackClerk), intent(inout)       :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen), intent(in)              :: name

    ! Assign name
    call self % setName(name)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(totalTrackClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

  end subroutine kill

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(totalTrackClerk),intent(in)               :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [trans_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(totalTrackClerk), intent(in)     :: self
    integer(shortInt)                 :: S

    S = 1

  end function getSize

  !!
  !! Process particle transition reports
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportTrans(self, p, xsData, mem)
    class(totalTrackClerk), intent(inout) :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    type(particleState)                   :: pBefore, pNow
    integer(longInt)                      :: adrr
    real(defReal)                         :: scoreVal, L
    character(100), parameter :: Here =' reportPath (totalTrackClerk_class.f90)'

    ! Get pre-transition particle state
    pBefore = p % preTransition

    ! Calculate bin address
    adrr = self % getMemAddress()

    ! tranfer information about Prestate material to a temporary particle
    pNow = p

    ! Calculate flux sample L = path travelled
    L = norm2(pNow % r - pBefore % r)

    ! Update tally
    scoreVal = L
    call mem % score(scoreVal, adrr)

  end subroutine reportTrans

  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(totalTrackClerk), intent(in) :: self
    type(scoreMemory), intent(in)      :: mem

    print *, 'totalTrackClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(totalTrackClerk), intent(in)         :: self
    class(outputFile), intent(inout)           :: outFile
    type(scoreMemory), intent(in)              :: mem
    real(defReal)                              :: val, std
    character(nameLen)                         :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    name ='tot_path'
    call mem % getResult(val, std, self % getMemAddress())
    call outFile % printResult(val, std, name)

    call outFile % endBlock()

  end subroutine print

end module totalTrackClerk_class
