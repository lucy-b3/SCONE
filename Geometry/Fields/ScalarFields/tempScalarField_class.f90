module tempScalarField_class

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use dictionary_class,  only : dictionary
  use dictParser_func,   only : fileToDict
  use particle_class,    only : particle, particleState
  use field_inter,       only : field
  use scalarField_inter, only : scalarField

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: tempScalarField_TptrCast

  !!
  !! Temperature Scalar Field
  !!
  !! Returns the temperature at a given location
  !!
  !! Sample Dictionary Input:
  !!   temperatures { type tempScalarField; file /home/temperatures.txt ;}
  !!
  !! Sample temperature file:
  !! map{}
  !! values{}
  !!
  !!
  !! Public Members:
  !!   net ->  map that lays over the geometry. Can be any type of Tally Map;
  !!   N   ->  total number of map bins
  !!
  !! Interface:
  !!   scalarField interface
  !!
  type, public, extends(scalarField) :: tempScalarField
    class(tallyMap), allocatable :: net
    integer(shortInt)            :: N
    real(defReal), dimension(:), allocatable   :: temperatures

  contains
    ! Superclass interface
    procedure :: init
    procedure :: kill
    procedure :: at
  end type tempScalarField

contains

  !!
  !! Initialise from dictionary
  !!
  !! See field_inter for details
  !!
  subroutine init(self, dict)
    class(tempScalarField), intent(inout) :: self
    class(dictionary), intent(in)            :: dict
    type(dictionary)                         :: dict2
    character(pathLen)                       :: path
    class(tallyMap), allocatable             :: tempNet
    integer(shortInt), parameter  :: ALL = 0
    character(100), parameter :: Here = 'init (tempScalarField_class.f90)'

    if (allocated(self % net)) deallocate(self % net)
    call dict % get(path,'file')

    ! Load dictionary
    call fileToDict(dict2, path)

    ! Initialise overlay map
    call new_tallyMap(self % net, dict2 % getDictPtr('map'))
    !call new_tallyMap(tempNet, dict2 % getDictPtr('map'))
    !self % net = tempNet
    self % N = self % net % bins(ALL)

    ! Read values from file
    call dict2 % get(self % temperatures, 'values')

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(tempScalarField), intent(inout) :: self

    call self % net % kill()
    if (allocated(self % net)) deallocate(self % net)
    self % N = 0

  end subroutine kill

  !!
  !! Get value of the scalar field given the phase-space location of a particle
  !!
  !! See scalarField_inter for details
  !!
  function at(self, p) result(val)
    class(tempScalarField), intent(in)    :: self
    class(particleState), intent(inout)   :: p
    real(defReal)                         :: val
    type(particleState)                   :: state
    integer(shortInt)                     :: binIdx

    ! Get current particle state
    state = p

    ! Read map bin index
    binIdx = self % net % map(state)

    ! Return if invalid bin index
    if (binIdx == 0) then
      val = ZERO
      return
    end if

    val = self % temperatures(binIdx)

  end function at

  !!
  !! Cast field pointer to tempScalarField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null if source is not of tempScalarField
  !!   Pointer to source if source is tempScalarField type
  !!
  pure function tempScalarField_TptrCast(source) result(ptr)
    class(field), pointer, intent(in) :: source
    type(tempScalarField), pointer :: ptr

    select type (source)
      type is (tempScalarField)
        ptr => source

      class default
        ptr => null()
    end select

  end function tempScalarField_TptrCast


end module tempScalarField_class
