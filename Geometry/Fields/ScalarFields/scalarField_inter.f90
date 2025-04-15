module scalarField_inter

  use numPrecision
  use field_inter,    only : field
  use particle_class, only : particleState, particle
  use dictionary_class,  only : dictionary
  
  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: scalarField_CptrCast

  !!
  !! Simple Real Scalar Field
  !!
  !! Access to field is via coordList to allow more fancy fields to be defined
  !! (e.g. assign value to each uniqueID etc.)
  !!
  !! Interface:
  !!   field interface
  !!   at -> Return scalar value given position coordinates
  !!
  type, public, abstract, extends(field) :: scalarField
  contains
    procedure(at), deferred :: at
    procedure(init), deferred :: init
    procedure(kill), deferred :: kill
  end type scalarField

  abstract interface

    !!
    !! Get value of the scalar field at the co-ordinate point
    !!
    !! Args:
    !!   coords [in] -> Coordinates of the position in the geometry
    !!
    !! Result:
    !!   Value of the scalar field. Real number.
    !!
    function at(self, p) result(val)
      import :: scalarField, particleState, defReal
      class(scalarField), intent(in) :: self
      class(particleState), intent(inout)   :: p
      real(defReal)                  :: val
    end function at

    
    !!
    !! Initialise from dictionary
    !!
    !! See field_inter for details
    !!
    subroutine init(self, dict)
      import :: scalarField, dictionary
      class(scalarField), intent(inout) :: self
      class(dictionary), intent(in)     :: dict
    end subroutine init

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      import :: scalarField
      class(scalarField), intent(inout) :: self
    end subroutine kill

  end interface

contains

  !!
  !! Cast field pointer to scalarField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of scalarField
  !!   Pointer to source if source is scalarField class
  !!
  pure function scalarField_CptrCast(source) result(ptr)
    class(field), pointer, intent(in) :: source
    class(scalarField), pointer       :: ptr

    select type (source)
    class is (scalarField)
        ptr => source

      class default
        ptr => null()
    end select

  end function scalarField_CptrCast

end module scalarField_inter
