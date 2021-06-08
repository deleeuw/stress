  subroutine FOR_SUB(i,x,y,ny)
    implicit none
    integer :: i, ny
    double precision :: x
    double precision :: y(ny)
    print *, "This is in Fortran routine..."
    print *, "i = ", i, ", x = ", x
    print *, "y = ", y(1:ny)
  end subroutine FOR_SUB
