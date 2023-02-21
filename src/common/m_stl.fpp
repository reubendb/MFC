!>
!! @file m_stl.fpp
!! @brief Contains module m_stl

module m_stl

    ! Dependencies =============================================================
    use m_mpi_proxy

    use iso_c_binding, only: c_int16_t, c_int32_t, c_float, c_char

    ! ==========================================================================

    implicit none

    private

    type :: t_stl_triangle
        real(kind(0d0)) :: v(3, 3) ! Vertices of triangle
        real(kind(0d0)) :: n(3)    ! Normal vector
    end type t_stl_triangle

    type :: t_stl_ray
        real(kind(0d0)) :: o(3) ! Ray origin
        real(kind(0d0)) :: d(3) ! Ray direction
    end type t_stl_ray

    public :: s_stl_read, s_stl_write, t_stl_triangle, f_stl_is_inside

contains

    function crossprod(a, b) RESULT(r)
        real(kind(0d0)), dimension(1:3), INTENT(in) :: a, b
        real(kind(0d0)), dimension(1:3) :: r

        r(1) = a(2) * b(3) - a(3) * b(2)
        r(2) = a(3) * b(1) - a(1) * b(3)
        r(3) = a(1) * b(2) - a(2) * b(1)
    end function crossprod

    !> This procedure reads a binary STL file. The caller is responsible for
    !! freeing the memory allocated for the triangles array.
    !! @param filepath Path to STL file.
    subroutine s_stl_read(filepath, trs)
    
        character(LEN=*),                                intent(IN)  :: filepath
        type(t_stl_triangle), dimension(:), allocatable, intent(OUT) :: trs

        integer :: iunit, iostat
        integer :: i, j

        character(kind=c_char, len=80) :: header
        integer  (kind=c_int32_t)      :: nTrs

        real   (kind=c_float)          :: normal(3), v(3, 3)
        integer(kind=c_int16_t)        :: attribute

        open(newunit=iunit,      file=filepath, action='READ', &
             form='UNFORMATTED', status='OLD',  iostat=iostat, &
             access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open STL file ", filepath
            
            call s_mpi_abort()
        end if

        read(iunit, iostat=iostat) header, nTrs

        if (iostat /= 0) then
            print *, "Error: could not read header from STL file ", filepath

            call s_mpi_abort()
        end if

        allocate(trs(nTrs))

        do i = 1, nTrs
            ! Skip the normal vector as it is not reliable.
            read(iunit) normal(1), normal(2), normal(3)

            do j = 1, 3
                read(iunit) v(j,1), v(j,2), v(j,3)
            end do

            ! We compute the normal vector to ensure that it is correct
            trs(i)%v = v
            trs(i)%n = crossprod(&
                trs(i)%v(2,1:3) - trs(i)%v(1,1:3),&
                trs(i)%v(3,1:3) - trs(i)%v(1,1:3))

            read(iunit) attribute
        end do

        close(iunit)

    end subroutine s_stl_read

    subroutine s_stl_write(filepath, triangles)
    
        character(LEN=*),     intent(IN) :: filepath
        type(t_stl_triangle), intent(IN) :: triangles(:)

        integer :: iunit, iostat
        integer :: i, j

        character(kind=c_char, len=80) :: header
        integer  (kind=c_int32_t)      :: nTriangles

        real(kind=c_float) :: normal(3), v(3, 3)

        integer(kind=c_int16_t), parameter :: attribute = 0

        open(newunit=iunit,      file=filepath, action='WRITE', &
             form='UNFORMATTED', iostat=iostat, access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open STL file ", filepath
            
            call s_mpi_abort()
        end if

        header     = "STL file written by MFC."
        ntriangles = size(triangles)

        write(iunit, iostat=iostat) header, ntriangles

        if (iostat /= 0) then
            print *, "Error: could not write header to STL file ", filepath

            call s_mpi_abort()
        end if

        do i = 1, ntriangles
            normal = triangles(i)%n

            write(iunit) normal

            do j = 1, 3
                v = triangles(i)%v
                write(iunit) v(j,1), v(j,2), v(j,3)
            end do

            write(iunit) attribute
        end do

        close(iunit)

    end subroutine s_stl_write

    function f_intersects(ray, tr) result(intersects)

        ! Adapted from https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution.html

        type(t_stl_ray),      intent(in) :: ray
        type(t_stl_triangle), intent(in) :: tr

        logical :: intersects

        real(kind(0d0)) :: P(3), C(3)
        real(kind(0d0)) :: t, ndotrdir

        integer :: i ! Loop index

        intersects = .false.
        ndotrdir   = dot_product(tr%n(:), ray%d(:))

        if (abs(ndotrdir) .lt. verysmall) then; return; endif

        t = ( dot_product(tr%n(:), tr%v(1,:)) &
             -dot_product(tr%n(:), ray%o(:)) ) / ndotrdir
    
        if (t .lt. 0) then; return; endif

        P = ray%o + t * ray%d

        do i = 1, 3
            C = crossprod(&
                    tr%v(1 + mod(i, 3),:) - tr%v(i,:),&
                    P - tr%v(i,:))
            
            if (dot_product(tr%n(:), C(:)) .lt. 0d0) then; return; endif
        end do

        intersects = .true.
        
    end function f_intersects

    function f_stl_is_inside(point, triangles) result(inorout)

        real(kind(0d0)),      intent(in) :: point(3)
        type(t_stl_triangle), intent(IN) :: triangles(:)

        logical :: inorout

        type(t_stl_ray) :: ray
        integer         :: i, hits

        hits  = 0
        ray%o = point
        ray%d = [0d0, 0d0, 1d0] ! This is arbitrary, but must be non-zero
        
        do i = 1, size(triangles)
             if (f_intersects(ray, triangles(i))) then
                hits = hits + 1
            end if
        end do

        inorout = (mod(hits, 2) .eq. 1)

    end function f_stl_is_inside

end module m_stl
