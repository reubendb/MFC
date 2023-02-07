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
        real(kind(0d0)) :: origin(3)
        real(kind(0d0)) :: direction(3)
    end type t_stl_ray

    public :: s_stl_read, s_stl_write, t_stl_triangle, f_stl_is_inside

contains

    !> This procedure reads a binary STL file. The caller is responsible for
    !! freeing the memory allocated for the triangles array.
    !! @param filepath Path to STL file.
    subroutine s_stl_read(filepath, triangles)
    
        character(LEN=*),                                intent(IN)  :: filepath
        type(t_stl_triangle), dimension(:), allocatable, intent(OUT) :: triangles

        integer :: iunit, iostat
        integer :: i, j

        character(kind=c_char, len=80) :: header
        integer  (kind=c_int32_t)      :: nTriangles

        real   (kind=c_float)          :: normal(3), v(3, 3)
        integer(kind=c_int16_t)        :: attribute

        open(newunit=iunit,      file=filepath, action='READ', &
             form='UNFORMATTED', status='OLD',  iostat=iostat, &
             access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open STL file ", filepath
            
            call s_mpi_abort()
        end if

        read(iunit, iostat=iostat) header, ntriangles

        if (iostat /= 0) then
            print *, "Error: could not read header from STL file ", filepath

            call s_mpi_abort()
        end if

        allocate(triangles(ntriangles))

        do i = 1, ntriangles
            read(iunit) normal

            triangles(i)%n = normal

            do j = 1, 3
                read(iunit) v(j,1), v(j,2), v(j,3)
            end do

            triangles(i)%v = v

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

        real   (kind=c_float)          :: normal(3), v(3, 3)
        integer(kind=c_int16_t)        :: attribute

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

            attribute = 0
            write(iunit) attribute
        end do

        close(iunit)

    end subroutine s_stl_write

    function cross(a, b) result(c)
        !$acc routine seq

        real(kind(0d0)), dimension(3), intent(in) :: a, b
        real(kind(0d0)), dimension(3)             :: c

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function cross

    ! From https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution.html
    function f_intersects(ray, triangle) result(intersects)
        !$acc routine seq

        type(t_stl_ray),      intent(in) :: ray
        type(t_stl_triangle), intent(in) :: triangle

        logical :: intersects

        real(kind(0d0)) :: v0v1(3), v0v2(3), N(3), P(3), C(3), edge(3), vp(3)
        real(kind(0d0)) :: area2, d, t, NdotRayDirection

        intersects = .false.

        N     = triangle%n
        area2 = sqrt(sum(N(:) * N(:)))
        
        NdotRayDirection = sum(N(:) * ray%direction(:))

        if (abs(NdotRayDirection) .lt. 0.0000001) then
            return
        end if

        d = - sum(N(:) * triangle%v(1,:))
        t = - (sum(N(:) * ray%origin(:)) + d) / NdotRayDirection
    
        if (t .lt. 0) then
            return
        end if

        P = ray%origin + t * ray%direction

        edge = triangle%v(2,:) - triangle%v(1,:)
        vp   = P - triangle%v(1,:)
        C = cross(edge, vp)
        if (sum(N(:) * C(:)) .lt. 0) then
            return
        end if

        edge = triangle%v(3,:) - triangle%v(2,:)
        vp   = P - triangle%v(2,:)
        C = cross(edge, vp)
        if (sum(N(:) * C(:)) .lt. 0) then
            return
        end if

        edge = triangle%v(1,:) - triangle%v(3,:)
        vp   = P - triangle%v(3,:)
        C = cross(edge, vp)
        if (sum(N(:) * C(:)) .lt. 0) then
            return
        end if

        intersects = .true.
        
    end function f_intersects

    function f_stl_is_inside(point, triangles) result(inorout)
        !$acc routine seq

        real(kind(0d0)),      intent(in) :: point(3)
        type(t_stl_triangle), intent(IN) :: triangles(:)

        logical :: inorout

        type(t_stl_ray) :: ray
        integer         :: i, hits

        hits          = 0
        ray%origin    = point
        ray%direction = [0d0, 0d0, 1d0] ! This is arbitrary, but must be non-zero
        
        do i = 1, size(triangles)
             if (f_intersects(ray, triangles(i))) then
                hits = hits + 1
            end if
        end do

        inorout = (mod(hits, 2) .eq. 1)

    end function f_stl_is_inside

end module m_stl
