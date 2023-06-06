!>
!! @file m_stl.fpp
!! @brief Contains module m_stl

module m_stl

    ! Dependencies =============================================================
    use m_helper
    use m_mpi_proxy
    use m_derived_types

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

    type :: t_stl_bounding_box
        real(kind(0d0)) :: min(3)
        real(kind(0d0)) :: max(3)
    end type t_stl_bounding_box

    type :: t_stl_mesh
        type(t_stl_triangle), dimension(:), allocatable :: triangles
        type(t_stl_bounding_box) :: bounding_box
    end type t_stl_mesh

    public :: f_stl_read, s_stl_write, t_stl_triangle, &
              t_stl_bounding_box, t_stl_mesh, f_stl_is_inside

contains

    function f_stl_transform_matrix(stl) result(out_matrix)

        type(stl_parameters) :: stl

        real(kind(0d0)), dimension(1:4,1:4) :: sc, rz, rx, ry, tr, out_matrix

        sc = transpose(reshape([ &
            stl%scale(1), 0d0,          0d0,          0d0, &
            0d0,          stl%scale(2), 0d0,          0d0, &
            0d0,          0d0,          stl%scale(3), 0d0, &
            0d0,          0d0,          0d0,          1d0  ], shape(sc)))

        rz = transpose(reshape([ &
            cos(stl%rotate(3)), -sin(stl%rotate(3)), 0d0, 0d0, &
            sin(stl%rotate(3)),  cos(stl%rotate(3)), 0d0, 0d0, &
            0d0,                 0d0,                1d0, 0d0, &
            0d0,                 0d0,                0d0, 1d0  ], shape(rz)))

        rx = transpose(reshape([ &
            1d0, 0d0,                 0d0,                0d0, &
            0d0, cos(stl%rotate(1)), -sin(stl%rotate(1)), 0d0, &
            0d0, sin(stl%rotate(1)),  cos(stl%rotate(1)), 0d0, &
            0d0, 0d0,                 0d0,                1d0  ], shape(rx)))

        ry = transpose(reshape([ &
            cos(stl%rotate(2)),  0d0, sin(stl%rotate(2)), 0d0, &
            0d0,                 1d0, 0d0,                0d0, &
            -sin(stl%rotate(2)), 0d0, cos(stl%rotate(2)), 0d0, &
            0d0,                 0d0, 0d0,                1d0  ], shape(ry)))

        tr = transpose(reshape([ &
            1d0, 0d0, 0d0, stl%offset(1), &
            0d0, 1d0, 0d0, stl%offset(2), &
            0d0, 0d0, 1d0, stl%offset(3), &
            0d0, 0d0, 0d0, 1d0            ], shape(tr)))

        out_matrix = matmul(tr, matmul(ry, matmul(rx, matmul(rz, sc))))

    end function f_stl_transform_matrix

    subroutine s_stl_transform_triangle(triangle, matrix)

        type(t_stl_triangle),                intent(inout) :: triangle
        real(kind(0d0)), dimension(1:4,1:4), intent(in)    :: matrix

        integer :: i

        real(kind(0d0)), dimension(1:4) :: tmp

        do i = 1, 3
            tmp = matmul(matrix, [ triangle%v(i,1:3), 1d0 ])
            triangle%v(i, 1:3) = tmp(1:3)
        end do

    end subroutine s_stl_transform_triangle

    !> This procedure reads a binary STL file. The caller is responsible for
    !! freeing the memory allocated for the triangles array.
    !! @param filepath Path to STL file.
    function f_stl_read(params) result(mesh)
    
        type(stl_parameters), intent(IN)  :: params

        type(t_stl_mesh) :: mesh

        integer :: iunit, iostat
        integer :: i, j

        real(kind(0d0)), dimension(1:4,1:4) :: transform

        character(kind=c_char, len=80) :: header
        integer  (kind=c_int32_t)      :: nTriangles

        real   (kind=c_float)          :: normal(3), v(3, 3)
        integer(kind=c_int16_t)        :: attribute

        transform = f_stl_transform_matrix(params)

        open(newunit=iunit,      file=params%filepath, action='READ', &
             form='UNFORMATTED', status='OLD',         iostat=iostat, &
             access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open STL file ", params%filepath
            
            call s_mpi_abort()
        end if

        read(iunit, iostat=iostat) header, nTriangles

        if (iostat /= 0) then
            print *, "Error: could not read header from STL file ", params%filepath

            call s_mpi_abort()
        end if

        ! Initialize bounding box with extreme values (assuming a smaller mesh).
        mesh%bounding_box%min = [  1d30,  1d30,  1d30 ]
        mesh%bounding_box%max = [ -1d30, -1d30, -1d30 ] 

        allocate(mesh%triangles(ntriangles))

        do i = 1, nTriangles
            read(iunit) normal

            mesh%triangles(i)%n = normal

            do j = 1, 3
                read(iunit) v(j,1), v(j,2), v(j,3)
            end do

            mesh%triangles(i)%v = v

            call s_stl_transform_triangle(mesh%triangles(i), transform)

            mesh%bounding_box%min = min(mesh%bounding_box%min, &
                                        minval(mesh%triangles(i)%v, dim=1))
            
            mesh%bounding_box%max = max(mesh%bounding_box%max, &
                                        maxval(mesh%triangles(i)%v, dim=1))

            read(iunit) attribute
        end do

        close(iunit)

    end function f_stl_read

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

    ! From https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution.html
    function f_intersects(ray, triangle) result(intersects)

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

    function f_stl_is_inside(point, mesh, point_spacing) result(inorout)

        integer, parameter :: resolution = 10

        real(kind(0d0)),  intent(in) :: point(1:3)
        type(t_stl_mesh), intent(in) :: mesh
        real(kind(0d0)),  intent(in) :: point_spacing(1:3)

        logical :: inorout

        type(t_stl_ray) :: ray
        integer         :: i, j, hits, totalInCount
    
        real(kind(0d0)), dimension(1:resolution, 1:3) :: ray_origins, ray_dirs
        
        call random_seed()

        do i = 1, resolution
            call random_number(ray_origins(i,:))
            ray_origins(i,:) = point + (ray_origins(i,:) - 0.5) * point_spacing(:)

            call random_number(ray_dirs(i,:))
            ray_dirs(i,:) = ray_dirs(i,:) - 0.5
            ray_dirs(i,:) = ray_dirs(i,:) / sqrt(sum(ray_dirs(i,:) * ray_dirs(i,:)))
        end do

        totalInCount = 0
       
        do i = 1, resolution
            ray%origin    = ray_origins(i,:)
            ray%direction = ray_dirs(i,:)
            hits          = 0
            do j = 1, size(mesh%triangles)
                 if (f_intersects(ray, mesh%triangles(j))) then
                    hits = hits + 1
                end if
            end do

            totalInCount = totalInCount + mod(hits, 2)
        end do

        ! TODO: totalInCount / precision can be used for smoothing for grid
        ! points partially inside and outside the model.
        inorout = (totalInCount / 2d0) .ge. resolution

    end function f_stl_is_inside

end module m_stl
