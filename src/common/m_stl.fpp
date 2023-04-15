!>
!! @file   m_stl.fpp
!! @author Henry Le Berre <hberre3@gatech.edu>
!! @brief  Contains module m_stl

module m_stl

    ! Dependencies =============================================================
    use m_helper
    use m_mpi_proxy
    use m_derived_types

    use iso_c_binding, only: c_int16_t, c_int32_t, c_float, c_char

    ! ==========================================================================

    implicit none

    private

    type :: t_triangle
        real(kind(0d0)) :: v(3, 3) ! Vertices of triangle
        real(kind(0d0)) :: n(3)    ! Normal vector
    end type t_triangle

    type :: t_ray
        real(kind(0d0)) :: o(3) ! Origin
        real(kind(0d0)) :: d(3) ! Direction
    end type t_ray

    type :: t_bbox
        real(kind(0d0)) :: min(3)
        real(kind(0d0)) :: max(3)
    end type t_bbox

    type :: t_stl
        type(t_triangle), dimension(:), allocatable :: trs
        type(t_bbox) :: bbox
    end type t_stl

    type :: t_octree
        type(t_triangle), pointer, dimension(:)     :: trs
        type(t_octree),   pointer, dimension(:,:,:) :: children
        type(t_bbox) :: bbox
    end type

    public :: t_bbox, t_stl, t_triangle, t_octree, &
              f_stl_read, s_stl_write, s_stl_free, &
              f_octree_create, s_octree_free, f_stl_is_inside

contains

    !> This procedure creates a transformation matrix from the STL parameters.
    !! @param  stl Parameters for the transformation.
    !! @return Transformation matrix.
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

    !> This procedure transforms a triangle by a matrix, one vertex at a time.
    !! @param triangle Triangle to transform.
    !! @param matrix   Transformation matrix.
    subroutine s_stl_transform_triangle(triangle, matrix)

        type(t_triangle),                intent(inout) :: triangle
        real(kind(0d0)), dimension(1:4,1:4), intent(in)    :: matrix

        integer :: i

        real(kind(0d0)), dimension(1:4) :: tmp

        do i = 1, 3
            tmp = matmul(matrix, [ triangle%v(i,1:3), 1d0 ])
            triangle%v(i, 1:3) = tmp(1:3)
        end do

    end subroutine s_stl_transform_triangle

    !> This procedure reads a binary STL file.
    !! @param params STL Parameters (filepath & transformation).
    !! @return       STL mesh.
    function f_stl_read(params) result(mesh)
    
        type(stl_parameters), intent(IN)  :: params

        type(t_stl) :: mesh

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
        mesh%bbox%min = [  1d30,  1d30,  1d30 ]
        mesh%bbox%max = [ -1d30, -1d30, -1d30 ] 

        allocate(mesh%trs(ntriangles))

        do i = 1, nTriangles
            read(iunit) normal

            mesh%trs(i)%n = normal

            do j = 1, 3
                read(iunit) v(j,1), v(j,2), v(j,3)
            end do

            mesh%trs(i)%v = v

            call s_stl_transform_triangle(mesh%trs(i), transform)

            mesh%bbox%min = min(mesh%bbox%min, &
                                        minval(mesh%trs(i)%v, dim=1))
            
            mesh%bbox%max = max(mesh%bbox%max, &
                                        maxval(mesh%trs(i)%v, dim=1))

            read(iunit) attribute
        end do

        close(iunit)

    end function f_stl_read

    !> This procedure writes a binary STL file.
    !! @param filepath  Path to the file to write.
    !! @param triangles Triangles to write.
    subroutine s_stl_write(filepath, triangles)
    
        character(LEN=*),     intent(IN) :: filepath
        type(t_triangle), intent(IN) :: triangles(:)

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

    !> This procedure frees the memory allocated for an STL mesh.
    subroutine s_stl_free(mesh)

        type(t_stl), intent(INOUT) :: mesh

        deallocate(mesh%trs)

    end subroutine s_stl_free

    ! From https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution.html
    !> This procedure checks if a ray intersects a triangle.
    !! @param ray      Ray.
    !! @param triangle Triangle.
    !! @return         True if the ray intersects the triangle, false otherwise.
    function f_intersects_triangle(ray, triangle) result(intersects)

        type(t_ray),      intent(in) :: ray
        type(t_triangle), intent(in) :: triangle

        logical :: intersects

        real(kind(0d0)) :: v0v1(3), v0v2(3), N(3), P(3), C(3), edge(3), vp(3)
        real(kind(0d0)) :: area2, d, t, NdotRayDirection

        intersects = .false.

        N     = triangle%n
        area2 = sqrt(sum(N(:) * N(:)))
        
        NdotRayDirection = sum(N(:) * ray%d(:))

        if (abs(NdotRayDirection) .lt. 0.0000001) then
            return
        end if

        d = - sum(N(:) * triangle%v(1,:))
        t = - (sum(N(:) * ray%o(:)) + d) / NdotRayDirection
    
        if (t .lt. 0) then
            return
        end if

        P = ray%o + t * ray%d

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
        
    end function f_intersects_triangle

    !> This procedure swaps two real numbers.
    !! @param lhs Left-hand side.
    !! @param rhs Right-hand side.
    subroutine s_swap(lhs, rhs)
        real(kind(0d0)), intent(inout) :: lhs, rhs
        real(kind(0d0))                :: ltemp
        ltemp = lhs; lhs = rhs; rhs = ltemp
    end subroutine s_swap

    !> This procedure checks if a ray intersects a bounding box.
    !! @param ray  Ray.
    !! @param bbox Bounding box.
    !! @return     True if the ray intersects the bounding box, false otherwise.
    function f_ray_intersects_bbox(ray, bbox) result(intersects)
        
        type(t_ray),  intent(in) :: ray
        type(t_bbox), intent(in) :: bbox

        logical :: intersects

        integer :: i
        
        real(kind(0d0)) :: tmin, tmax, tymin, tymax, tzmin, tzmax

        tmin = (bbox%min(1) - ray%o(1)) / ray%d(1)
        tmax = (bbox%max(1) - ray%o(1)) / ray%d(1)

        if (tmin .gt. tmax) then
            call s_swap(tmin, tmax)
        end if

        tymin = (bbox%min(2) - ray%o(2)) / ray%d(2)
        tymax = (bbox%max(2) - ray%o(2)) / ray%d(2)

        if (tymin .gt. tymax) then
            call s_swap(tymin, tymax)
        end if

        if (tmin .gt. tymax .or. tymin .gt. tmax) then
            intersects = .false.
            return
        end if

        if (tymin .gt. tmin) then
            tmin = tymin
        end if

        if (tymax .lt. tmax) then
            tmax = tymax
        end if

        tzmin = (bbox%min(3) - ray%o(3)) / ray%d(3)
        tzmax = (bbox%max(3) - ray%o(3)) / ray%d(3)

        if (tzmin .gt. tzmax) then
            call s_swap(tzmin, tzmax)
        end if

        if (tmin .gt. tzmax .or. tzmin .gt. tmax) then
            intersects = .false.
            return
        end if

        intersects = .true.

    end function f_ray_intersects_bbox

    function f_triangle_intersects_bbox(tr, bbox) result(intersects)

        type(t_triangle), intent(in) :: tr
        type(t_bbox),     intent(in) :: bbox

        logical :: intersects

        intersects = .true.

    end function f_triangle_intersects_bbox

    !> This procedure builds, recursively, an octree from a set of triangles.
    !! @param node  Octree node to build the children of.
    !! @param bbox  Bounding box of the node.
    !! @param trs   Set of triangles in bbox.
    !! @param ntrs  Number of triangles in trs.
    !! @param depth Current depth of the octree.
    recursive subroutine s_octree_impl(node, bbox, trs, ntrs, depth)
        
        !FIXME: These values sound *reasonable*, we should maybe find better ones.
        integer, parameter :: max_depth     = 1
        integer, parameter :: max_triangles = 10

        type(t_octree), intent(inout) :: node
        type(t_bbox), intent(in) :: bbox
        type(t_triangle), allocatable, dimension(:), intent(in) :: trs
        integer, intent(in) :: ntrs
        integer, intent(in) :: depth

        type(t_bbox) :: child_bbox
        type(t_triangle), allocatable, dimension(:) :: child_trs
        integer :: child_ntrs

        logical :: intersects

        integer :: i, j, k, l, q

        node%bbox     =  bbox
        node%trs      => null()
        node%children => null()

        if ((depth .ge. max_depth) .or. (ntrs .le. max_triangles)) then
            if (ntrs .ne. 0) then
                allocate(node%trs(1:ntrs))
                node%trs(1:ntrs) = trs(1:ntrs)
            end if

            return
        end if

        allocate(node%children(1:2, 1:2, 1:2))
        allocate(child_trs(1:ntrs))
        
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    child_ntrs = 0
                    
                    child_bbox%min = bbox%min + &
                        0.5 * (bbox%max - bbox%min) * (/ i - 1, j - 1, k - 1 /)
                    child_bbox%max = bbox%min + &
                        0.5 * (bbox%max - bbox%min) * (/ i, j, k /)
                    
                    do l = 1, ntrs
                        if (f_triangle_intersects_bbox(trs(l), child_bbox)) then
                            child_ntrs = child_ntrs + 1
                            child_trs(child_ntrs) = trs(l)
                        end if
                    end do

                    call s_octree_impl(node%children(i,j,k), child_bbox, &
                        child_trs, child_ntrs, depth + 1)
                end do
            end do
        end do

        deallocate(child_trs)

    end subroutine s_octree_impl

    !> This procedure creates an octree from a mesh.
    !! @param mesh Mesh to create the octree from.
    !! @return Octree created from the mesh.
    function f_octree_create(mesh) result(root)

        type(t_stl), intent(in) :: mesh
        type(t_octree)          :: root

        root%bbox = mesh%bbox

        call s_octree_impl(root, root%bbox, mesh%trs, size(mesh%trs, 1), 1)

    end function f_octree_create

    !> This procedure frees the memory associated with an octree.
    !! @param node Octree to free.
    recursive subroutine s_octree_free(node)

        type(t_octree), intent(inout) :: node

        integer :: i, j, k

        if (associated(node%trs)) then
            deallocate(node%trs)
        elseif (associated(node%children)) then
            do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        call s_octree_free(node%children(i,j,k))
                    end do
                end do
            end do

            deallocate(node%children)
        end if

    end subroutine s_octree_free

    !> This procedure, recursively, counts the number of triangles intersected
    !! by a ray in an octree.
    !! @param ray    Ray.
    !! @param octree Octree to search in.
    !! @return Number of triangles intersected by the ray.
    recursive function f_octree_count_hits(ray, octree) result(hits)

        type(t_ray),    intent(in) :: ray
        type(t_octree), intent(in) :: octree

        integer :: hits

        integer :: i, j, k, l, q

        hits = 0

        if (.not. f_ray_intersects_bbox(ray, octree%bbox)) then
            return
        end if

        if (associated(octree%trs)) then
            do i = 1, size(octree%trs, 1)
                if (f_intersects_triangle(ray, octree%trs(i))) then
                    hits = hits + 1
                end if
            end do
        elseif (associated(octree%children)) then
            do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        hits = hits + f_octree_count_hits(ray, octree%children(i,j,k))
                    end do
                end do
            end do
        end if

    end function f_octree_count_hits

    !> This procedure, recursively, finds whether a point is inside an octree.
    !! @param point         Point to search for.
    !! @param octree        Octree to search in.
    !! @param point_spacing Space around the point to search in (grid spacing).
    !! @return True if the point is inside the octree, false otherwise.
    function f_stl_is_inside(point, octree, point_spacing) result(inorout)

        integer, parameter :: resolution = 10

        real(kind(0d0)), intent(in) :: point(1:3)
        type(t_octree),  intent(in) :: octree
        real(kind(0d0)), intent(in) :: point_spacing(1:3)

        logical :: inorout

        type(t_ray) :: ray
        integer     :: i, j, totalInCount

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
            ray%o = ray_origins(i,:)
            ray%d = ray_dirs(i,:)
            totalInCount = totalInCount + &
                mod(f_octree_count_hits(ray, octree), 2)
        end do

        ! TODO: totalInCount / precision can be used for smoothing for grid
        ! points partially inside and outside the model.
        inorout = (totalInCount / 2d0) .ge. resolution

    end function f_stl_is_inside

end module m_stl
