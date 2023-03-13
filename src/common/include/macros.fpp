#:def LOG(expr)
#ifdef MFC_DEBUG
print *, '${_FILE_.split('/')[-1]}$:${_LINE_}$: ', ${expr}$
#endif
#:enddef


#:def ALLOCATE(*args)
@:LOG({'@:ALLOCATE(${', '.join(args)}$)'})
allocate(${', '.join(args)}$)
#:enddef ALLOCATE


#:def DEALLOCATE(*args)
@:LOG({'@:DEALLOCATE(${', '.join(args)}$)'})
deallocate(${', '.join(args)}$)
#:enddef DEALLOCATE


#:def ACC_SETUP_VFs(*args)
block
    integer :: macros_setup_vfs_i

    @:LOG({'@:ACC_SETUP_VFs(${', '.join(args)}$)'})

    #:for arg in args
        !$acc enter data copyin(${arg}$)
        !$acc enter data copyin(${arg}$%vf)
        if (allocated(${arg}$%vf)) then 
            do macros_setup_vfs_i = lbound(${arg}$%vf, 1), ubound(${arg}$%vf, 1)
                if (associated(${arg}$%vf(macros_setup_vfs_i)%sf)) then
                    !$acc enter data copyin(${arg}$%vf(macros_setup_vfs_i))
                    !$acc enter data create(${arg}$%vf(macros_setup_vfs_i)%sf)
                end if
            end do
        end if
    #:endfor
end block
#:enddef

