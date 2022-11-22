#:def ALLOCATE(*args)
#ifdef MFC_DEBUG
    print*, '@:ALLOCATE(${', '.join(args)}$)'
#endif
    allocate(${', '.join(args)}$)
#ifdef MFC_SIMULATION
    !$acc enter data create(${', '.join(args)}$)
#endif
#:enddef ALLOCATE

#:def DEALLOCATE(*args)
#ifdef MFC_DEBUG
    print*, '@DEALLOCATE(${', '.join(args)}$)'
#endif
    deallocate(${', '.join(args)}$)
#ifdef MFC_SIMULATION
    !$acc exit data delete(${', '.join(args)}$)
#endif
#:enddef DEALLOCATE
