      module timep
        real, allocatable :: tprev(:)
      end module
        
      ! HDHM model S-function
      subroutine HDHM_sfunc(t, n, ind, bin, hout, dhout, ihy, iup)

        use hydatas
        use timep
        implicit none
        integer ihy, n, ind(n), iup
        double precision hvec(n), dhvec(n), dble
 
        ! The inputs and ouputs have to be real
        real t, hout(n), dhout(n), bin(n)

        ! Update only when time step has changed
        if(.not.allocated(tprev)) then
          allocate(tprev(info(1).nmodel))
          tprev(:) = 0.
        endif
        if(t.gt.tprev(ind(1))) then                 ! Assumes currently that n = 1
          iup = 1
          tprev(ind(1)) = t
        else
          iup = 0
        endif
        
        
        ! Evaluate hysteresis model
        if(ihy.ge.1) then
          call evaluate_HDHM(n,ind,dble(bin),hvec,dhvec,iup,0,0)
        
        ! Evaluate only single-valued model
        else
          call evaluate_sv(n,dble(bin),hvec,dhvec,0)
        endif
        hout(:)=hvec(:)
        dhout(:)=dhvec(:)
        
      return
      end


      ! Initialization function (from file)
      subroutine HDHM_sfunc_init(nmodel, nrevmax)
        use hydatas
        use timep
        implicit none
        integer nmodel,nrevmax

        if(allocated(tprev)) deallocate(tprev)

        ! Initialize model structures from file
        call deallocate_HDHM
        call initFromFile(nmodel, nrevmax,
     &   './hyst/hydata_nippon_35h300.txt')
!        
      return
      end
        
      ! Initialization function (from file)
      subroutine HDHM_sfunc_terminate
        use hydatas
        use timep
        implicit none
        integer nmodel,nrevmax

        if(allocated(tprev)) deallocate(tprev)
        call deallocate_HDHM

      return
      end
