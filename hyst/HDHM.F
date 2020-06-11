      ! Datas for the model, everything has to be allocatable in order to keep the variables static
      module hydatas
      
        ! Number of models and maximum number of reversal points
        type :: fo
          integer :: nmodel, nrevmax
        end type
        type(fo), allocatable :: info(:)
      
        ! Model datas
        type :: hy
          integer nm,nav
          double precision, allocatable :: ba_data(:),ha_data(:),
     &     dha_data(:),bav_data(:),hav_data(:),dhav_data(:)
          double precision bs
        end type
        type(hy), allocatable :: hydata(:)
        
        ! History structure
        type :: hi
          double precision, allocatable :: br(:),bt(:),dbout(:)
          integer, allocatable :: drev(:),irevs(:,:)
          integer dp,nrev
          double precision :: bp,hp
          double precision, allocatable ::
     &      has(:,:),hds(:,:),dhas(:,:),dhds(:,:),bas(:,:),bds(:,:)
        end type
        type(hi), allocatable :: hist(:)
      end module

! ##################################################################################################################################################################################################################
      
      ! Routine for evaluating the models
      subroutine evaluate_HDHM(n,ind,bk,h,dhdb,iup,iinv,ider)
      use hydatas

        implicit none
        integer n, iup, iinv, ind(n), it, i, iup2, ider
        double precision bk(n), h(n), dhdb(n), bin, res, href, jac, db

        ! Loop models
        do i=1,n
      
          ! Don't invert
          if(iinv.eq.0) then
          
            call HDHM(ind(i),bk(i),h(i),dhdb(i),iup,ider)
          
          ! Invert
          else
          
            ! Get previous H value
            if(abs(hist(ind(i))%hp-9e6).lt.1e-6) then
              call HDHM(ind(i),hist(ind(i))%bp,hist(ind(i))%hp,
     &                   dhdb(i),0)
            endif

            ! Initial B value based on if H increases or decreases
            href = bk(i)
            if(href.gt.hist(ind(i))%hp) then
              bin = hist(ind(i))%bp + 1d-6
            else
              bin = hist(ind(i))%bp - 1d-6
            endif
           
            ! Iterate for B
            do it = 1,100
              iup2 = 0
              if(it.eq.1) iup2 = 1
              call HDHM(ind(i),bin,h(i),dhdb(i),iup2)
              res = h(i)-href
              jac = dhdb(i)
              db = -res/jac
              if((abs(db/bin).lt.1d-6).or.(abs(res)<1d-6)) exit
              bin = bin + db
            enddo
            
            ! Update model
            call HDHM(ind(i),bin,h(i),dhdb(i),iup)
            dhdb(i) = 1./dhdb(i)
            h(i) = bin
            
          endif
        enddo
      return
      end

      ! Single-valued curve only
      subroutine evaluate_sv(n,bk,h,dhdb,iinv)
      use hydatas

        implicit none
        integer n,iinv,i
        double precision bk(n), h(n), dhdb(n), hav

        ! Loop models
        do i=1,n
      
          ! Don't invert
          if(iinv.eq.0) then
          
            h(n)    = dsign(1d0,bk(n))*hav(abs(bk(n)),0);
            dhdb(n) = hav(abs(bk(n)),1);
          
          ! Invert
          else
             call mexprintf('Single-valued inverted
     &               model not implemented yet.\n')
          endif
        enddo
      return
      end
      
      
! ##################################################################################################################################################################################################################
      
      ! Hysteresis model
      ! Based on S. E. Zirka, Y. I. Moroz, R. G. Harrison, N. Chisea, "Inverse Hysteresis Models for Transient Simulation,"
      ! IEEE Trans. Power Deliv., Vol. 29, No. 2, pp. 552-559, April 2014.
      subroutine HDHM(ind,bk,h,dhdb,iup,ider)
        use hydatas
        implicit none
        integer d,dd,nmodel,nrev,nrevmax,ider
        double precision bk,h,dhdb,db
        double precision ddhdb,ddhdx,dedx,dh,e,dhout,x,dhrev,dxdb,dbrev
        double precision apar,bpar,cc,beta
        double precision exp,ha,hd,hav,haa,hdd
        integer min,i,ii,j,iup,ind
        integer id(info(1)%nrevmax),bri(info(1)%nrevmax)
        
        nmodel = info(1)%nmodel
        nrevmax = info(1)%nrevmax
        
        if(ind.gt.nmodel) then
          call mexPrintf('Model index too large')
          return
        endif

        ! Keep history
        if(iup.eq.0) then
          hist(nmodel+1)%br(:)      = hist(ind)%br(:)
          hist(nmodel+1)%bt(:)      = hist(ind)%bt(:)
          hist(nmodel+1)%dbout(:)   = hist(ind)%dbout(:)
          hist(nmodel+1)%drev(:)    = hist(ind)%drev(:)
          hist(nmodel+1)%irevs(:,:) = hist(ind)%irevs(:,:)
          hist(nmodel+1)%dp         = hist(ind)%dp
          hist(nmodel+1)%nrev       = hist(ind)%nrev
          hist(nmodel+1)%bp         = hist(ind)%bp
          hist(nmodel+1)%hp         = hist(ind)%hp
        endif
      
        ! B and direction
        if(bk-hist(ind) % bp.gt.1e-6) then
          d=1
        elseif(hist(ind) % bp-bk.gt.1e-6) then
          d=-1
        else
          d=0
        endif

        ! If no change, keep sign
        if(d.eq.0) then
          d = hist(ind) % dp;
        endif
  
        ! Above saturation
        if(abs(bk).ge.hydata(1)%bs) then
          hist(ind) % nrev = 0;
          h         = dsign(1d0,bk)*hav(abs(bk),0);
          if(ider.eq.1) dhdb      = hav(abs(bk),1);
          
        ! Below saturation
        else
      
          ! Check if reversal occurs
      if((abs(hist(ind) % bp).lt.hydata(1)%bs).and.
     &(d*hist(ind) % dp.lt.0).and.(hist(ind) % nrev.lt.nrevmax)) then
  
            hist(ind) % nrev = hist(ind) % nrev+1;
            hist(ind) % drev(hist(ind) % nrev) = -d;
            hist(ind) % br(hist(ind) % nrev) = hist(ind) % bp;

            ! Turns descending
            if(d.lt.0) then
              ! Tip and height of outer loop
              if(hist(ind) % nrev.eq.1) then
                hist(ind) % bt(hist(ind) % nrev) = -hydata(1)%bs
                hist(ind) % dbout(hist(ind) % nrev) = 2*hydata(1)%bs
                hist(ind) % irevs(hist(ind) % nrev,:) = 1
              else
                hist(ind) % bt(hist(ind) % nrev) 
     &        = hist(ind) % br(hist(ind) % nrev-1)
                hist(ind) % dbout(hist(ind) % nrev)
     &        = hist(ind) % br(hist(ind) % nrev)
     &        - hist(ind) % br(hist(ind) % nrev-1)
              hist(ind) % irevs(hist(ind) % nrev,1) = hist(ind) % nrev   ! Turns descending -> Outer ascending curve changes
              hist(ind) % irevs(hist(ind) % nrev,2) = hist(ind) % nrev-1 ! Turns descending -> Outer descending curve same as before
              call update_ha(ind,hist(ind) % nrev)
              !call update_hd(ind,hist(ind) % nrev-1)
              endif

            ! Turns ascending
            else
              ! Tip and height of outer loop
              if(hist(ind) % nrev.eq.1) then
                hist(ind) % bt(hist(ind) % nrev) = hydata(1)%bs
                hist(ind) % dbout(hist(ind) % nrev) = 2*hydata(1)%bs
                hist(ind) % irevs(hist(ind) % nrev,:) = 1
              else
                hist(ind) % bt(hist(ind) % nrev) 
     &        = hist(ind) % br(hist(ind) % nrev-1)
                hist(ind) % dbout(hist(ind) % nrev)
     &        = hist(ind) % br(hist(ind) % nrev-1)
     &        - hist(ind) % br(hist(ind) % nrev)
              hist(ind) % irevs(hist(ind) % nrev,1) = hist(ind) % nrev-1 ! Turns ascending -> Outer ascending curve same as before
              hist(ind) % irevs(hist(ind) % nrev,2) = hist(ind) % nrev   ! Turns ascending -> Outer descending curve changes
              call update_hd(ind,hist(ind) % nrev)
              !call update_ha(ind,hist(ind) % nrev-1)
              endif
            endif
          endif

          ! Check if reversal points are wiped out
          if(hist(ind) % nrev.ge.2) then
            ii=0
            do i=1,hist(ind) % nrev
              if(hist(ind) % drev(i).eq.d) then
                ii=ii+1
                id(ii)=i
              endif
            enddo
            j=nrevmax+1
            if(ii.ge.1) then
              if(d.gt.0) then
                do i=1,ii
                  if(bk.ge.hist(ind) % br(id(i))) j=min(j,i)
                enddo
              else
                do i=1,ii
                  if(bk.le.hist(ind) % br(id(i))) j=min(j,i)
                enddo
              endif
              if(j.lt.nrevmax+1) then
                j=id(j);
                hist(ind) % nrev = j-1;
              endif
            endif;
          endif

          ! On major loop
          if(hist(ind) % nrev.eq.0) then
            if(d.gt.0) then
              h    = ha(ind,0,bk,0);
              if(ider.eq.1) dhdb = ha(ind,0,bk,1);
            else
              h    = hd(ind,0,bk,0);
              if(ider.eq.1) dhdb = hd(ind,0,bk,1);
            endif

          ! On reversal curves
          else

            ! If maxmimum reversals reached, keep last sign
            dd=d
            if(hist(ind) % nrev.eq.nrevmax) then
            dd=-hist(ind) %drev(nrevmax)
            endif
          
            ! Ascending
            if(dd.gt.0) then
              db = hist(ind) % bt(hist(ind) % nrev)-bk
              dbrev = hist(ind) % bt(hist(ind) % nrev)
     &               -hist(ind) % br(hist(ind) % nrev)
              dxdb = -1/dbrev;

            ! Descending
            else
              db = bk-hist(ind) % bt(hist(ind) % nrev)
              dbrev = hist(ind) % br(hist(ind) % nrev)
     &               -hist(ind) % bt(hist(ind) % nrev)
              dxdb = 1/dbrev
            endif
            dhrev = ha(ind,hist(ind) % irevs(hist(ind) % nrev,1),
     &                     hist(ind) % br(hist(ind) % nrev),0)
     &             -hd(ind,hist(ind) % irevs(hist(ind) % nrev,2),
     &                     hist(ind) % br(hist(ind) % nrev),0)
            x = db/dbrev
            beta = dbrev/hist(ind) % dbout(hist(ind) % nrev)
            
            ! Parameters
            apar = dbrev*(7.73+2.76*beta-28.63*beta**2+28.36*beta**3);
            bpar = 0.22*(1-beta); cc = 0.125;
            
            ! H and its derivative
            haa = ha(ind,hist(ind)%irevs(hist(ind)%nrev,1),bk,0)
            hdd = hd(ind,hist(ind)%irevs(hist(ind)%nrev,2),bk,0) 
            dhout = haa-hdd

            e = exp(-apar*(1-x))
            dh = dhrev*(1-bpar)*x*e + dhout*bpar*x**cc
            dedx = apar*exp(-apar*(1-x))
            ddhdx = dhrev*(1-bpar)*(e+x*dedx) + dhout*bpar*cc*x**(cc-1)
            ddhdb = ddhdx*dxdb
      if(dd.gt.0) then
        h    = ha(ind,hist(ind) % irevs(hist(ind) % nrev,1),bk,0)-dh
        if(ider.eq.1) dhdb = 
     &    ha(ind,hist(ind) % irevs(hist(ind) % nrev,1),bk,1)-ddhdb
      else
        h    = hd(ind,hist(ind) % irevs(hist(ind) % nrev,2),bk,0)+dh
        if(ider.eq.1) dhdb = 
     &  hd(ind,hist(ind) % irevs(hist(ind) % nrev,2),bk,1)+ddhdb
        endif
          endif
        endif
        

        hist(ind) % bp = bk;
        hist(ind) % hp = h;
        if(d.ne.0) then
          hist(ind) % dp = d;
        endif

        ! Don't update
        if(iup.eq.0) then
          hist(ind)%br(:)      = hist(nmodel+1)%br(:)
          hist(ind)%bt(:)      = hist(nmodel+1)%bt(:)
          hist(ind)%dbout(:)   = hist(nmodel+1)%dbout(:)
          hist(ind)%drev(:)    = hist(nmodel+1)%drev(:)
          hist(ind)%irevs(:,:) = hist(nmodel+1)%irevs(:,:)
          hist(ind)%dp         = hist(nmodel+1)%dp
          hist(ind)%nrev       = hist(nmodel+1)%nrev
          hist(ind)%bp         = hist(nmodel+1)%bp
          hist(ind)%hp         = hist(nmodel+1)%hp
        endif
        
      return
      end
      
! ##################################################################################################################################################################################################################

      ! Ascending loop function
      double precision function ha(ind,irev,b,ider)
      use hydatas
      implicit double precision(a-h,o-z)
        if(irev.le.1) then
          do i=1,hydata(1)%nm-1
            if(hydata(1)%ba_data(i+1).gt.b) exit
          enddo
          if(i.ge.hydata(1)%nm) i = hydata(1)%nm-1
          db=b-hydata(1)%ba_data(i)
          if(ider.eq.0) then
            ha=hydata(1)%ha_data(i)
     &          +(hydata(1)%ha_data(i+1)-hydata(1)%ha_data(i))
     &          /(hydata(1)%ba_data(i+1)-hydata(1)%ba_data(i))*db
          else
            ha=hydata(1)%dha_data(i)
     &          +(hydata(1)%dha_data(i+1)-hydata(1)%dha_data(i))
     &          /(hydata(1)%ba_data(i+1)-hydata(1)%ba_data(i))*db
          endif
        else
          do i=1,hydata(1)%nm-1
            if(hist(ind)%bas(irev,i+1).gt.b) exit
          enddo
          if(i.ge.hydata(1)%nm) i = hydata(1)%nm-1
          db=b-hist(ind)%bas(irev,i)
          if(ider.eq.0) then
            ha=hist(ind)%has(irev,i)
     &          +(hist(ind)%has(irev,i+1) - hist(ind)%has(irev,i))
     &          /(hist(ind)%bas(irev,i+1) - hist(ind)%bas(irev,i))*db
          else
            ha=hist(ind)%dhas(irev,i)
     &          +(hist(ind)%dhas(irev,i+1)-hist(ind)%dhas(irev,i))
     &          /(hist(ind)%bas(irev,i+1)-hist(ind)%bas(irev,i))*db
          endif
        endif
      return
      end

      
! ##################################################################################################################################################################################################################

      
      ! Descending loop function
      double precision function hd(ind,irev,b,ider)
      use hydatas
      implicit double precision(a-h,o-z)
        if(irev.le.1) then
          hd=-ha(ind,irev,-b,ider)
          if(ider.eq.1) hd=-hd
        else
          do i=1,hydata(1)%nm-1
            if(hist(ind)%bds(irev,i+1).gt.b) exit
          enddo
          if(i.ge.hydata(1)%nm) i = hydata(1)%nm-1
          db=b-hist(ind)%bds(irev,i)
          if(ider.eq.0) then
            hd=hist(ind)%hds(irev,i)
     &          +(hist(ind)%hds(irev,i+1) - hist(ind)%hds(irev,i))
     &          /(hist(ind)%bds(irev,i+1) - hist(ind)%bds(irev,i))*db
          else
            hd=hist(ind)%dhds(irev,i)
     &          +(hist(ind)%dhds(irev,i+1)-hist(ind)%dhds(irev,i))
     &          /(hist(ind)%bds(irev,i+1)-hist(ind)%bds(irev,i))*db
          endif
        endif
      return
      end
      
      
      ! Updating ascending curves
      subroutine update_ha(ind,irev)
      use hydatas
      implicit double precision(a-h,o-z)
        if(irev.gt.info(1)%nrevmax) then
          call mexprintf('ha update error\n')
          stop
        endif
        deltab = hist(ind)%br(irev)-hist(ind)%br(irev-1)
        ddb = deltab/(hydata(1)%nm-1)
        do i = 0,hydata(1)%nm-1
          b = hist(ind)%br(irev-1)+i*ddb
          hist(ind)%bas(irev,i+1) = b
        
          db=hist(ind) % bt(irev-1)-b
          dbrev = hist(ind) % bt(irev-1)-hist(ind) % br(irev-1)
          dxdb = -1/dbrev
          dhrev = ha(ind,hist(ind) % irevs(irev-1,1),
     &                   hist(ind) % br(irev-1),0)
     &           -hd(ind,hist(ind) % irevs(irev-1,2),
     &                   hist(ind) % br(irev-1),0)
          x = db/dbrev;
          beta = dbrev/hist(ind) % dbout(irev-1)
          apar = dbrev*(7.73+2.76*beta-28.63*beta**2+28.36*beta**3)
          bpar = 0.22*(1-beta); cc = 0.125;
          dhout = ha(ind,hist(ind) % irevs(irev-1,1),b,0)
     &           -hd(ind,hist(ind) % irevs(irev-1,2),b,0)
          e = exp(-apar*(1-x))
          dh =  dhrev*(1-bpar)*x*e + dhout*bpar*x**cc
          do ider = 0,1
            if(ider.eq.0) then
              hist(ind)%has(irev,i+1) = ha(ind,hist(ind) %
     %                      irevs(irev-1,1),b,0)-dh
            else
              dedx = apar*exp(-apar*(1-x))
              ddhdx =  dhrev*(1-bpar)*(e+x*dedx)
     &               + dhout*bpar*cc*x**(cc-1)
              ddhdb =  ddhdx*dxdb
              hist(ind)%dhas(irev,i+1) = ha(ind,hist(ind) %
     &                       irevs(irev-1,1),b,1)-ddhdb
            endif
          enddo
        enddo
      return
      end
      
      ! Updating desceding curves
      subroutine update_hd(ind,irev)
      use hydatas
      implicit double precision(a-h,o-z)
        if(irev.gt.info(1)%nrevmax) then
          call mexprintf('hd update error\n')
          stop
        endif

        deltab = hist(ind)%br(irev-1)-hist(ind)%br(irev)
        ddb = deltab/(hydata(1)%nm-1)
        do i = 0,hydata(1)%nm-1
          b = hist(ind)%br(irev)+i*ddb
          hist(ind)%bds(irev,i+1) = b

          db =  b-hist(ind) % bt(irev-1);
          dbrev = hist(ind) % br(irev-1)-hist(ind) % bt(irev-1);
          dxdb = 1/dbrev;
          dhrev = ha(ind,hist(ind) % irevs(irev-1,1),
     &                   hist(ind) % br(irev-1),0)
     &           -hd(ind,hist(ind) % irevs(irev-1,2),
     &                   hist(ind) % br(irev-1),0)
          x = db/dbrev;
          beta = dbrev/hist(ind) % dbout(irev-1);
          apar = dbrev*(7.73+2.76*beta-28.63*beta**2+28.36*beta**3);
          bpar = 0.22*(1-beta); cc = 0.125;
          dhout = ha(ind,hist(ind) % irevs(irev-1,1),b,0)
     &           -hd(ind,hist(ind) % irevs(irev-1,2),b,0)
          e = exp(-apar*(1-x))
          dh =  dhrev*(1-bpar)*x*e + dhout*bpar*x**cc
          do ider = 0,1
            if(ider.eq.0) then
              hist(ind)%hds(irev,i+1) =  hd(ind,hist(ind) %
     &                       irevs(irev-1,2),b,0)+dh;
            else
              dedx  = apar*exp(-apar*(1-x))
              ddhdx = dhrev*(1-bpar)*(e+x*dedx)
     &                + dhout*bpar*cc*x**(cc-1)
              ddhdb = ddhdx*dxdb; 
              hist(ind)%dhds(irev,i+1) = hd(ind,hist(ind) %
     &         irevs(irev-1,2),b,1)+ddhdb;
            endif
          enddo
        enddo
      return
      end
      

      ! Single-valued curve
      double precision function hav(b,ider)
        use hydatas        
        implicit double precision(a-h,o-z)
        do i=1,hydata(1)%nav-1
          if(hydata(1)%bav_data(i+1).gt.b) exit
        enddo
        if(i.ge.hydata(1)%nav) i = hydata(1)%nav-1
        db=b-hydata(1)%bav_data(i)
        if(ider.eq.0) then
          hav=hydata(1)%hav_data(i)+
     &        (hydata(1)%hav_data(i+1)-hydata(1)%hav_data(i))
     &       /(hydata(1)%bav_data(i+1)-hydata(1)%bav_data(i))*db
        else
          hav=hydata(1)%dhav_data(i)+
     &         (hydata(1)%dhav_data(i+1)-hydata(1)%dhav_data(i))
     &        /(hydata(1)%bav_data(i+1)-hydata(1)%bav_data(i))*db
        endif
      return
      end

! ##################################################################################################################################################################################################################
      
	  ! Ascending loop branch as b(h)
      double precision function ba(h)
        use hydatas
        implicit double precision(a-h,o-z)
		
        do i=1,hydata(1)%nm-1
          if(hydata(1)%ha_data(i+1).gt.h) exit
        enddo
        if(i.ge.hydata(1)%nm) i = hydata(1)%nm-1
        dh=h-hydata(1)%ha_data(i)
        ba=hydata(1)%ba_data(i)
     &       +(hydata(1)%ba_data(i+1)-hydata(1)%ba_data(i))
     &       /(hydata(1)%ha_data(i+1)-hydata(1)%ha_data(i))*dh
	     ba=0d0
      return
      end

      ! Deallocation routines
      subroutine deallocate_info
        use hydatas
        if(allocated(info)) deallocate(info)
      return
      end
      subroutine deallocate_hist
        use hydatas
        if(allocated(hist)) then
          do i = 1,size(hist)
            if(allocated(hist(i)%br))    deallocate(hist(i)%br)
            if(allocated(hist(i)%bt))    deallocate(hist(i)%bt)
            if(allocated(hist(i)%dbout)) deallocate(hist(i)%dbout)
            if(allocated(hist(i)%drev))  deallocate(hist(i)%drev)
            if(allocated(hist(i)%irevs)) deallocate(hist(i)%irevs)
            
            if(allocated(hist(i)%has))   deallocate(hist(i)%has)
            if(allocated(hist(i)%hds))   deallocate(hist(i)%hds)
            if(allocated(hist(i)%dhas))  deallocate(hist(i)%dhas)
            if(allocated(hist(i)%dhds))  deallocate(hist(i)%dhds)
            if(allocated(hist(i)%bas))   deallocate(hist(i)%bas)
            if(allocated(hist(i)%bds))   deallocate(hist(i)%bds)
          enddo
          deallocate(hist)
        endif
      return
      end
      subroutine deallocate_hydata
        use hydatas
        if(allocated(hydata)) then
          do i = 1,size(hydata)
            if(allocated(hydata(i)%ba_data))
     &        deallocate(hydata(i)%ba_data)
            if(allocated(hydata(i)%ha_data))
     &        deallocate(hydata(i)%ha_data)
            if(allocated(hydata(i)%dha_data))
     &        deallocate(hydata(i)%dha_data)
            if(allocated(hydata(i)%bav_data))
     &        deallocate(hydata(i)%bav_data)
            if(allocated(hydata(i)%hav_data))
     &        deallocate(hydata(i)%hav_data)
            if(allocated(hydata(i)%dhav_data))
     &        deallocate(hydata(i)%dhav_data)
          enddo
          deallocate(hydata)
        endif
      return
      end
      subroutine deallocate_HDHM
        call deallocate_hist
        call deallocate_hydata
        call deallocate_info
      return
      end
      
      

      ! Routine for initializing hysteresis hydata structure from file
      subroutine initFromFile(nmodel, nrevmax, fname)
        use hydatas
        implicit none
        character(len=200) fname
        integer nmodel,nrevmax,nm,nav,i,ndm,ind
        double precision b,fi,x
      
        ! Number of models
        if(allocated(info)) call deallocate_info
        if(allocated(hydata)) call deallocate_hydata
        allocate(info(1),hydata(1))
        info(1)%nmodel = nmodel
        info(1)%nrevmax = nrevmax
  
        ! Read hydata-structure from file
        open(unit=10,file=trim(fname))
          read(10,'(2i20)') nm,nav
          allocate(hydata(1)%ba_data(nm))
          allocate(hydata(1)%ha_data(nm))
          allocate(hydata(1)%dha_data(nm))
          allocate(hydata(1)%bav_data(nav))
          allocate(hydata(1)%hav_data(nav))
          allocate(hydata(1)%dhav_data(nav))
          read(10,'(<nm>f20.8)') hydata(1)%ba_data(1:nm)
          read(10,'(<nm>f20.8)') hydata(1)%ha_data(1:nm)
          read(10,'(<nm>f20.8)') hydata(1)%dha_data(1:nm)
          read(10,'(<nav>f20.8)') hydata(1)%bav_data(1:nav)
          read(10,'(<nav>f20.8)') hydata(1)%hav_data(1:nav)
          read(10,'(<nav>f20.8)') hydata(1)%dhav_data(1:nav)
        close(unit=10)
        hydata(1)%bs = hydata(1)%ba_data(nm)
        hydata(1)%nm = nm
        hydata(1)%nav = nav
        
        ! Initialize history
        call initHist
                
      return
      end

      ! Routine for initializing history structure
      subroutine initHist
        use hydatas
        implicit none
        integer nrevmax,i,j
        double precision b,x,fi,h

        nrevmax = info(1)%nrevmax

          ! Allocate, deallocate if needed first
          if(allocated(hist)) call deallocate_hist
          allocate(hist(info(1)%nmodel+1))

          ! Allocate space for arrays
          do i = 1,info(1)%nmodel+1
            allocate(hist(i)%br(nrevmax))
            allocate(hist(i)%bt(nrevmax))
            allocate(hist(i)%dbout(nrevmax))
            allocate(hist(i)%drev(nrevmax))
            allocate(hist(i)%irevs(nrevmax,2))

            allocate(hist(i)%has(nrevmax,hydata(1)%nm))
            allocate(hist(i)%hds(nrevmax,hydata(1)%nm))
            allocate(hist(i)%dhas(nrevmax,hydata(1)%nm))
            allocate(hist(i)%dhds(nrevmax,hydata(1)%nm))
            allocate(hist(i)%bas(nrevmax,hydata(1)%nm))
            allocate(hist(i)%bds(nrevmax,hydata(1)%nm))
            
            ! Initial values
            hist(i)%br    = 0
            hist(i)%bt    = 0
            hist(i)%dbout = 0
            hist(i)%drev  = 0
            hist(i)%irevs = 0
            hist(i)%bp    = 0
            hist(i)%hp    = 9e6
            hist(i)%dp    = -1
            hist(i)%nrev  = 0
          enddo
        return
        end