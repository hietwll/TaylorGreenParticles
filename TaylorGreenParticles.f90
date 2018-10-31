module taylor
implicit none
integer :: ng, np
real,allocatable :: u(:,:),v(:,:),om(:,:)
real :: dt,rhop,dp,s,vis
real :: pi=3.141592653,dx,t
real,allocatable :: x(:),y(:),px(:),py(:),pu(:),pv(:)
integer iend,isave,iprint
contains

subroutine init()
    implicit none
    integer i,j
    real tmp


    open(unit=133,file='input.txt')
        read(133,*) ng
        read(133,*) np
        read(133,*) iend
        read(133,*) isave
        read(133,*) iprint
        read(133,*) dt
        read(133,*) vis
        read(133,*) rhop
        read(133,*) dp
    close(133)


    write(*,*) 'read input ok!'
    allocate(u(ng,ng));allocate(v(ng,ng));allocate(om(ng,ng))
    allocate(x(ng));allocate(y(ng))
    allocate(px(np));allocate(py(np))
    allocate(pu(np));allocate(pv(np))
    
    write(*,*) 'allocate ok'

    dx = 2*pi/real(ng-1)

    do i = 1,ng
        x(i) = (i-1)*dx
        y(i) = (i-1)*dx
    enddo

    do j = 1,ng
        do i = 1, ng
            u(i,j) = sin(x(i))*cos(y(j))
            v(i,j) = -cos(x(i))*sin(y(j))
            om(i,j) = sin(x(i))*sin(y(j))
        enddo
    enddo

    call random_seed()

    do i = 1, np
        call random_number(tmp)
        px(i) = tmp*2.0*pi
        call random_number(tmp)
        py(i) = tmp*2.0*pi
    enddo

    pu = 0.0
    pv = 0.0

    write(*,*) 'Init ok'
end subroutine init

subroutine finalize()
    implicit none

    deallocate(u);deallocate(v);deallocate(om)
    deallocate(x);deallocate(y)
    deallocate(px);deallocate(py)
    deallocate(pu);deallocate(pv)


end subroutine finalize

subroutine update_fluid()
    implicit none
    integer i,j
    do i = 1, ng
        do j = 1 , ng
            u(i,j) = sin(x(i))*cos(y(j))*exp(-2.0*vis*t)
            v(i,j) = -cos(x(i))*sin(y(j))*exp(-2.0*vis*t)
            om(i,j) = sin(x(i))*sin(y(j))*exp(-2.0*vis*t)
        enddo
    enddo
end subroutine update_fluid

subroutine update_particle()

    implicit none
    integer i 
    real u_f,v_f,beta,rep,taop,puo,pvo

    do i = 1, np
        u_f = sin(px(i))*cos(py(i))*exp(-2.0*vis*t)
        v_f = -cos(px(i))*sin(py(i))*exp(-2.0*vis*t)
        rep = sqrt((u_f-pu(i))**2+(v_f-pv(i))**2)*dp/vis   
        beta = 1+0.15*rep**0.687
        taop = rhop*dp**2/(18.0*vis)
        puo = pu(i)
        pvo = pv(i)
        pu(i) = pu(i)+dt*beta*(u_f-pu(i))/taop
        pv(i) = pv(i)+dt*beta*(v_f-pv(i))/taop
        px(i) = px(i)+0.5*dt*(pu(i)+puo)
        py(i) = py(i)+0.5*dt*(pv(i)+pvo)
        if(px(i)>2*pi) then
            px(i)=px(i)-2*pi
        elseif(px(i)<0) then
            px(i)=px(i)+2*pi
        endif
        if(py(i)>2*pi) then
            py(i)=py(i)-2*pi
        elseif(py(i)<0) then
            py(i)=py(i)+2*pi
        endif
    enddo
end subroutine


subroutine output(step)
    implicit none 
    integer i,j
    character*100 filename1,filename2
    character*14 zone1,zone2
    character*8 charstep
    integer step 

    write(charstep,'(I8.8)') step

    filename1 = charstep//'fluid_'//'.dat'
    filename2 = charstep//'particle_'//'.csv'
    zone1 = charstep//'fluid_'
    zone2 = charstep//'parti_'



    ! write fluid data
    open(unit=101,file=filename1)
    
    write(101,*) 'VARIABLES = "x","y","u","v","om"'
    write(101,"('Zone T=""',A14,'"" I=',I8,',J=',I8,',F=POINT')") zone1,ng,ng

    do i = 1, ng
        do j = 1, ng
            write(101,"(5F15.7)") x(i),y(j),u(i,j),v(i,j),om(i,j)
        enddo
    enddo

    close(101)


    ! write particle data
    open(unit=102,file=filename2)
    write(102,*) 'x,y,u,v,dp'!paraviewÖÐparticle¸ñÊ½
    !write(102,*) 'VARIABLES = "x","y","u","v","dp"'
    !write(102,"('Zone T=""',A14,'"" I=',I8,',F=POINT')") zone2,np

    do i = 1, np
        write(102,"(F15.7,',',F15.7,',',F15.7,',',F15.7,',',F15.7)") px(i),py(i),pu(i),pv(i),dp    
    enddo

    close(102)



end subroutine output

end module


program main
    use taylor
    implicit none
    integer istep

    call init()

    do istep = 1,iend,1
        t = dt*istep
        call update_fluid()
        !write(*,*) 'update fluid ok'
        call update_particle()
        !write(*,*) 'update par ok'
        if(mod(istep,iprint).eq.0) then
          write(*,*) 'current step ', istep
        endif

        if(mod(istep,isave).eq.0) then
            write(*,*) 'writing data ',istep
            call output(istep)
        endif
    enddo

    call finalize()

end program
