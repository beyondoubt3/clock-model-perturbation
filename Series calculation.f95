!This algorithm calculates the ground state energy E（G_ns） of a
!chain with length ns. The model is described in the paper
! 'Perturbative study of the one dimensional quantum clock model'.
!Please run the file 'States generation.nb' first before running
! this one.                         Author:B.Zhang  Date:06/2020
program series_calculation
  implicit integer(a-z)
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer ns,nsta,ord,lendim,temp,lev,i1,i2,i3,maxupdis,maxdndis
  integer,allocatable::v(:,:),states(:,:),state(:),lenlist(:),mytag(:)
  real(dp):: alp,temp1,temp2,pi
  real(dp),allocatable:: elist(:),etot(:),c(:,:)
  open(unit=1,file='States.txt',status='old')
  open(unit=2,file='Lengthlist.txt',status='old')
  open(unit=3,file='Parameters.txt',status='old')
  open(unit=4,file='Result.txt',status='new')
  read(3,*) ns,nsta,ord,alp,lendim,maxupdis,maxdndis
  allocate(lenlist(lendim))
  allocate(elist(nsta))
  allocate(states(ns,nsta))
  allocate(v(nsta,2*ns))
  allocate(etot(ord+1))
  allocate(c(nsta,ord+1))
  allocate(mytag(nsta))

  read(1,*) states
  read(2,*) lenlist

  pi = 3.1415926535897932
  elist=0
  etot=0
  c=0
  c(1,1)=1
  v=0
  temp1=0
  temp2=0
  mytag=0

!Compute the unperturbed energy for every state that will be used
  do i=1,nsta
    state=states(:,i)
    if(ns.eq.1) go to 10
    do j=1,ns-1
        temp1=temp1+2-2*cos((state(j)-state(j+1))*alp)
    end do
10  temp1=temp1+2-2*cos(state(1)*alp)+2-2*cos(state(ns)*alp)
    !The two ends of a cluster are adjacent to '0'.
    elist(i)=temp1
    temp1=0
  end do
  temp1=0
  write(*,*)'elist calculated'



!This part computes the adjacent matrix v between states that are
!in different levels, where the level measures the state's distance
!to the ground state. When N>2, every state is only adjacent to two
!other states, so we can store v with a nsta*2 matrix.
20  Do lev=0,lendim-2
    temp=0
    if(lev.eq.0) go to 30
    do i3=1,lev
        temp=temp+lenlist(i3)
    end do
30  Do i1=temp+1,temp+lenlist(lev+1)
        Do i2=temp+lenlist(lev+1)+1,temp+lenlist(lev+1)+lenlist(lev+2)
             Do i3=1,ns
                temp2=temp2+min(min(abs(states(i3,i1)-states(i3,i2)),&
                 abs(states(i3,i1)-states(i3,i2)-maxupdis-maxdndis-1)),&
                 abs(states(i3,i1)-states(i3,i2)+maxupdis+maxdndis+1))
             End Do
             if(temp2.eq.1) then
                mytag(i1)=mytag(i1)+1
                mytag(i2)=mytag(i2)+1
                v(i1,mytag(i1))=i2
                v(i2,mytag(i2))=i1
             endif
             temp2=0
        End Do
    End Do
    write(*,*) 'Level',lev,'inter-level coupling computed'
  End Do

!This part computes the adjacent matrix elements between states in the
! same level. For example, {1,1,0} is adjacent to {-1,1,0}. When N is
!even, this part can be ignored to speed up the calculation.
    temp=0
    temp2=0
    do lev=2,lendim
    do i3=1,lev-1
        temp=temp+lenlist(i3)
    end do
    Do i1=temp+1,temp+lenlist(lev)
        Do i2=i1,temp+lenlist(lev)
             Do i3=1,ns
                temp2=temp2+min(min(abs(states(i3,i1)-states(i3,i2)),&
                abs(states(i3,i1)-states(i3,i2)-maxupdis-maxdndis-1)),&
                abs(states(i3,i1)-states(i3,i2)+maxupdis+maxdndis+1))
             End Do
             if(temp2.eq.1) then
                mytag(i1)=mytag(i1)+1
                mytag(i2)=mytag(i2)+1
                v(i1,mytag(i1))=i2
                v(i2,mytag(i2))=i1
             endif
             temp2=0
        End Do
    End Do
    temp=0
    write(*,*) 'Level',lev,'intra-level coupling calculated'
    End do



  write(*,*)'v calculated'
  temp2=0



 !This part calculates the ground state energy density series
 ! using iteration
 50 do r=2,ord+1
    do j=1,2*ns
        if(v(1,j).eq.0) go to 60
        etot(r)=etot(r)-c(v(1,j),r-1)
    end do
60  do k=2,nsta
        do j=1,2*ns
            if(v(k,j).eq.0)go to 70
            temp1=temp1-c(v(k,j),r-1)

        end do
70      do s=1,r
            temp2=temp2+etot(r-s+1)*c(k,s)
        end do
        c(k,r)=(-temp1+temp2)/elist(k)
        temp1=0
        temp2=0
    end do
    write(*,*) 'r=',r
  end do



100  write(4,*) 'ns=',ns,'ord=',ord,'N=',2*pi/alp
  write(4,*) 'etot=',etot
  write(*,*) 'ns=',ns,'ord=',ord,'N=',2*pi/alp
  write(*,*) 'etot=',etot!Print the ground state energy density series
  deallocate(elist)
  deallocate(v,etot,c,lenlist,states,mytag)
end
