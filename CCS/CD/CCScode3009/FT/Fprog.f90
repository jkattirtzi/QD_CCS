Program FTtest
use FT1
Implicit None
character(Len=10)::Line1
real(kind=8),dimension(:),allocatable::tarray
complex(kind=8),dimension(:),allocatable::ACFarray
! Time variables
real(kind=8)::starttime, endtime
!character(LEN=8)::date
character(LEN=10)::time
CALL CPU_TIME(starttime)
Print *, 'Starting Fprog'
Print *, 'calling ReadFT1'
call ReadFT1(Line1)
Print *, 'done...calling ReadFT2'
call ReadFT2(Line1,tarray,ACFarray)
Print *, 'done...calling FTRtn'
call FTRtn(tarray,acfarray)
Print *,'everything worked ok...well done!'
CALL CPU_TIME(endtime)
Print *, 'Time Taken:', (endtime-starttime)/60.0d0, 'mins'
End Program FTtest
