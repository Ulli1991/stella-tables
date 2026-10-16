c Driver: real lookup_temp_gamma from calc_temp.F at given (abh2, abhe, ekn).
#include "GIZMO_config.h"
      program test_lookup
      implicit none
#include "cool.h"
      REAL abh2, abhe_in, ekn, temp, gamma
      integer ios
      call init_temperature_lookup
      open(10, file='points.txt', status='old')
      open(11, file='lookup_out.txt', status='replace')
 10   read(10, *, iostat=ios) abh2, abhe_in, ekn
      if (ios .ne. 0) goto 20
      abhe = abhe_in
      call lookup_temp_gamma(abh2, ekn, temp, gamma)
      write(11, '(5es25.16)') abh2, abhe_in, ekn, temp, gamma
      goto 10
 20   close(10)
      close(11)
      end
