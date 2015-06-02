program test2

  real, dimension(10) :: test
  
  test=(/ 1,2,3,4,5,6,7,8,9,10 /)
  write(*,*) test(1:10:2)

end program test2
