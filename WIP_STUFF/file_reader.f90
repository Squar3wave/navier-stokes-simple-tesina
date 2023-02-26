CHARACTER(128) :: buffer

integer strlen, rows, cols
real, dimension(:,:), allocatable :: x

OPEN (1, file = 'matrix.txt', status='old', action='read')

!Count the number of columns

read(1,'(a)') buffer !read first line WITH SPACES INCLUDED
REWIND(1) !Get back to the file beginning

strlen = len(buffer) !Find the REAL length of a string read
do while (buffer(strlen:strlen) == ' ')
  strlen = strlen - 1
enddo

cols=0 !Count the number of spaces in the first line
do i=0,strlen
  if (buffer(i:i) == ' ') then
    cols=cols+1
  endif
enddo

cols = cols+1

!Count the number of rows

rows = 0 !Count the number of lines in a file
DO
  READ(1,*,iostat=io)
  IF (io/=0) EXIT
  rows = rows + 1
END DO

REWIND(1)

print*, 'Number of rows:', rows
print*, 'Number of columns:', cols

allocate(x(rows,cols))

do I=1,rows,1
  read(1,*) x(I,:)
  write(*,*) x(I,:)
enddo

CLOSE (1)
