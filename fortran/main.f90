module hpmodel
    implicit none
    type protein
        logical, allocatable :: hp(:)
        integer, allocatable :: x(:), y(:), local(:)
    end type
end module hpmodel


program main
  use hpmodel
  implicit none

  integer, parameter     :: input_file      = 10
  integer, parameter     :: energy_file     = 20
  integer, parameter     :: trajectory_file = 30
  integer, parameter     :: structure_file  = 40
  integer, parameter     :: log_file        = 50
  real,    parameter     :: kB              = 1.0

  type(protein)            :: chain, pre_structure, min_structure
  integer,   allocatable   :: seed(:), corner_idx(:), crank_idx(:)
  logical                  :: reject
  integer                  :: num_residue, num_local, num_corner, num_crank
  integer                  :: E, Epre, Emin, i, step, timestep, seedsize
  integer                  :: c_rotation, c_corner, c_crank, c_total
  real                     :: T, dT, rnd

  open(log_file, file="hpmodel.log")

  ! read input file
  open(input_file, file="input.dat", action="read")
    call read_file(input_file, chain, num_residue, timestep)
  close(input_file)

  ! generate seed
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  do i=1, seedsize
    call system_clock(count=seed(i))
  end do
  write(log_file,*) "seed size = ", seedsize, ", seed = ", seed(:)
  call random_seed(put=seed)

  ! initialize variable
  T = 1.0
  dT = 0.99995

  num_local  = num_residue - 2
  num_corner = num_local - 2
  num_crank  = num_local - 3
  allocate(corner_idx(num_corner))
  allocate(crank_idx(num_crank))

  chain%local = 0
  chain%x     = 0
  do i=1, size(chain%y)
    chain%y(i) = i-1
  end do
  pre_structure = chain

  E = 0
  Epre = 0
  Emin = 0

  ! main loop
  open(energy_file, file="energy.dat")
  open(trajectory_file, file="trajectory.dat")
  call write_structure(trajectory_file, chain)

  do step = 1, timestep

    ! choose trial move
    call count_rotation(chain, c_rotation)
    call count_corner(chain, c_corner, corner_idx, num_corner)
    call count_crank(chain, c_crank, crank_idx, num_crank)
    c_total = c_rotation + c_corner + c_crank

    call random_number(rnd)
    i = rnd * c_total + 1
    if(i <= c_rotation) then
      call rotation(chain, i)
    else if(i <= c_rotation + c_corner) then
      i = i - c_rotation
      call corner_flip(chain, i, corner_idx, num_corner)
    else if(i <= c_total) then
      i = i - c_rotation - c_corner
      call crank_shaft(chain, i, crank_idx, num_crank)
    else
      write(*, *) "invalid index", i
      stop
    end if

    call update_structure(chain) !< transform local->xy

    call calc_energy(chain, E)

    reject = .false.
    if(E > Epre) then
      if(E == 1) then !< collision
        reject = .true.
      else
        call random_number(rnd)
        reject = (rnd > exp((Epre - E) / (kB * T)))
      endif
    endif

    if(reject) then
      chain = pre_structure
      E = Epre
    else
      pre_structure = chain
      Epre = E

      if(Emin > E) then
        Emin = E
        min_structure = chain
      endif
    end if

    write(energy_file, *) step, E
    call write_structure(trajectory_file, chain)

    T = T * dT
  end do

  write(log_file, *) "T_end = ", T

  write(*, *) "minimum_energy = ", Emin
  open(structure_file, file="min_structure.dat")
    write(structure_file, *) "minimum_energy = ", Emin
    write(structure_file, *) ""
    write(structure_file, *) ""
    call write_structure(structure_file, min_structure)
  close(structure_file)

  close(log_file)
  close(energy_file)
  close(trajectory_file)

  deallocate(chain%x)
  deallocate(chain%y)
  deallocate(chain%local)
  deallocate(chain%hp)
  deallocate(corner_idx)
  deallocate(crank_idx)

end program main

subroutine read_file(input_file, chain, num_residue, total_step)
  use hpmodel
  implicit none

  type(protein), intent(inout) :: chain
  integer,       intent(inout) :: total_step, num_residue
  integer,       intent(in)    :: input_file

  integer                      :: i, j, num_local
  character(len=80)            :: buffer, valstr, seq

  do while(.true.)
    read(input_file, '(a80)', err=100, end=200) buffer

    if(buffer(1:6) .eq. "length") then
        i = index(buffer, "=")

        valstr = buffer(i+1:)
        read(valstr, *) num_residue
        num_local  = num_residue - 2

    else if(buffer(1:10) .eq. "total_step") then
        i = index(buffer, "=")

        valstr = buffer(i+1:)
        read(valstr, *) total_step

    else if(buffer(1:8) .eq. "sequence") then
        i = index(buffer, "H")
        j = index(buffer, "P")
        if(j < i) i = j
        seq = buffer(i:)

    else
        write(*, *) "WARNING: unknown input line: ", buffer
    end if

  end do

100 continue !err
    write(*, *) "file read error"
    stop

200 continue !EOF

    allocate(chain%hp(num_residue))
    allocate(chain%x(num_residue))
    allocate(chain%y(num_residue))
    allocate(chain%local(num_residue - 2))

    do i=1, num_residue
      if (seq(i:i) == 'H') then
        chain%hp(i) = .true.
      else if (seq(i:i) == 'P') then
        chain%hp(i) = .false.
      else
        write(*,*) "invalid sequence: seq(", i, ") = ", seq(i:i)
        stop
      end if
    end do

    return
end subroutine read_file

subroutine rotation(chain, idx)
  use hpmodel
  implicit none
  type(protein), intent(inout) :: chain
  integer,       intent(in)    :: idx
  real                         :: rnd

  if(idx > size(chain%local)) then
    write(*,*) "invalid rotation index: ", idx
    stop
  end if

  call random_number(rnd)
  select case (chain%local(idx))
    case (1)
      if(rnd < 0.5) then
        chain%local(idx) = 0
      else
        chain%local(idx) = -1
      end if
    case (0)
      if(rnd < 0.5) then
        chain%local(idx) = 1
      else
        chain%local(idx) = -1
      end if
    case (-1)
      if(rnd < 0.5) then
        chain%local(idx) = 1
      else
        chain%local(idx) = 0
      end if
    case default
      write(*, *) "invalid local value", chain%local
      stop
  end select

  return
end subroutine rotation

subroutine count_rotation(chain, c)
  use hpmodel
  implicit none

  type(protein), intent(in)    :: chain
  integer,       intent(out)   :: c

  c = size(chain%local)
  return
end subroutine count_rotation

subroutine corner_flip(chain, idx, indices, n)
  use hpmodel
  implicit none
  integer,       intent(in)    :: idx, n
  integer,       intent(in)    :: indices(n)
  type(protein), intent(inout) :: chain
  integer                      :: local_index, i, j, k

  if (idx > n) then
    write(*,*) "invalid corner index: ", idx
    stop
  end if

  local_index = indices(idx)

  if (local_index > size(chain%local)-2 .or. local_index == 0) then
    write(*,*) "invalid corner local index: ", local_index
    stop
  end if

  i = chain%local(local_index)
  j = chain%local(local_index+1)
  k = chain%local(local_index+2)

  chain%local(local_index)   = i+j
  chain%local(local_index+1) = -j
  chain%local(local_index+2) = j+k

  return
end subroutine corner_flip

subroutine count_corner(chain, c, indices, n)
  use hpmodel
  implicit none
  integer,       intent(in)    :: n
  integer,       intent(inout) :: indices(n)
  type(protein), intent(inout) :: chain
  integer,       intent(out)   :: c
  integer                      :: i, num_local

  indices = 0

  c = 0
  do i = 1, n
    if (chain%local(i+1) == 0) then 
      cycle
    else if (abs(chain%local(i)   + chain%local(i+1)) == 2) then
      cycle
    else if (abs(chain%local(i+1) + chain%local(i+2)) == 2) then
      cycle
    else
      c = c + 1
      indices(c) = i
    end if
  end do

  return
end subroutine count_corner

subroutine crank_shaft(chain, idx, indices, n)
  use hpmodel
  implicit none
  integer,       intent(in)    :: idx, n
  integer,       intent(in)    :: indices(n)
  type(protein), intent(inout) :: chain
  integer                      :: local_index

  if (idx > n) then
    write(*,*) "invalid crank index: ", idx
    stop
  end if
  local_index = indices(idx)

  if (local_index > size(chain%local)-3 .or. local_index == 0) then
    write(*,*) "invalid crank local index: ", local_index
    stop
  end if

  chain%local(local_index)   = -chain%local(local_index)
  chain%local(local_index+1) = -chain%local(local_index+1)
  chain%local(local_index+2) = -chain%local(local_index+2)
  chain%local(local_index+3) = -chain%local(local_index+3)

  return
end subroutine crank_shaft

subroutine count_crank(chain, c, indices, n)
  use hpmodel
  implicit none

  integer,       intent(in)    :: n
  type(protein), intent(in)    :: chain
  integer,       intent(inout) :: indices(n)
  integer,       intent(out)   :: c
  integer                      :: i

  indices = 0
  c = 0
  do i = 1, n
    if (chain%local(i)   ==  1 .and. chain%local(i+1) == -1 .and. &
        chain%local(i+2) == -1 .and. chain%local(i+3) ==  1 ) then
      c = c + 1
      indices(c) = i
    else if (chain%local(i)   == -1 .and. chain%local(i+1) ==  1 .and. &
             chain%local(i+2) ==  1 .and. chain%local(i+3) == -1 ) then
      c = c + 1
      indices(c) = i
    end if
  end do
  return
end subroutine count_crank

subroutine update_structure(chain)
  use hpmodel
  implicit none

  type(protein), intent(inout) :: chain
  integer, allocatable         :: vec(:)
  integer                      :: i, num_residue

  num_residue = size(chain%x)
  allocate(vec(num_residue-1))

  vec(1) = 0
  do i=2, num_residue-1
    vec(i) = vec(i-1) + chain%local(i-1)
    if(vec(i) < 0) then
        vec(i) = vec(i) + 4
    else if(vec(i) > 3) then
        vec(i) = vec(i) - 4
    end if
  end do

  do i=2, num_residue
    select case(vec(i-1))
      case (0)
        chain%x(i) = chain%x(i-1)
        chain%y(i) = chain%y(i-1) + 1
      case (1)
        chain%x(i) = chain%x(i-1) + 1
        chain%y(i) = chain%y(i-1)
      case (2)
        chain%x(i) = chain%x(i-1)
        chain%y(i) = chain%y(i-1) - 1
      case (3)
        chain%x(i) = chain%x(i-1) - 1
        chain%y(i) = chain%y(i-1)
      case default
        write(*,*) "invalid vector value: ", vec(i-1)
        stop
    end select
  end do

  return
end subroutine update_structure

subroutine calc_energy(chain, E)
  use hpmodel
  implicit none
  type(protein), intent(in)    :: chain
  integer,       intent(out)   :: E
  integer                      :: num_residue, i, j

  num_residue = size(chain%x)

  E = 0
  collision:do i=1, num_residue - 4
    do j=i+4, num_residue, 2
      if(abs(chain%x(i) - chain%x(j)) + abs(chain%y(i) - chain%y(j)) == 0) then
        E = 1
        exit collision
      end if
    end do
  end do collision
  if (E == 1) return

  do i=1, num_residue-3
    do j=i+3, num_residue, 2
      if (chain%hp(i) .and. chain%hp(j)) then
        if (abs(chain%x(i) - chain%x(j)) + abs(chain%y(i) - chain%y(j)) == 1) then
          E = E-1
        end if
      end if
    end do
  end do

  return
end subroutine calc_energy

subroutine write_structure(trajectory_file, chain)
  use hpmodel
  implicit none
  type(protein), intent(in) :: chain
  integer,       intent(in) :: trajectory_file
  integer                   :: num_residue, i

  num_residue = size(chain%x)
  do i=1, num_residue
    write(trajectory_file, *) chain%x(i), chain%y(i)
  end do
  write(trajectory_file, *) ""
  write(trajectory_file, *) ""

  return
end subroutine write_structure
