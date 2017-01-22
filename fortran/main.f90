program hpmodel
  implicit none

  type protein
    logical, allocatable :: hp(:)
    integer, allocatable :: x(:), y(:), local(:)
  end type protein

  integer, parameter     :: input_file      = 10
  integer, parameter     :: energy_file     = 20
  integer, parameter     :: trajectory_file = 30
  integer, parameter     :: structure_file  = 40
  integer, parameter     :: log_file        = 50
  type(protein)          :: chain, min_structure
  character, allocatable :: seq(:)
  integer  , allocatable :: seed(:)
  logical                :: reject
  integer                :: num_residue, timestep, seedsize
  integer                :: trial_move, random_index, direction
  integer                :: E, Epre, Emin, i, step
  real, parameter        :: kB = 1.0
  real                   :: T, dT, rnd

  open(log_file, file="hpmodel.log")

  ! read input file
  open(input_file, file="input.dat", action="read")
    read(input_file, '(i)') num_residue
    write(log_file, *) "num_residue = ", num_residue

    allocate(chain%hp(num_residue))
    allocate(chain%x(num_residue))
    allocate(chain%y(num_residue))
    allocate(chain%local(num_residue-2))

    read(input_file, '(i)') timestep
    write(log_file, *) "step = ", timestep

    allocate(seq(num_residue))
    read(input_file, '(a)') seq
    write(log_file, *) "seq  = ", seq(1:num_residue)

    do i=1, num_residue
      if (seq(i) == 'H') then
        chain%hp(i) = .true.
      else if (seq(i) == 'P') then
        chain%hp(i) = .false.
      else
        write(*,*) "invalid sequence: ", seq(i)
        stop
      end if
    end do

    deallocate(seq)
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

  chain%local = 0
  do i=1, num_residue
    chain%x(i) = 0
    chain%y(i) = i-1
  end do

  E = 0
  Epre = 0
  Emin = 0

  ! main loop
  open(energy_file, file="energy.dat")
  open(trajectory_file, file="trajectory.dat")
  call write_structure(trajectory_file, chain)

  do step = 1, timestep
    ! trial move
!     call random_number(rnd)
!     trial_move = rnd * 3
!     select case(trial_move)
!       case (0)
        call rotation(chain, random_index, direction)
!       case (1)
!         call corner_flip(chain, random_index)
!       case (2)
!         call crank_shaft(chain, random_index)
!       case default
!         write(*, *) "invalid random number"
!         stop
!     end select
    call update_structure(chain) !< transform local->xy

    ! Metropolis-Hastings update
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
!         select case(trial_move)
!           case (0)
      call restore_rotation(chain, random_index, direction)
!           case (1)
!             call restore_corner_flip(chain, random_index)
!           case (2)
!             call restore_crank_shaft(chain, random_index)
!           case default
!             write(*, *) "invalid trial move"
!             stop
!         end select
      call update_structure(chain)
      E = Epre
    else
      Epre = E !< update
    end if

    if(Emin > E) then
      Emin = E
      min_structure = chain
    endif

    write(energy_file, *), step, E
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

end program hpmodel

subroutine rotation(chain, random_index, direction)
  type protein
    logical, allocatable :: hp(:)
    integer, allocatable :: x(:), y(:), local(:)
  end type protein

  type(protein), intent(inout) :: chain
  integer, intent(out)         :: random_index, direction
  integer                      :: num_local, trial
  real                         :: rnd

  num_local = size(chain%local)

  call random_number(rnd)
  random_index = rnd * num_local + 1
  direction = chain%local(random_index)

  call random_number(rnd)
  select case (direction)
    case (1)
      if(rnd < 0.5) then
        chain%local(random_index) = 0
      else
        chain%local(random_index) = -1
      end if
    case (0)
      if(rnd < 0.5) then
        chain%local(random_index) = 1
      else
        chain%local(random_index) = -1
      end if
    case (-1)
      if(rnd < 0.5) then
        chain%local(random_index) = 1
      else
        chain%local(random_index) = 0
      end if
    case default
      write(*, *) "invalid local value", direction
      stop
  end select

  return
end subroutine rotation

subroutine restore_rotation(chain, random_index, direction)
  type protein
    logical, allocatable :: hp(:)
    integer, allocatable :: x(:), y(:), local(:)
  end type protein

  type(protein), intent(inout) :: chain
  integer, intent(in)          :: random_index, direction
  
  chain%local(random_index) = direction

  return
end subroutine restore_rotation

subroutine update_structure(chain)
  type protein
    logical, allocatable :: hp(:)
    integer, allocatable :: x(:), y(:), local(:)
  end type protein

  type(protein), intent(inout) :: chain
  integer, allocatable :: vec(:)
  integer              :: num_residue, sum_local

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
  type protein
    logical, allocatable :: hp(:)
    integer, allocatable :: x(:), y(:), local(:)
  end type protein

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
  type protein
    logical, allocatable :: hp(:)
    integer, allocatable :: x(:), y(:), local(:)
  end type protein

  type(protein), intent(in) :: chain
  integer,       intent(in) :: trajectory_file
  integer                   :: num_residue, i

  num_residue = size(chain%x)
  do i=1, num_residue
    write(trajectory_file, *), chain%x(i), chain%y(i)
  end do
  write(trajectory_file, *), ""
  write(trajectory_file, *), ""

  return
end subroutine write_structure
