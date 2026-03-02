PROGRAM particles2D
    use collisions                              ! Import collision module
    implicit none                               ! Enforce explicit typing

    ! --- Constants and Parameters ---
    integer, parameter :: npoints = 1000000     ! Total number of time steps
    integer, parameter :: npar = 5              ! Number of particles in system
    integer, parameter :: ncola = 10            ! Max number of collisions
    double precision :: dt, g, k, xmax, ymax    ! Physics parameters
    double precision :: c, overlap, e           ! Damping and restitution
    double precision :: ddelta1, ddelta2        ! Overlap derivatives
    double precision :: dist, delta             ! Distance and overlap
    double precision :: total_energy            ! Total system energy
    double precision :: kin_e, pot_e            ! Kinetic and potential energy
    
    integer :: i, j, nc, it, p1, p2             ! Loop counters and particle IDs
    integer :: nframes, delay_val               ! Animation frame counters

    ! --- Particle Arrays ---
    double precision :: r(npar), m(npar)        ! Radius and Mass (Fixed properties)
    double precision :: x(npar), y(npar)        ! Position X and Y (Dynamic)
    double precision :: vx(npar), vy(npar)      ! Velocity X and Y (Dynamic)
    double precision :: Fx(npar), Fy(npar)      ! Forces X and Y (Dynamic)
    
    ! --- Tracking Variables ---
    double precision :: deltas(npar)            ! Previous overlap (Memory for damping)

    ! --- Collision Arrays ---
    double precision :: distcol(ncola)          ! Distances for detected collisions
    integer :: cols(ncola, 2)                   ! List of colliding pairs
    
    CHARACTER(LEN=1000) :: cmd                  ! Command string for system calls

    ! --- Setup Output ---
    CALL system("mkdir -p frames")              ! Create folder for images
    OPEN(10, FILE='plot2D.dat', STATUS='REPLACE') ! Open trajectory output file
    OPEN(20, FILE='energy.dat', STATUS='REPLACE') ! Open energy output file
    WRITE(20, *) "# Time(s)    Total_Energy(J)" ! Write energy header

    ! --- Initialization ---
    dt = 1e-6                                   ! Time step [s]
    g = 9.81d0                                  ! Gravity [m/s^2]
    k = 10e4                                    ! Stiffness [N/m]
    e = 0.7                                     ! Coefficient of restitution
    xmax = 0.01d0                               ! Domain width [m]
    ymax = 0.01d0                               ! Domain height [m]
    deltas(:) = 0.0d0                           ! Initialize overlap memory
    nframes = 100                               ! Number of frames for animation

    ! --- Particle Setup ---
    r = (/0.6d0, 0.5d0, 0.8d0, 0.5d0, 0.5d0/)/1000  ! Set radii [m]
    m = 2700 * 1.33 * 3.14 * r**3               ! Calculate mass from density
    write(*,*) "m", m                           ! Print masses to check
    x = (/ 3.0d0, 3.0d0, 5.0d0, 7.0d0, 9.0d0 /)/1000 ! Initial X positions
    y = (/ 5.0d0, 8.0d0, 8.0d0, 8.0d0, 8.0d0 /)/1000 ! Initial Y positions
    vx = (/ 1.0d0, 1.0d0, 2.0d0, 2.5d0, -0.5d0 /)/30 ! Initial X velocities
    vy = 0.0d0                                  ! Initial Y velocities

    ! --- Main Time Loop ---
    do i = 1, npoints                           ! Start integration loop
        
        ! 1. Reset Forces 
        Fy(:) = -m(:) * g                       ! Apply gravity
        Fx(:) = 0.0d0                           ! Reset horizontal force
        
        ! 2. Detect Collisions
        CALL collisionlist(x, y, r, npar, ncola, cols, nc, distcol) ! Find neighbors
        
        ! 3. Resolve Particle-Particle Collisions
        do it = 1, nc                           ! Loop over detected contacts
            p1 = cols(it, 1)                    ! First particle ID
            p2 = cols(it, 2)                    ! Second particle ID
            dist = sqrt(distcol(it))            ! Distance between centers
            delta = r(p1) + r(p2) - dist        ! Calculate overlap
            
            ddelta1 = (delta - deltas(p1))/dt   ! Overlap rate (Part 1)
            ddelta2 = (delta - deltas(p2))/dt   ! Overlap rate (Part 2)
            deltas(p1) = delta                  ! Update memory
            deltas(p2) = delta                  ! Update memory

            ! Calculate Damping Coefficient (Analytic)
            c = -log(e)*2*(k*(1/(1/m(p1)+1/m(p2))))**0.5/3.141
            
            if (dist > 0.000001d0) then         ! Avoid division by zero
                overlap = 1 - dist/(r(p1)+r(p2)) ! Check relative overlap
                if (overlap > 0.01) then        ! Warning if overlap is huge
                    write (*,*) "too much overlap", i, p1, p2, overlap  
                end if
                ! Apply Forces (Spring + Dashpot)
                Fy(p1) = Fy(p1) - (k*delta + c*ddelta1)*(y(p2)-y(p1))/dist
                Fy(p2) = Fy(p2) - (k*delta + c*ddelta2)*(y(p1)-y(p2))/dist
                
                Fx(p1) = Fx(p1) - (k*delta + c*ddelta1)*(x(p2)-x(p1))/dist
                Fx(p2) = Fx(p2) - (k*delta + c*ddelta2)*(x(p1)-x(p2))/dist
            end if
        end do
        
        ! 4. Resolve Wall Collisions
        do it = 1, npar                         ! Loop over all particles
            c = -log(e)*2*(k*m(it))**0.5/3.141  ! Wall damping
            
            ! Bottom Wall (y < r) 
            if (y(it) - r(it) < 0.0d0) then
                Fy(it) = Fy(it) - k*(y(it)-r(it)) - c*vy(it) ! Force UP
            end if
            
            ! Top Wall (y + r > ymax) 
            if (y(it) + r(it) > ymax) then
                Fy(it) = Fy(it) - k*(y(it)+r(it)-ymax) - c*vy(it) ! Force DOWN
            end if
            
            ! Left Wall (x < r) 
            if (x(it) - r(it) < 0.0d0) then
                Fx(it) = Fx(it) - k*(x(it)-r(it)) - c*vx(it) ! Force RIGHT
            end if
            
            ! Right Wall (x + r > xmax) 
            if (x(it) + r(it) > xmax) then
                Fx(it) = Fx(it) - k*(x(it)+r(it)-xmax) - c*vx(it) ! Force LEFT
            end if
        end do
        
        ! 5. Time Integration (Explicit Euler)
        vy(:) = vy(:) + (1.0d0/m(:)) * Fy(:) * dt ! Update velocities
        vx(:) = vx(:) + (1.0d0/m(:)) * Fx(:) * dt ! Update velocities
        
        y(:) = y(:) + dt * vy(:)                ! Update Y positions
        x(:) = x(:) + dt * vx(:)                ! Update X positions

        ! 6. Output Data (Sub-sampling)
        if ((i*nframes)/npoints /= ((i-1)*nframes)/npoints) then
            write(10,*) y, x, Fx, Fy, vy        ! Write trajectory data
            
            ! Calculate Energy
            total_energy = 0.0d0                ! Reset energy
            do j = 1, npar                      ! Loop over particles
                kin_e = 0.5d0 * m(j) * (vx(j)**2 + vy(j)**2) ! Kinetic
                pot_e = m(j) * g * y(j)         ! Potential
                total_energy = total_energy + kin_e + pot_e
            end do
            write(20, *) real(i)*dt, total_energy ! Write energy data
        end if
    end do
    
    ! --- Finalize ---
    CLOSE(10)                                   ! Close trajectory file
    CLOSE(20)                                   ! Close energy file

    ! --- Share Parameters for Gnuplot ---
    open(unit=8, file="params.gp", status="replace") ! Open params file
        write(8,*) "npar=", npar                ! Write N particles
        write(8,*) "xmax=", xmax                ! Write X limit
        write(8,*) "ymax=", ymax                ! Write Y limit
        write(8,*) "dt=", dt                    ! Write Time step
        write(8,*) "npoints=", npoints          ! Write Total points
        write(8,*) "nframes=", nframes          ! Write Total frames
        write(8,*) "array r[npar]"              ! Declare array
    do i = 1, npar                              ! Loop to write radii
        write(8,*) "r[", i, "]=", r(i)          ! Write radius
    end do
    close(8)                                    ! Close params file
    close(10)                                   ! Safety close
    
    ! --- Generate Animation ---
    CALL system('gnuplot generate_frames.gp')   ! Run Gnuplot script
    delay_val = int(nframes * 100 / 1000)       ! Calculate delay
    write(cmd, *) "convert -delay ", delay_val, " frames/frames_*.png animation.gif"
    call system(trim(cmd))                      ! Convert frames to GIF

END PROGRAM particles2D