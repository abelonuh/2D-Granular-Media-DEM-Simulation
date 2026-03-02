MODULE collisions
    IMPLICIT NONE 

    CONTAINS

    ! ==========================================
    ! Subroutine to detect colliding pairs
    ! ==========================================
    SUBROUTINE collisionlist(x, y, r, npar, ncola, cols, nc, distcol)
        IMPLICIT NONE
        
        ! --- Inputs ---
        INTEGER, INTENT(IN) :: npar, ncola             ! Number of particles, Max collisions allowed
        DOUBLE PRECISION, INTENT(IN) :: x(npar), y(npar) ! Particle positions (X, Y)
        DOUBLE PRECISION, INTENT(IN) :: r(npar)        ! Particle radii

        ! --- Outputs ---
        DOUBLE PRECISION, INTENT(OUT) :: distcol(ncola)! Squared distances of colliding pairs
        INTEGER, INTENT(OUT) :: cols(ncola, 2)         ! List of colliding pairs (ID1, ID2)
        INTEGER, INTENT(OUT) :: nc                     ! Number of detected collisions

        ! --- Local Variables ---
        INTEGER :: i, j                                ! Loop counters
        DOUBLE PRECISION :: dx, dy, dist2              ! Distance components and squared distance
        DOUBLE PRECISION :: rsum                       ! Sum of radii

        ! --- Initialization ---
        cols = 0                                       ! Reset collision list
        nc = 0                                         ! Reset collision counter

        ! --- Loop over particle pairs ---
        DO i = 1, npar - 1                             ! Loop over first particle
            DO j = i + 1, npar                         ! Loop over second particle (unique pairs)
                
                ! Calculate distance
                dx = x(i) - x(j)                       ! X-distance
                dy = y(i) - y(j)                       ! Y-distance
                rsum = r(i) + r(j)                     ! Sum of radii (contact threshold)
                dist2 = dx*dx + dy*dy                  ! Squared distance (avoids sqrt for speed)

                ! Check for overlap
                IF (dist2 <= rsum*rsum) THEN           ! If distance < sum of radii -> Contact
                    
                    ! Store collision if space exists
                    IF (nc < ncola) THEN               
                        nc = nc + 1                    ! Increment counter
                        cols(nc, 1) = i                ! Store ID of first particle
                        cols(nc, 2) = j                ! Store ID of second particle
                        distcol(nc) = dist2            ! Store squared distance
                    END IF
                    
                END IF
            END DO
        END DO

    END SUBROUTINE collisionlist

END MODULE collisions