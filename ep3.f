        program ep3F
        implicit none
        integer i, j, k, n, m, NMAX
        integer qrdecomp, lssolve, hue
        parameter (NMAX = 700)
        character fInput*128
        double precision A(NMAX, NMAX), b(NMAX)
        integer p(NMAX)
        
        write (*, *) 'Digite o nome do arquivo:'
        read (*, *) fInput
        open(7, FILE = fInput)
        read (7, *) n, m
        do  k = 1, n*m
            read (7, *) i, j, A(i, j)
        end do
        do k = 1, m
            read (7, *) i, b(i)
        end do
        hue = qrdecomp(n, m, NMAX, A, p)
        stop
        end

        integer function qrdecomp(n, m, lda, A, p)
        implicit none

C       SCALAR ARGUMENTS

        integer n, m, lda

C       ARRAY ARGUMENTS

        double precision A(lda, m)
        integer p(m)

C       LOCAL SCALARS

        double precision innerprod, t, gama
        integer i, j, k, maxv, iMax
        double precision norm(m)
        double precision EPS
        parameter (EPS = 0.0000001)

        qrdecomp = 0

C       CALCULO DAS NORMAS
        do j = 1, m
            maxv = 0
            do i = 1, n
                if (abs(A(i, j)) > maxv) then
                    maxv = A(i, j)
                end if
            end do
            do i = 1, n
                norm(j) = norm(j) + A(i, j)*A(i, j)/(maxv*maxv)
            end do
        end do

C       DECOMPOSICAO QR
        do k = 1, m

C       ACHAR O REFLETOR
            t = 0        
            maxv = 0
            do i = k, n
                if (abs(A(i, k)) > maxv) then
                    maxv = abs(A(i, k))
                end if
            end do
            if (maxv == 0) then
                gama = 0
            else 
                do i = k, n
                    A(i, k) = A(i, k)/maxv
                end do
                do i = k, n
                    t = t + A(i, k)*A(i, k)
                end do
                t = sqrt(t)
                if (A(k, k) < 0) then
                    t = -t
                end if
                A(k, k) = A(k, k) + t
                gama = A(k, k)/t
                do i = (k + 1), n
                    A(i, k) = A(i, k)/A(k, k)
                end do
                A(k, k) = 1
                t = t*maxv
            end if
C       FAZER Q^(K)*A^(K)
        end do
        return
        end
