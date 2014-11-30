        program ep3F
        implicit none
        integer i, j, k, n, m, NMAX
        integer qrdecomp, lssolve, r
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
        do k = 1, n
            read (7, *) i, b(i)
        end do
        r = qrdecomp(n, m, NMAX, A, p)
        write (*, *) 'Pronto'
        write (*, *) 'Posto de A: ', r
        do i = 1, n
            do j = 1, m
                write (*, *) A(i, j)
            end do
        end do
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

        double precision innerprod, t, gama, maxv, aux
        integer i, j, k, iMax
        double precision norm(m)
        double precision EPS
        parameter (EPS = 0.0000001)

        qrdecomp = m 

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

C       PERMUTACOES
            do j = k , m
                maxv = 0
                if (abs(norm(j)) > maxv) then
                    maxv = abs(norm(j))
                    iMax = j
                end if
            end do

            do i = 1, n
                aux = A(i, k)
                A(i, k) = A(i, iMax)
                A(i, iMax) = aux
            end do

            aux = norm(k)
            norm(k) = norm(iMax)
            norm(iMax) = aux

            aux = p(k)
            p(k) = p(iMax)
            p(iMax) = aux

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
                if (abs(t*maxv) < EPS) then
                    do i = k, n
                        A(i, k) = A(i, k)*maxv
                    end do
                    qrdecomp = k - 1
                    exit
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
            do j = (k + 1), m
                innerprod = 0
                do i = k, n
                    write (*, *) 'A(i, j), i, j, k, A(i, k)'
                    write (*, *) A(i, j), i, j, k, A(i, k)
                    innerprod = innerprod + A(i, j)*A(i, k)
                end do
                innerprod = innerprod*gama
                write (*, *) innerprod
                do i = k, n
                    A(i, j) = A(i, j) - A(i, k)*innerprod
                end do
                do i = k, n
                    write (*, *) A(i , j)
                end do
            end do
            A(k, k) = -t
C       NORMAS ATUALIZADAS
            do j = (k + 1), m
                norm(j) = norm(j) - A(k, j) ** 2
            end do
        end do
        return
        end
