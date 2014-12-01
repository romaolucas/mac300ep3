        program ep3F
        implicit none
        integer i, j, k, n, m, NMAX
        integer qrdecomp, lssolve, r, fim
        parameter (NMAX = 700)
        character fInput*128
        double precision A(NMAX, NMAX), b(NMAX), gama(NMAX)
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
        r = qrdecomp(n, m, NMAX, A, p, gama)
        write (*, *) 'Pronto'
        write (*, *) 'Posto de A: ', r
        do i = 1, n
            do j = 1, m
                write (*, *) A(i, j)
            end do
        end do

        fim = lssolve(n, m, NMAX, r, A, b, p, gama)
        stop
        end

        integer function qrdecomp(n, m, lda, A, p, gama)
        implicit none

C       SCALAR ARGUMENTS

        integer n, m, lda

C       ARRAY ARGUMENTS

        double precision A(lda, m), gama(m)
        integer p(m)

C       LOCAL SCALARS

        double precision innerprod, t, maxv, aux
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

            p(k) = iMax

C       ACHAR O REFLETOR
            t = 0        
            maxv = 0
            do i = k, n
                if (abs(A(i, k)) > maxv) then
                    maxv = abs(A(i, k))
                end if
            end do
            if (maxv == 0) then
                gama(k) = 0
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
                gama(k) = A(k, k)/t
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
                innerprod = innerprod*gama(k)
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


        integer function lssolve(n, m, lda, rank, A, b, p, gama)
        implicit none

C       SCALAR ARGUMENTS

        integer n, m, lda, rank

C       ARRAY ARGUMENTS

        double precision A(lda, m), gama(m), b(n)
        integer p(m)

C       LOCAL SCALARS

        double precision innerprod, d, aux
        integer i, j, k

        lssolve = 0
C       FAZER Q^(K)*A^(K)
        do k = 1, rank
            innerprod = b(k)
            do i = (k + 1), n
                write (*, *) 'b(i), i, k, A(i, k)'
                write (*, *) b(i), i, k, A(i, k)
                innerprod = innerprod + b(i)*A(i, k)
            end do
            innerprod = innerprod*gama(k)
            write (*, *) innerprod
            b(k) = b(k) - innerprod
            do i = (k + 1), n
                b(i) = b(i) - A(i, k)*innerprod
            end do
            do i = k, n
                write (*, *) b(i)
            end do
        end do

C       TRANSFORMAR A(1:RANK, 1:RANK) EM R_11
        write(*,*), "A(1:rank, 1:rank):"
        do i = 1, rank
            do j = 1, rank
                if (i > j) then
                    A(i, j) = 0
                end if
                write(*,*), A(i,j)
            end do
        end do
        d = 0
        do i = (rank + 1), n
            d = d + b(i) ** 2
            b(i) = 0
        end do
C       RESOLVER O SISTEMA
        do j = rank, 1, -1
            b(j) = b(j) / A(j, j)
            do i = 1, j-1
                b(i) = b(i) - A(i, j) * b(j)
            end do
        end do
        write (*, *) p(1:m)
        do i = m, 1, -1
            aux = b(i)
            b(i) = b(p(i))
            b(p(i)) = aux
        end do
        write(*,*), "x = ", b(1:rank)
        write (*, *) "d = ", d

        return
        end
