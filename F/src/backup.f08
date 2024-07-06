!> @file main.f08
!> @brief This file provides you with the original implementation of pagerank.
!> Your challenge is to optimise it using OpenMP and/or MPI.
!> @author Ludovic Capelli (l.capelli@epcc.ed.ac.uk)

PROGRAM main
    USE omp_lib
    USE mpi_f08

    IMPLICIT NONE

    !> The number of vertices in the graph.
    INTEGER, PARAMETER :: GRAPH_ORDER = 1000
    !> Parameters used in pagerank convergence, do not change.
    REAL(KIND=8), PARAMETER :: DAMPING_FACTOR = 0.85
    !> The number of seconds to not exceed forthe calculation loop.
    INTEGER, PARAMETER :: MAX_TIME = 10
    REAL(KIND=8) :: start
    REAL(KIND=8) :: end
    !> The array in which each vertex pagerank is stored.
    REAL(KIND=8), DIMENSION(0:GRAPH_ORDER-1) :: pagerank
    ! Calculates the sum of all pageranks. It should be 1.0, so it can be used as a quick verification.
    REAL(KIND=8) :: sum_ranks = 0.0
    INTEGER :: i

    !> @brief Indicates which vertices are connected.
    !> @details If an edge links vertex A to vertex B, then adjacency_matrix[A][B]
    !> will be 1.0. The absence of edge is represented with value 0.0.
    !> Redundant edges are still represented with value 1.0.
    REAL(KIND=8), DIMENSION(0:GRAPH_ORDER-1,0:GRAPH_ORDER-1) :: adjacency_matrix

    ! Get the time at the very start.
    start = omp_get_wtime()
    CALL generate_graph()
    CALL calculate_pagerank(pagerank)
    end = omp_get_wtime()

    WRITE(*, '(A,F0.2,A)') 'Total time taken: ', end - start, ' seconds.'
    
    CONTAINS

        SUBROUTINE initialize_graph()
            INTEGER :: i
            INTEGER :: j

            DO i = 0, GRAPH_ORDER - 1
                DO j = 0, GRAPH_ORDER - 1
                    adjacency_matrix(j,i) = 0.0
                END DO
            END DO
            RETURN
        END

    !> @brief Calculates the pagerank of all vertices in the graph.
    !> @param pagerank The array in which store the final pageranks.
    SUBROUTINE calculate_pagerank(pagerank)
        REAL(KIND=8), DIMENSION(0:GRAPH_ORDER-1) :: pagerank
        REAL(KIND=8), DIMENSION(0:GRAPH_ORDER-1) :: new_pagerank
        REAL(KIND=8) :: pagerank_total = 1.0
        REAL(KIND=8) :: initial_rank = 1.0 / GRAPH_ORDER;
        REAL(KIND=8) :: damping_value = (1.0 - DAMPING_FACTOR) / GRAPH_ORDER
        REAL(KIND=8) :: diff = 0.0
        REAL(KIND=8) :: diff_iteration = 0.0
        REAL(KIND=8) :: min_diff = 0.0
        REAL(KIND=8) :: max_diff = 0.0
        INTEGER(KIND=8) :: iteration = 0
        REAL(KIND=8) :: start
        REAL(KIND=8) :: elapsed
        REAL(KIND=8) :: time_per_iteration = 0.0
        REAL(KIND=8) :: iteration_start
        REAL(KIND=8) :: iteration_end
        INTEGER :: outdegree
        INTEGER :: i
        INTEGER :: j
        INTEGER :: k
    
        ! Initialise all vertices to 1/n.
        DO i = 0, GRAPH_ORDER - 1
            pagerank(i) = initial_rank
        END DO
    
        start = omp_get_wtime()
        elapsed = omp_get_wtime() - start
        DO i = 0, GRAPH_ORDER - 1
            new_pagerank(i) = 0.0
        END DO

        ! If we exceeded the MAX_TIME seconds, we stop. If we typically spend X seconds on an iteration, and we are less than X seconds away from MAX_TIME, we stop.
        DO WHILE (elapsed .LT. MAX_TIME .AND. (elapsed + time_per_iteration) .LT. MAX_TIME .AND. pagerank_total .EQ. 1.0)
            iteration_start = omp_get_wtime();
    
            DO i = 0, GRAPH_ORDER - 1
                new_pagerank(i) = 0.0
            END DO
    
            DO i = 0, GRAPH_ORDER - 1
                DO j = 0, GRAPH_ORDER - 1
                    IF (adjacency_matrix(j, i) .EQ. 1.0) THEN
                        outdegree = 0
                        DO k = 0, GRAPH_ORDER - 1
                            IF (adjacency_matrix(k, j) == 1.0) THEN
                                outdegree = outdegree + 1
                            END IF
                        END DO
                        new_pagerank(i) = new_pagerank(i) + pagerank(j) / outdegree
                    END IF
                END DO
            END DO
    
            DO i = 0, GRAPH_ORDER - 1
                new_pagerank(i) = DAMPING_FACTOR * new_pagerank(i) + damping_value
            END DO
    
            diff_iteration = 0.0
            DO i = 0, GRAPH_ORDER - 1
                diff_iteration = diff_iteration + ABS(new_pagerank(i) - pagerank(i))
            END DO
            diff = diff + diff_iteration
            max_diff = MAX(max_diff, diff_iteration)
            min_diff = MIN(max_diff, diff_iteration)
    
            DO i = 0, GRAPH_ORDER - 1
                pagerank(i) = new_pagerank(i)
            END DO
            
            pagerank_total = 0.0
            DO i = 0, GRAPH_ORDER - 1
                pagerank_total = pagerank_total + pagerank(i)
            END DO
            IF (pagerank_total .NE. 1.0) THEN
                WRITE(*, '(A,I0,A,F0.12)') '[ERROR] Iteration ', iteration, ': sum of all pageranks is not 1 but ', pagerank_total
            END IF
    
            iteration_end = omp_get_wtime()
            elapsed = omp_get_wtime() - start
            iteration = iteration + 1
            time_per_iteration = elapsed / iteration
        END DO
        IF (pagerank_total .EQ. 1.0) THEN
            diff = diff / iteration
            WRITE(*, '(I0,A,F0.2,A,F0.6)') iteration, ' iterations achieved in ', elapsed, ' seconds. Average diff = ', diff,', max diff = ', max_diff, ' and min diff = ', min_diff
        END IF
        RETURN
    END

    !> @brief Populates the edges in the graph.
    SUBROUTINE generate_graph()
        IMPLICIT NONE

        INTEGER :: i
        INTEGER :: j
        INTEGER :: source
        INTEGER :: destination
        REAL(KIND=8) :: start

        start = omp_get_wtime()
        CALL initialize_graph()
        DO i = 0, GRAPH_ORDER - 1
            DO j = 0, GRAPH_ORDER - i
                IF (i .NE. j) THEN
                    source = i
                    destination = j
                    adjacency_matrix(destination, source) = 1.0
                END IF
            END DO
        END DO
        WRITE(*, '(F0.2,A)') omp_get_wtime() - start, ' seconds to generate the graph.'
        RETURN
    END
END PROGRAM main