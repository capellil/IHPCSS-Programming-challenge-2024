/**
 * @file main.f08
 * @brief This file provides you with the original implementation of pagerank.
 * Your challenge is to optimise it using OpenMP and/or MPI.
 * @author Ludovic Capelli (l.capelli@epcc.ed.ac.uk)
 **/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

/// The number of vertices in the graph.
#define GRAPH_ORDER 1000
/// Parameters used in pagerank convergence, do not change.
#define DAMPING_FACTOR 0.85
/// The number of seconds to not exceed forthe calculation loop.
#define MAX_TIME 10

/**
 * @brief Indicates which vertices are connected.
 * @details If an edge links vertex A to vertex B, then adjacency_matrix[A][B]
 * will be 1.0. The absence of edge is represented with value 0.0.
 * Redundant edges are still represented with value 1.0.
 */
double adjacency_matrix[GRAPH_ORDER][GRAPH_ORDER];
double max_diff = 0.0;
double min_diff = 1.0;
double total_diff = 0.0;
 
void initialize_graph(void)
{
    for(int i = 0; i < GRAPH_ORDER; i++)
    {
        for(int j = 0; j < GRAPH_ORDER; j++)
        {
            adjacency_matrix[i][j] = 0.0;
        }
    }
}

/**
 * @brief Calculates the pagerank of all vertices in the graph.
 * @param pagerank The array in which store the final pageranks.
 */
void calculate_pagerank(double pagerank[])
{
    double initial_rank = 1.0 / GRAPH_ORDER;
 
    // Initialise all vertices to 1/n.
    for(int i = 0; i < GRAPH_ORDER; i++)
    {
        pagerank[i] = initial_rank;
    }
 
    double damping_value = (1.0 - DAMPING_FACTOR) / GRAPH_ORDER;
    double diff = 1.0;
    size_t iteration = 0;
    double start = omp_get_wtime();
    double elapsed = omp_get_wtime() - start;
    double time_per_iteration = 0;
    double new_pagerank[GRAPH_ORDER];
    for(int i = 0; i < GRAPH_ORDER; i++)
    {
        new_pagerank[i] = 0.0;
    }

    // If we exceeded the MAX_TIME seconds, we stop. If we typically spend X seconds on an iteration, and we are less than X seconds away from MAX_TIME, we stop.
    while(elapsed < MAX_TIME && (elapsed + time_per_iteration) < MAX_TIME)
    {
        double iteration_start = omp_get_wtime();
 
        for(int i = 0; i < GRAPH_ORDER; i++)
        {
            new_pagerank[i] = 0.0;
        }
 
		for(int i = 0; i < GRAPH_ORDER; i++)
        {
			for(int j = 0; j < GRAPH_ORDER; j++)
            {
				if (adjacency_matrix[j][i] == 1.0)
                {
					int outdegree = 0;
				 
					for(int k = 0; k < GRAPH_ORDER; k++)
                    {
						if (adjacency_matrix[j][k] == 1.0)
                        {
							outdegree++;
						}
					}
					new_pagerank[i] += pagerank[j] / (double)outdegree;
				}
			}
		}
 
        for(int i = 0; i < GRAPH_ORDER; i++)
        {
            new_pagerank[i] = DAMPING_FACTOR * new_pagerank[i] + damping_value;
        }
 
        diff = 0.0;
        for(int i = 0; i < GRAPH_ORDER; i++)
        {
            diff += fabs(new_pagerank[i] - pagerank[i]);
        }
        max_diff = (max_diff < diff) ? diff : max_diff;
        total_diff += diff;
        min_diff = (min_diff > diff) ? diff : min_diff;
 
        for(int i = 0; i < GRAPH_ORDER; i++)
        {
            pagerank[i] = new_pagerank[i];
        }
            
        double pagerank_total = 0.0;
        for(int i = 0; i < GRAPH_ORDER; i++)
        {
            pagerank_total += pagerank[i];
        }
        if(fabs(pagerank_total - 1.0) >= 1E-12)
        {
            printf("[ERROR] Iteration %zu: sum of all pageranks is not 1 but %.12f.\n", iteration, pagerank_total);
        }
 
		double iteration_end = omp_get_wtime();
		elapsed = omp_get_wtime() - start;
		iteration++;
		time_per_iteration = elapsed / iteration;
    }
    
    printf("%zu iterations achieved in %.2f seconds\n", iteration, elapsed);
}

/**
 * @brief Populates the edges in the graph for testing.
 **/
void generate_nice_graph(void)
{
    printf("Generate a graph for testing purposes (i.e.: a nice and conveniently designed graph :) )\n");
    double start = omp_get_wtime();
    initialize_graph();
    for(int i = 0; i < GRAPH_ORDER; i++)
    {
        for(int j = 0; j < GRAPH_ORDER; j++)
        {
            int source = i;
            int destination = j;
            if(i != j)
            {
                adjacency_matrix[source][destination] = 1.0;
            }
        }
    }
    printf("%.2f seconds to generate the graph.\n", omp_get_wtime() - start);
}

/**
 * @brief Populates the edges in the graph for the challenge.
 **/
void generate_sneaky_graph(void)
{
    printf("Generate a graph for the challenge (i.e.: a sneaky graph :P )\n");
    double start = omp_get_wtime();
    initialize_graph();
    for(int i = 0; i < GRAPH_ORDER; i++)
    {
        for(int j = 0; j < GRAPH_ORDER - i; j++)
        {
            int source = i;
            int destination = j;
            if(i != j)
            {
                adjacency_matrix[source][destination] = 1.0;
            }
        }
    }
    printf("%.2f seconds to generate the graph.\n", omp_get_wtime() - start);
}

int main(int argc, char* argv[])
{
    // We do not need argc, this line silences potential compilation warnings.
    (void) argc;
    // We do not need argv, this line silences potential compilation warnings.
    (void) argv;

    printf("This program has two graph generators: generate_nice_graph and generate_sneaky_graph. If you intend to submit, your code will be timed on the sneaky graph, remember to try both.\n");

    // Get the time at the very start.
    double start = omp_get_wtime();
    
    generate_nice_graph();
 
    /// The array in which each vertex pagerank is stored.
    double pagerank[GRAPH_ORDER];
    calculate_pagerank(pagerank);
 
    // Calculates the sum of all pageranks. It should be 1.0, so it can be used as a quick verification.
    double sum_ranks = 0.0;
    for(int i = 0; i < GRAPH_ORDER; i++)
    {
        if(i % 100 == 0)
        {
            printf("PageRank of vertex %d: %.6f\n", i, pagerank[i]);
        }
        sum_ranks += pagerank[i];
    }
    printf("Sum of all pageranks = %.12f, total diff = %.12f, max diff = %.12f and min diff = %.12f.\n", sum_ranks, total_diff, max_diff, min_diff);
    double end = omp_get_wtime();
 
    printf("Total time taken: %.2f seconds.\n", end - start);
 
    return 0;
}
