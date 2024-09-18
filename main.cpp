#include <iostream>
#include <chrono>
#include <omp.h>
#include "Maze.cpp"
#include "StatsUtil.cpp"


using namespace std;
using namespace chrono;

#define FIXED_SEED 5
#define NUMBER_TEST 25


void sequential_solver(Maze& maze, vector<Particle>& particles);

void parallel_solver(Maze& maze, vector<Particle>& particles, int threads_number);

vector<Particle> get_particles(Point start, int num_pariceles);

vector<Maze> get_mazes(vector<pair<int, int>> maze_dimensions);


int main()
{
    srand(FIXED_SEED);

    vector<int> NUM_PARTICLES = {100, 250, 500, 1000, 1500, 2000};
    vector<int> NUM_THREADS = {2, 4, 6, 8, 16, 32, 64, 128};
    vector<pair<int, int>> MAZE_DIM = {{5, 5}, {10, 10}, {20, 20}, {30, 30}, {40, 40}, {50, 50}};

    //vector<int> NUM_PARTICLES = {500};
    //vector<pair<int, int>> MAZE_DIM = {{15, 15}};
    //vector<int> NUM_THREADS = {6, 8};


    vector<Maze> mazes = get_mazes(MAZE_DIM);

    StatsUtil stats_util;
    cout << "------------------------------------" << endl
        << "-------- SEQUENTIAL VERSION --------" << endl
        << "------------------------------------" << endl;
    for (auto original_maze : mazes)
    {
        for (auto num_par : NUM_PARTICLES)
        {
            cout << num_par << " particles" << endl;
            for (int i = 0; i < NUMBER_TEST; i++)
            {
                vector<Particle> particles = get_particles(original_maze.start, num_par);
                Maze maze = original_maze;

                double start_time = omp_get_wtime();
                sequential_solver(maze, particles);
                //auto seq_time_elapsed = duration_cast<milliseconds>(omp_get_wtime() - start_time).count();
                auto seq_time_elapsed = (omp_get_wtime() - start_time) * 1000;

                stats_util.addSequentialTime(original_maze.m, num_par, seq_time_elapsed);
            }
            stats_util.printSummarySequentialStats({original_maze.m, num_par});
        }
    }
    cout << "-------- END SEQUENTIAL VERSION --------\n" << endl << endl;


    cout << "-------------------------------------" << endl
        << "--------- PARALLEL VERSION ----------" << endl
        << "-------------------------------------" << endl;
    for (auto original_maze : mazes)
    {
        for (auto num_par : NUM_PARTICLES)
        {
            cout << num_par << " particles" << endl;
            for (auto num_thread : NUM_THREADS)
            {
                for (int i = 0; i < NUMBER_TEST; i++)
                {
                    vector<Particle> particles = get_particles(original_maze.start, num_par);
                    Maze maze = original_maze;

                    auto start_time = omp_get_wtime();
                    parallel_solver(maze, particles, num_thread);
                    //auto par_time_elapsed = duration_cast<milliseconds>(system_clock::now() - start_time).count();
                    auto par_time_elapsed = (omp_get_wtime() - start_time) * 1000;
                    
                    stats_util.addParallelTime(original_maze.m, num_par, num_thread, par_time_elapsed);
                }
                stats_util.printSummaryParallelStats({original_maze.m, num_par, num_thread});
            }
        }
    }
    cout << "-------- END PARALLEL VERSION --------\n" << endl;

    stats_util.printAllStatistics();

    stats_util.exportStatisticsToCSV("C:/Users/plato/Downloads/maze_solver_statistics.csv");


    git config --global user.name "PlatRama"
    git config --global user.email "plator97@gmail.com"


    return 0;
}

/**
 * @brief Solves the maze by moving particles sequentially until one finds the exit.
 *
 * This function iteratively moves each particle in the maze until one of them finds the exit.
 * When a solution is found, the path is drawn on the maze.
 *
 * @param maze The maze to be solved.
 * @param particles A vector of particles that attempt to solve the maze.
 * @return A pointer to the maze array with the solution path drawn, if found.
 */
void sequential_solver(Maze& maze, vector<Particle>& particles)
{
    Particle* winning_particle = nullptr;

    // Main loop: continue until a solution is found
    while (!winning_particle)
    {
        for (Particle& particle : particles)
        {
            if (winning_particle)
                continue;

            maze.randomMove(particle);

            // Check if the current particle has found the exit
            if (maze.isExitFound(particle.currentPosition))
                winning_particle = &particle;
        }
    }

    // Draw the solution if a winning particle was found
    if (winning_particle)
        maze.drawPath(*winning_particle);
}

/**
 * @brief Solves the maze using multiple particles in parallel.
 *
 * This method uses OpenMP to distribute the particles among several threads.
 * Each thread processes a subset of particles, moving them until one finds the exit.
 * Once a particle finds the exit, the solution is drawn in the maze.
 *
 * @param maze The maze in which the particles are moving.
 * @param particles A vector of particles that attempt to solve the maze.
 * @param threads_number The number of threads to use for parallel processing.
 */
void parallel_solver(Maze& maze, vector<Particle>& particles, int threads_number)
{
    int partition_size = particles.size() / threads_number;
    Particle* winning_particle = nullptr;
    bool exit_found = false;

#pragma omp parallel num_threads(threads_number) shared(particles, maze, partition_size)
    {
        int threadId = omp_get_thread_num(); // Get the current thread ID
        int start_partition = partition_size * threadId; // Calculate the start index for this thread's partition
        int end_partition = (threadId == threads_number - 1) ? particles.size() : start_partition + partition_size;
        // Calculate the end index for this partition

        // Main loop: continues until a particle finds the exit
        while (!exit_found)
        {
            for (int i = start_partition; i < end_partition; i++)
            {
                // Skip iteration if a winning particle is already found
                if (exit_found) continue;

                // Move the current particle in a random direction
                maze.randomMove(particles[i]);

                // Check if the current particle has found the exit
                if (maze.isExitFound(particles[i].currentPosition))
                {
                    // Critical section: ensures only one thread can assign winning_particle
#pragma omp critical
                    {
                        if (!exit_found)
                        {
                            // Additional check to avoid multiple threads assigning the pointer
                            winning_particle = &particles[i];
                            exit_found = true;
                        }
                    }
                }
            }
        }
    }

    // Draw the solution in the maze if a winning particle was found
    if (winning_particle)
        maze.drawPath(*winning_particle);
}

vector<Particle> get_particles(Point start, int num_pariceles)
{
    vector<Particle> particles;
    for (int i = 0; i < num_pariceles; i++)
    {
        Particle p = {start};
        p.path.push_back(p.currentPosition);
        particles.push_back(p);
    }
    return particles;
}

vector<Maze> get_mazes(vector<pair<int, int>> maze_dimensions)
{
    vector<Maze> mazes;
    for (auto dim : maze_dimensions)
    {
        mazes.push_back(Maze(dim.first, dim.second));
    }
    return mazes;
}
