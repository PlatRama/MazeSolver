//
// Created by plato on 11/09/2024.
//
#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <tuple>
#include <fstream>

using namespace std;

// Structure for sequential configuration
struct SequentialConfig {
    int maze_size;
    int num_particles;

    // Operator to uniquely compare configurations
    bool operator<(const SequentialConfig& other) const {
        return tie(maze_size, num_particles) < tie(other.maze_size, other.num_particles);
    }
};

// Structure for parallel configuration (maze size + particles + threads)
struct ParallelConfig {
    int maze_size;
    int num_particles;
    int num_threads;

    // Operator to uniquely compare configurations
    bool operator<(const ParallelConfig& other) const {
        return tie(maze_size, num_particles, num_threads) < tie(other.maze_size, other.num_particles, other.num_threads);
    }
};

class StatsUtil {
private:
    map<SequentialConfig, vector<double>> sequential_times;
    map<ParallelConfig, vector<double>> parallel_times;

public:
    // Add an execution time for the sequential algorithm
    void addSequentialTime(int maze_size, int num_particles, double time) {
        SequentialConfig config = {maze_size, num_particles};
        sequential_times[config].push_back(time);
    }

    // Add an execution time for the parallel algorithm
    void addParallelTime(int maze_size, int num_particles, int numThreads, double time) {
        ParallelConfig config = {maze_size, num_particles, numThreads};
        parallel_times[config].push_back(time);
    }

    // Calculate the mean of a vector of times
    double calculateMean(const vector<double>& times) {
        double total = 0.0;
        for (double time : times) {
            total += time;
        }
        return total / times.size();
    }

    // Calculate the standard deviation of a vector of times
    double calculateStdDev(const vector<double>& times) {
        double mean = calculateMean(times);
        double sum = 0.0;
        for (double time : times) {
            sum += (time - mean) * (time - mean);
        }
        return sqrt(sum / times.size());
    }

    // Calculate the average speedup for a specific configuration (requires both sequential and parallel times)
    double calculateSpeedup(const SequentialConfig& seqConfig, const ParallelConfig& parConfig) {
        double meanSeq = calculateMean(sequential_times[seqConfig]);
        double meanPar = calculateMean(parallel_times[parConfig]);
        return meanSeq / meanPar;
    }

    // Calculate the average efficiency for a specific parallel configuration
    double calculateEfficiency(const SequentialConfig& seqConfig, const ParallelConfig& parConfig) {
        return calculateSpeedup(seqConfig, parConfig) / parConfig.num_threads;
    }

    // Print the statistics for a sequential configuration
    void printStatisticsForSequentialConfig(const SequentialConfig& config) {
        cout << "\nStatistics for maze size " << config.maze_size << "x" << config.maze_size
                  << " with " << config.num_particles << " particles (Sequential):\n";
        cout << "Number of sequential runs: " << sequential_times[config].size() << endl;
        cout << "Mean sequential time: " << calculateMean(sequential_times[config]) << " ms" << endl;
        cout << "Sequential time standard deviation: " << calculateStdDev(sequential_times[config]) << " ms" << endl;
    }

    // Print the statistics for a parallel configuration
    void printStatisticsForParallelConfig(const ParallelConfig& config) {
        cout << "\nStatistics for maze size " << config.maze_size << "x" << config.maze_size
                  << " with " << config.num_particles << " particles and " << config.num_threads << " threads (Parallel):\n";
        cout << "Number of parallel runs: " << parallel_times[config].size() << endl;
        cout << "Mean parallel time: " << calculateMean(parallel_times[config]) << " ms" << endl;
        cout << "Parallel time standard deviation: " << calculateStdDev(parallel_times[config]) << " ms" << endl;
    }

    // Print the speedup and efficiency for a given configuration (matching both sequential and parallel)
    void printSpeedupAndEfficiency(const SequentialConfig& seqConfig, const ParallelConfig& parConfig) {
        cout << "Average speedup: " << calculateSpeedup(seqConfig, parConfig) << endl;
        cout << "Average efficiency: " << calculateEfficiency(seqConfig, parConfig) << endl;
    }

    // Print all statistics
    void printAllStatistics() {
        for (const auto& seqEntry : sequential_times) {
            const SequentialConfig& seqConfig = seqEntry.first;
            printStatisticsForSequentialConfig(seqConfig);

            // Print parallel stats that match this sequential config (for all threads)
            for (const auto& parEntry : parallel_times) {
                if (parEntry.first.maze_size == seqConfig.maze_size && parEntry.first.num_particles == seqConfig.num_particles) {
                    const ParallelConfig& parConfig = parEntry.first;
                    printStatisticsForParallelConfig(parConfig);
                    printSpeedupAndEfficiency(seqConfig, parConfig);
                }
            }
        }
    }

    void printSummarySequentialStats(const SequentialConfig sequential_config)
    {
        double time_ms = calculateMean(sequential_times[sequential_config]);
        cout << "Sequential time: " << time_ms << " ms("<< time_ms / 1000 << " s), "
                   << "Number of particles: " << sequential_config.num_particles << ", "
                   << "Maze dimensions: " << sequential_config.maze_size << " x " << sequential_config.maze_size
                   << endl << endl;
    }

    void printSummaryParallelStats(const ParallelConfig parallel_config)
    {
        double time_ms = calculateMean(parallel_times[parallel_config]);
        cout << "Parallel time: " << time_ms << " ms(" << time_ms / 1000 << " s), "
                   << "Number of thread: " << parallel_config.num_threads << ", "
                   << "Number of particles: " << parallel_config.num_particles << ", "
                   << "Maze dimensions: " << parallel_config.maze_size << " x " << parallel_config.maze_size
                   << endl << endl;
    }

    // Method to export all statistics to a CSV file
    void exportStatisticsToCSV(const string& filename)
    {
        ofstream file;
        file.open(filename);

        // Write header
        file << "MazeSize,NumParticles,NumThreads,SeqMeanTime,SeqStdDev,ParMeanTime,ParStdDev,Speedup,Efficiency\n";

        // Loop through all sequential configurations
        for (const auto& seqEntry : sequential_times)
        {
            const SequentialConfig& seqConfig = seqEntry.first;
            double seqMeanTime = calculateMean(sequential_times[seqConfig]);
            double seqStdDev = calculateStdDev(sequential_times[seqConfig]);

            // For each sequential configuration, find corresponding parallel configurations (same maze size and particles)
            for (const auto& parEntry : parallel_times)
            {
                if (parEntry.first.maze_size == seqConfig.maze_size && parEntry.first.num_particles == seqConfig.
                    num_particles)
                {
                    const ParallelConfig& parConfig = parEntry.first;
                    double parMeanTime = calculateMean(parallel_times[parConfig]);
                    double parStdDev = calculateStdDev(parallel_times[parConfig]);
                    double speedup = calculateSpeedup(seqConfig, parConfig);
                    double efficiency = calculateEfficiency(seqConfig, parConfig);

                    // Write data to CSV
                    file << seqConfig.maze_size << ","
                        << seqConfig.num_particles << ","
                        << parConfig.num_threads << ","
                        << seqMeanTime << ","
                        << seqStdDev << ","
                        << parMeanTime << ","
                        << parStdDev << ","
                        << speedup << ","
                        << efficiency << "\n";
                }
            }
        }

        file.close();
        cout << "Statistics exported to " << filename << "\n";
    }

    void reset() {
        sequential_times.clear();
        parallel_times.clear();
    }
};
