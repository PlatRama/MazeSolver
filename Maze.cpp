//
// Created by Plator Rama on 03/09/2024.
//

#include <iostream>
#include <stack>
#include <vector>
#include <random>

using namespace std;

struct Point
{
    int x, y;
};

struct Particle
{
    Point currentPosition;
    vector<Point> path;
};

class Maze
{
    char wall_character = '#';
    char not_wall_character = '.';
    char start_character = 'S';
    char exit_character = 'E';
    char path_character = '*';

public:
    char** grid;
    int m;
    int M;
    int n;
    int N;
    Point start{};
    Point exit{};

    explicit Maze(const int _m = 4, const int _n = 4)
    {
        m = _m;
        n = _n;

        M = 2 * m + 1;
        N = 2 * n + 1;

        grid = new char*[M];

        start = {1, 1};
        exit = {2 * m, 2 * n - 1};

        initMaze();
        generateMaze();
    }

    ~Maze()
    {
        // Free dynamically allocated memory for each row of the maze
        for (int i = 0; i < M; ++i)
            delete[] grid[i];
        // Free the array of pointers
        delete[] grid;
    }

    // Copy constructor: creates a deep copy of the maze, including all members and dynamically allocated grid
    Maze(const Maze& other)
    {
        copyFrom(other);
    }

    // Assignment operator: performs a deep copy of the maze, handling self-assignment and memory deallocation
    Maze& operator=(const Maze& other)
    {
        if (this == &other)
            return *this;

        this->~Maze();

        // Copy data from the other object
        copyFrom(other);
        return *this;
    }


    // Utility function to perform the deep copy
    void copyFrom(const Maze& other)
    {
        m = other.m;
        M = other.M;
        n = other.n;
        N = other.N;
        start = other.start;
        exit = other.exit;

        grid = new char*[M];
        for (int i = 0; i < M; ++i)
        {
            grid[i] = new char[N];
            for (int j = 0; j < N; ++j)
                grid[i][j] = other.grid[i][j];
        }
    }

    void initMaze()
    {
        for (int i = 0; i < M; i++)
            grid[i] = new char[N];

        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
                //Assigns '#' if either i or j is even, otherwise assigns ' '.
                grid[i][j] = !(i & 1 && j & 1) ? wall_character : not_wall_character;
        }

        for (int i = 1; i < M; i += 2)
            for (int j = 1; j < N; j += 2)
                grid[i][j] = not_wall_character;
    }

    /**
     * @brief Generates a maze using the Depth-First Search (DFS) algorithm.
     *
     * This method generates a maze by starting from a random cell and performing
     * a Depth-First Search (DFS) to carve out a path through the grid. The grid
     * is represented as a 2D matrix, where walls and paths are defined.
     *
     * The algorithm continues to explore neighboring cells that haven't been visited,
     * and creates a path by removing walls between adjacent cells. If a cell has no
     * unvisited neighbors, the algorithm backtracks to the previous cell and continues
     * from there.
     *
     * The maze is generated with an entrance ('S') at the top and an exit ('E')
     * at the bottom. This method assumes that the grid dimensions (M, N) and the
     * maze matrix (`maze`) are defined and initialized beforehand.
     *
     * The method employs randomization to ensure that each generated maze is unique.
     */
    void generateMaze() const
    {
        // Vector to store the list of all potential cells in the grid
        vector<pair<int, pair<int, int>>> cell_list;
        cell_list.reserve(m * n); // Preallocate memory for efficiency

        // Vector to keep track of visited cells in the maze
        vector<bool> visited_cells(m * n, false);

        // Stack to manage the path during the Depth-First Search
        stack<pair<int, pair<int, int>>> cell_stack;

        // Random number generator setup
        random_device rand_dev;
        mt19937 rng(rand_dev());
        uniform_int_distribution<mt19937::result_type> cell_dist(0, m * n - 1); // Fit within cell_list size

        int num_visited_cells = 0; // Counter for the number of visited cells
        int cell_index = 0; // Index for each cell in cell_list

        // Populate the list of cells with their respective indices and grid positions
        for (int row = 1; row < M; row += 2)
        {
            for (int col = 1; col < N; col += 2)
            {
                cell_list.emplace_back(cell_index, make_pair(row, col));
                cell_index++;
            }
        }

        // Start the DFS from a random cell
        int random_index = cell_dist(rng);
        cell_stack.push(cell_list[random_index]);
        visited_cells[random_index] = true;
        num_visited_cells++;

        // DFS loop: Continue until all cells have been visited
        while (num_visited_cells < m * n)
        {
            // Get the current cell from the top of the stack
            auto current_cell = cell_stack.top();
            vector<int> neighbors; // To store available neighbors' directions

            // Check each direction (North, East, South, West) for valid neighbors
            if (current_cell.second.first > 1)
                if (grid[current_cell.second.first - 2][current_cell.second.second]
                    && !visited_cells[getIdx(current_cell.second.first - 2, current_cell.second.second, cell_list)])
                    neighbors.push_back(0); // North

            if (current_cell.second.second < N - 2)
                if (grid[current_cell.second.first][current_cell.second.second + 2]
                    && !visited_cells[getIdx(current_cell.second.first, current_cell.second.second + 2, cell_list)])
                    neighbors.push_back(1); // East

            if (current_cell.second.first < M - 2)
                if (grid[current_cell.second.first + 2][current_cell.second.second]
                    && !visited_cells[getIdx(current_cell.second.first + 2, current_cell.second.second, cell_list)])
                    neighbors.push_back(2); // South

            if (current_cell.second.second > 1)
                if (grid[current_cell.second.first][current_cell.second.second - 2]
                    && !visited_cells[getIdx(current_cell.second.first, current_cell.second.second - 2, cell_list)])
                    neighbors.push_back(3); // West

            // If there are valid neighboring cells, choose one randomly to visit next
            if (!neighbors.empty())
            {
                int next_idx;
                int next_cell_direction = neighbors[cell_dist(rng) % neighbors.size()]; // Randomly select a direction

                // Create a path to the chosen neighbor
                switch (next_cell_direction)
                {
                case 0: // North
                    grid[current_cell.second.first - 1][current_cell.second.second] = not_wall_character;
                    next_idx = getIdx(current_cell.second.first - 2, current_cell.second.second, cell_list);
                    break;
                case 1: // East
                    grid[current_cell.second.first][current_cell.second.second + 1] = not_wall_character;
                    next_idx = getIdx(current_cell.second.first, current_cell.second.second + 2, cell_list);
                    break;
                case 2: // South
                    grid[current_cell.second.first + 1][current_cell.second.second] = not_wall_character;
                    next_idx = getIdx(current_cell.second.first + 2, current_cell.second.second, cell_list);
                    break;
                case 3: // West
                    grid[current_cell.second.first][current_cell.second.second - 1] = not_wall_character;
                    next_idx = getIdx(current_cell.second.first, current_cell.second.second - 2, cell_list);
                    break;
                default: next_idx = -1;
                }

                // Push the new cell onto the stack and mark it as visited
                cell_stack.push(cell_list[next_idx]);
                visited_cells[next_idx] = true;
                num_visited_cells++;
            }
            else
            {
                // If no neighbors are available, backtrack by popping the stack
                cell_stack.pop();
            }
        }

        // Define the maze entrance and exit
        grid[0][1] = start_character; // Start
        grid[2 * m][2 * n - 1] = exit_character; // Exit
    }


    static int getIdx(const int x, const int y, const vector<pair<int, pair<int, int>>>& cell_list)
    {
        for (int i = 0; i < cell_list.size(); i++)
        {
            if (cell_list[i].second.first == x && cell_list[i].second.second == y)
                return cell_list[i].first;
        }
        cout << "getIdx() couldn't find the index!" << endl;
        return -1;
    }

    /**
     * @brief Retrieves all valid moves for a given particle's current position.
     *
     * This method checks the four possible directions (right, left, down, up) from the particle's
     * current position and returns a vector containing all valid moves.
     *
     * @param particle The particle whose valid moves are to be determined.
     * @return A vector of valid moves, each represented as a Point.
     */
    vector<Point> getValidMoves(const Particle& particle)
    {
        // Possible moves: (dx, dy) represents directions: right, left, down, up
        const vector<Point> possibleMoves = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
        vector<Point> validMoves;

        const int x = particle.currentPosition.x;
        const int y = particle.currentPosition.y;

        // Check each possible move and add to validMoves if it's a valid move
        for (const auto& move : possibleMoves)
        {
            if (isValidMove(x + move.x, y + move.y))
                validMoves.push_back(move);
        }
        return validMoves;
    }

    /**
     * @brief Moves the particle to a random valid position within the maze.
     *
     * This method retrieves all valid moves for the particle's current position,
     * selects one at random, and updates the particle's position accordingly.
     * The maze is updated to reflect the particle's movement, and the particle's
     * path is recorded.
     *
     * @param particle The particle to move within the maze.
     */
    void randomMove(Particle& particle)
    {
        // Retrieve all valid moves
        vector<Point> validMoves = getValidMoves(particle);

        // Randomly select one of the valid moves
        const Point& move = validMoves[rand() % validMoves.size()];

        // Calculate the new position
        int newX = particle.currentPosition.x + move.x;
        int newY = particle.currentPosition.y + move.y;

        particle.currentPosition = {newX, newY}; // Update particle's position

        // Record the new position in the particle's path
        particle.path.push_back(particle.currentPosition);
    }

    /**
     * @brief Checks if a move to a given position in the maze is valid.
     *
     * This method verifies that the position (x, y) is within the bounds of the maze,
     * is not a wall, and is not the start position.
     *
     * @param row_index The row index in the maze.
     * @param column_index The column index in the maze.
     * @return true if the move is valid; false otherwise.
     */
    bool isValidMove(int row_index, int column_index)
    {
        return row_index >= 0 && row_index < M && column_index >= 0 && column_index < N // within the bounds
            && grid[row_index][column_index] != wall_character // is not a wall
            && grid[row_index][column_index] != start_character; // is not the start position
    }


    bool isExitFound(const Point& pos)
    {
        return pos.x == exit.x && pos.y == exit.y;
    }


    void displayMaze()
    {
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
                cout << grid[i][j] << " ";
            cout << endl;
        }
        cout << endl;
    }

    void drawPath(Particle winning_particle)
    {
        for (auto& [x, y] : winning_particle.path)
        {
            grid[x][y] = path_character;
        }
        grid[2 * m][2 * n - 1] = exit_character;
    }
};
