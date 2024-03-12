#include "common.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <mpi.h>
#include <thread>
#include <unordered_set>
#include <vector>
using namespace std;

// Create a vector of particles
std::vector<particle_t> particle_per_bin;
std::vector<particle_t> upper_halo;
std::vector<particle_t> lower_halo;

int num_intervals;

// Put any static global variables here that you will use throughout the simulation.
#define cutoff_squared (cutoff * cutoff)
#define minr_squared   (min_r * min_r)
#define INTERVAL                                                                                   \
    (1.05 * cutoff) // For correctness, this can't be less than cutoff. We add some leeway

double interval_;
double y_min;
double y_max;

int* send_counts = NULL;
int* send_displacements = NULL;

int* recv_counts = NULL;
int* recv_displacements = NULL;

int* part_sizes = NULL;
int* part_displacements = NULL;

std::vector<particle_t>* mini_bins = NULL;
int num_mini_bins;
int num_mini_intervals_y;
int num_mini_intervals_x;
int global_step = 0;
MPI_Comm mpi_active_world;

#define get_mini_row(y) static_cast<int>(std::floor(((y)-y_min + INTERVAL) / INTERVAL))
#define get_mini_col(x) static_cast<int>(std::floor((x) / INTERVAL));


#define BIN(r, c) mini_bins[(r) * num_mini_intervals_x + (c)]

void init_simulation(particle_t* parts, int num_parts, double size, int rank, int num_procs) {

    // Divide space into mini-bins:
    double mini_bin_size = INTERVAL;
    double num_bins_per_edge = std::ceil(static_cast<double>(size) / static_cast<double>(mini_bin_size));
    // Minimum bins per proc is 3 (or whatever -- arbitrarily chosen -- any value >= 1
    // should work. 3 just means we at least have bins dedicated to ghost regions WITHIN
    // the processor).
    //
    // Note: we can use the hardcoded 3 value to determine the minimum number of work
    // done by a processor, so we adjust the ratio of work per processor vs MPI
    // overhead. If we ahve time we can experiment with different values:
    const int MIN_BINS_PER_PROCESSOR = 3;
    int bins_per_proc = std::max(MIN_BINS_PER_PROCESSOR, static_cast<int>(std::ceil(num_bins_per_edge / num_procs)));
    // Number of horizontal sections
    num_intervals = std::ceil(static_cast<double>(num_bins_per_edge) / static_cast<double>(bins_per_proc));
    interval_ = size / num_intervals;

    // num_mini_intervals_y = current rank's bins + upper and lower halo from
    // neighboring processors. Note, each processor tracks a grid of mini-bins with
    // bin_row/bin_col relative to itself, it is not a global grid.
    num_mini_intervals_x = num_bins_per_edge;
    num_mini_intervals_y = bins_per_proc + 2;
    num_mini_bins = num_mini_intervals_x * num_mini_intervals_y;

    mini_bins = new std::vector<particle_t>[num_mini_bins];

    y_min = rank * interval_;
    y_max = (rank + 1) * interval_;

    send_counts = new int[num_intervals];
    send_displacements = new int[num_intervals];
    recv_counts = new int[num_intervals];
    recv_displacements = new int[num_intervals];

    // Create an MPI communicator consisting only of the processes that we will use.
    // This handles the case when num_procs > num_intervals, in which case we use a
    // world comm that only contains the processors/ranks that we are using
    int color = 0;
    if (rank >= num_intervals) {
        color = 1;
    }
    MPI_Comm_split(
        MPI_COMM_WORLD,
        color,
        rank,
        &mpi_active_world
    );


    if (rank == 0) {
        part_sizes = new int[num_intervals];
        part_displacements = new int[num_intervals];
    }

    //[DEBUG] Print out some basic info
    if (rank == 0) {
        std::cerr << "rank: " << rank << ", size: " << size << std::endl;
        std::cerr << "rank: " << rank << ", cutoff: " << cutoff << std::endl;
        std::cerr << "rank: " << rank << ", mini_bin_size: " << mini_bin_size << std::endl;
        std::cerr << "rank: " << rank << ", num_bins_per_edge: " << num_bins_per_edge << std::endl;
        std::cerr << "rank: " << rank << ", bins_per_proc: " << bins_per_proc << std::endl;
        std::cerr << "rank: " << rank << ", num_intervals: " << num_intervals << std::endl;
        std::cerr << "rank: " << rank << ", INTERVAL: " << INTERVAL << std::endl;
        std::cerr << "rank: " << rank << ", interval_: " << interval_ << std::endl;
        std::cerr << "rank: " << rank << ", num_procs: " << num_procs << std::endl;
        std::cerr << "rank: " << rank << ", color: " << color << std::endl;
        std::cerr << "rank: " << rank << ", num_mini_intervals_x: " << num_mini_intervals_x << std::endl;
        std::cerr << "rank: " << rank << ", num_mini_intervals_y: " << num_mini_intervals_y << std::endl;
        std::cerr << "rank: " << rank << ", num_mini_bins: " << num_mini_bins << std::endl;
        std::cerr << "rank: " << rank << ", y_min: " << y_min << ", y_max: " << y_max << std::endl;
    }

    // Every processor should go through and figure out what particles are in its domain
    //  then add these particels to the vector "particles_per_bin"
    for (int particle = 0; particle < num_parts; particle++) {
        int r = static_cast<int>(std::floor(parts[particle].y / interval_)); // row index
        parts[particle].ax = 0;
        parts[particle].ay = 0;
        if (r == rank) {
            particle_per_bin.push_back(
                parts[particle]); // push particle to bin's vector in the array
        }
    }
    global_step = 0;
}

double max_ax_seen = 0.0;
double max_ay_seen = 0.0;
// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff_squared)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;

    if (fabs(coef * dx) > max_ax_seen) {
        max_ax_seen = fabs(coef * dx);
    }

    if (fabs(coef * dy) > max_ay_seen) {
        max_ay_seen = fabs(coef * dy);
    }
}

int move(particle_t& p, double size, int rank) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
    int row = static_cast<int>(std::floor(p.y / interval_));

    if (p.id == 4) {
        LOG("[%04d,%03d] 2 move(): id:%d, (ax:%f, ay%f), (vx:%f, vy%f), (x:%f, y%f) row: %d", global_step, rank, p.id, p.ax, p.ay, p.vx, p.vy, p.x, p.y, row);
    }

    if (row != rank) {
        // LOG("A particle is leaving bin %d to %d, p.y = %lf interval = %lf, p.y < 0? %d", rank,
        // row, p.y, interval_, p.y < 0);
        return row; // Signal to move to different bin
    }
    else {
        return -1; // No movement == -1
    }
}

void apply_force_per_bin(int p0, std::vector<particle_t> bin, int bin_row, int bin_col) {
    if (particle_per_bin[p0].id == 4) {
        LOG("[%04d,???] 1b apply_force_per_bin(bin:[%d, %d]): id:%d, (ax:%f, ay:%f), (vx:%f, vy%f), (x:%f, y%f)", global_step, bin_row, bin_col, particle_per_bin[p0].id, particle_per_bin[p0].ax, particle_per_bin[p0].ay, particle_per_bin[p0].vx, particle_per_bin[p0].vy, particle_per_bin[p0].x, particle_per_bin[p0].y);
    }
    for (particle_t& p1 : bin) {
        if (p1.id == particle_per_bin[p0].id)
            continue; // Don't apply my force on myself

        // std::cerr << p0 << " " << p1 << std::endl;

        if (particle_per_bin[p0].id == 4) {
            LOG("[%04d,???] 1c applying_f with other particle: id:%d, (:%f, %f)", global_step, p1.id, p1.x, p1.y)
        }
        apply_force(particle_per_bin[p0], p1);
    }
}

// ##############################
//       Apply Force
// ##############################

// For each particle; calculate force with
//  (a) other particles in that bin
//  (b) particles in "upstairs" vector
//  (c) particles in "downstairs" vector
void apply_force_mini_grid(int rank, std::vector<particle_t>& upstairs,
    std::vector<particle_t>& downstairs) {
    for (int i = 0; i < num_mini_intervals_y; i++) {
        for (int j = 0; j < num_mini_intervals_x; j++) {
            BIN(i, j).clear();
        }
    }
    for (auto& p : particle_per_bin) { // looping over particles p in 
        int r = get_mini_row(p.y);    // row index
        int c = get_mini_col(p.x);    // column index
        BIN(r, c).push_back(p);      
    }
    for (auto& p : upstairs) {
        int r = get_mini_row(p.y);
        int c = get_mini_col(p.x);
        BIN(r, c).push_back(p);
    }
    for (auto& p : downstairs) {
        int r = get_mini_row(p.y);
        int c = get_mini_col(p.x);
        BIN(r, c).push_back(p);
    }

    for (int particle = 0; particle < particle_per_bin.size(); particle++) {

        particle_per_bin[particle].ax = particle_per_bin[particle].ay = 0.0;
        int r = get_mini_row(particle_per_bin[particle].y);
        int c = get_mini_col(particle_per_bin[particle].x);


        for (int dr = -1; dr <= 1; dr++) {
            for (int dc = -1; dc <= 1; dc++) {
                if (r + dr < 0 || r + dr >= num_mini_intervals_y // Shouldn't go out of bounds
                    || c + dc < 0 || c + dc >= num_mini_intervals_x)
                    continue;


                apply_force_per_bin(particle, BIN(r + dr, c + dc), r + dr, c + dc);
            }
        }
    }
}


void print_part(int rank, int print_num) {
    if (print_num != 6 and print_num != 0) {
        return;
    }
    for (int i = 0; i < particle_per_bin.size(); i++) {
        auto p = particle_per_bin[i];
        if (p.id == 4 or p.id == 3 or p.id == 6 or p.id == 7) {
            LOG("[%04d,%03d] %d id:%d, (ax:%f, ay:%f), (vx:%f, vy:%f), (x:%f, y:%f)", global_step, rank, print_num, p.id, p.ax, p.ay, p.vx, p.vy, p.x, p.y);
        }
    }
}


// up/down here is relative to the receiver
#define TAG_SIZE_UP   12
#define TAG_SIZE_DOWN 1
#define TAG_DATA_UP   11
#define TAG_DATA_DOWN 0
void simulate_one_step(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    if (rank >= num_intervals) {
        return;
    }

    upper_halo.clear();
    lower_halo.clear();
    for (int i = 0; i < particle_per_bin.size(); i++) {
        if (particle_per_bin[i].y - y_min < INTERVAL) {
            lower_halo.push_back(particle_per_bin[i]);
        }

        if (y_max - particle_per_bin[i].y < INTERVAL) {
            upper_halo.push_back(particle_per_bin[i]);
        }
    }

    int upper_halo_size = upper_halo.size();
    int lower_halo_size = lower_halo.size();
    int upstairs_size = 0;
    int downstairs_size = 0;

    // Send upper_halo_size to upstairs, receive upstairs_size from upstairs:
    if (rank < (num_intervals - 1)) {
        MPI_Sendrecv(&upper_halo_size,      // Send buffer
                    1,                      // Number of elements in send buffer
                    MPI_INT,                // Type of elements in send buffer
                    rank + 1,               // Destination rank
                    TAG_SIZE_DOWN,          // Tag
                    &upstairs_size,         // Receive buffer
                    1,                      // Number of elements in receive buffer
                    MPI_INT,                // Type of elements in receive buffer
                    rank + 1,               // Source rank
                    TAG_SIZE_UP,            // Tag
                    mpi_active_world,       // Communicator
                    MPI_STATUS_IGNORE);     // Status
    }

    // Send lower_halo_size to downstairs, receive downstairs_size from downstairs:
    if ((rank >= 1) && (rank < num_intervals)) {
        MPI_Sendrecv(&lower_halo_size,      // Send buffer
                       1,                   // Number of elements in send buffer
                       MPI_INT,             // Type of elements in send buffer
                       rank - 1,            // Destination rank
                       TAG_SIZE_UP,         // Tag
                       &downstairs_size,    // Receive buffer
                       1,                   // Number of elements in receive buffer
                       MPI_INT,             // Type of elements in receive buffer
                       rank - 1,            // Source rank
                       TAG_SIZE_DOWN,       // Tag
                       mpi_active_world,    // Communicator
                       MPI_STATUS_IGNORE);  // Status
    }

    std::vector<particle_t> upstairs(upstairs_size);
    std::vector<particle_t> downstairs(downstairs_size);

    int np2receive;    // Used by MPI Probe to check size of incoming message
    MPI_Status status; // Used by MPI to check status of incoming message

    // Bins[0:N-1]
    if (rank < (num_intervals - 1)) {
        int mpi_state = MPI_Sendrecv(upper_halo.data(), upper_halo_size, PARTICLE, rank + 1,
            TAG_DATA_DOWN, upstairs.data(), upstairs_size, PARTICLE,
            rank + 1, TAG_DATA_UP, mpi_active_world, &status);

        MPI_Get_count(&status, PARTICLE, &np2receive); // get size of incoming message
        upstairs.resize(np2receive);                   // resize the vector to the size of the incoming message
    }
    // Bins[1:N]
    if ((rank >= 1) && (rank < num_intervals)) {
        int mpi_state = MPI_Sendrecv(lower_halo.data(), lower_halo_size, PARTICLE, rank - 1,
            TAG_DATA_UP, downstairs.data(), downstairs_size, PARTICLE,
            rank - 1, TAG_DATA_DOWN, mpi_active_world, &status);
        if (mpi_state != MPI_SUCCESS) {
        }
        MPI_Get_count(&status, PARTICLE, &np2receive);  // get size of incoming message
        downstairs.resize(np2receive);                  // resize the vector to the size of the incoming message
    }

    // Probably best to make sure all messages have been sent & received
    // Send and Recv are blocking. No Barriers needed.
    apply_force_mini_grid(rank, upstairs, downstairs);

    // Move particles
    std::map<int, std::vector<particle_t>> scatter_parts; 
    std::unordered_set<uint64_t> remove_idx;
    for (int i = 0; i < particle_per_bin.size(); ++i) {
        int send2bin = move(particle_per_bin[i], size, rank);
        if (send2bin != -1) {
            scatter_parts[send2bin].push_back(particle_per_bin[i]);
            remove_idx.insert(particle_per_bin[i].id);
        }
    }

    particle_per_bin.erase(std::remove_if(particle_per_bin.begin(), particle_per_bin.end(),
        [&remove_idx](const particle_t& p) {
            return remove_idx.find(p.id) != remove_idx.end();
        }),
        particle_per_bin.end());

    std::vector<particle_t> send_buff;
    int curr_send_displacement = 0;
    for (int i = 0; i < num_intervals; i++) {
        auto it = scatter_parts.find(i);
        if (it == scatter_parts.end()) {
            send_counts[i] = 0;
            send_displacements[i] = curr_send_displacement;
            continue;
        }

        for (auto& x : it->second) {
            send_buff.push_back(x);
        }

        send_counts[i] = it->second.size();
        send_displacements[i] = curr_send_displacement;
        curr_send_displacement += it->second.size();
    }

    // First we tell each proc j how much it must be expecting from proc i
    // This is an all_to_all
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, mpi_active_world);

    int curr_recv_displacement = 0;
    for (int i = 0; i < num_intervals; i++) {
        recv_displacements[i] = curr_recv_displacement;
        curr_recv_displacement += recv_counts[i];
    }

    // curr_recv_displacements now is the total size of particles I am going to receive.
    std::vector<particle_t> recv_buff(curr_recv_displacement);

    // Now, let's exchange particles.
    MPI_Alltoallv(send_buff.data(), send_counts, send_displacements, PARTICLE, recv_buff.data(),
        recv_counts, recv_displacements, PARTICLE, mpi_active_world);

    for (auto& p : recv_buff) {
        particle_per_bin.push_back(p);
    }
}

bool compare_particles(const particle_t& a, const particle_t& b) { return a.id < b.id; }

void gather_for_save(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    // Write this function such that at the end of it, the master (rank == 0)
    // processor has an in-order view of all particles. That is, the array
    // parts is complete and sorted by particle id.
    if (rank >= num_intervals) {
        return;
    }

    int my_parts = particle_per_bin.size();

    MPI_Gather(&my_parts, 1, MPI_INT, part_sizes, 1, MPI_INT, 0, mpi_active_world);

    int curr_displacement = 0;
    if (rank == 0) {
        for (int i = 0; i < num_intervals; i++) {
            part_displacements[i] = curr_displacement;
            curr_displacement += part_sizes[i];
        }
    }

    MPI_Gatherv(particle_per_bin.data(), my_parts, PARTICLE, parts, part_sizes, part_displacements,
        PARTICLE, 0, mpi_active_world);

    if (rank == 0) {
        std::sort(parts, parts + num_parts, compare_particles);
    }


}

// Syntax of MPI

// MPI_Probe(
// int source,
// int tag,
// MPI_Comm comm,
// MPI_Status* status)

// MPI_Send(
//     void* data,
//     int count,
//     MPI_Datatype datatype,
//     int destination,
//     int tag,
//     MPI_Comm communicator)

// MPI_Recv(
//     void* data,
//     int count,
//     MPI_Datatype datatype,
//     int source,
//     int tag,
//     MPI_Comm communicator,
//     MPI_Status* status)