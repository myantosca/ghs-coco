#include <cstdio>
#include <cstdint>
#include <cstring>
#include <sys/time.h>
#include <random>
#include <map>
#include <vector>
#include <iostream>
#include <mpi.h>

// https://en.wikipedia.org/wiki/Mersenne_prime#List_of_known_Mersenne_primes
#define MERSENNE_61 2305843009213693951

typedef struct vertex
{
  uint32_t id;
  uint32_t parent;
  uint32_t group;
  uint32_t group_ct;
  uint32_t awaiting;
  uint32_t edge_ct;
  uint32_t edge_sz;
  uint32_t *neighbors;
} vertex_t;

void exchange(int rank, int machines,
	      int *send_counts, MPI_Datatype send_type, uint32_t **send_bufs,
	      int *recv_counts, int *recv_displs, int *recv_totals,
	      MPI_Datatype recv_type, uint32_t **recv_buf) {
  memset(recv_counts, 0, sizeof(int) * machines);
  memset(recv_displs, 0, sizeof(int) * machines);
  // Gather for each machine in turn from every other machine's
  // targeted send buffers.
  for (int machine = 0; machine < machines; machine++) {
    // Determine the receive sub-buffer element counts.
    MPI_Gather(&send_counts[machine], 1, MPI_INT,
	       recv_counts, 1, MPI_INT,
	       machine, MPI_COMM_WORLD);
    // Only one machine is target of the gather operation at a time.
    if (rank == machine) {
      for (int i = 1; i < machines; i++) {
	recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
      }
      // Calculate the total receive buffer size so each machine knows its extent.
      recv_totals[machine] = recv_counts[machines - 1] + recv_displs[machines - 1];
      // Allocate the receive buffer according to the calculated size.
      *recv_buf = (uint32_t *)malloc(recv_totals[machine] * sizeof(uint32_t));
      // Zero out the receive buffer as a matter of safety.
      memset(*recv_buf, 0, recv_totals[machine] * sizeof(uint32_t));
    }
    // Gather all the targeted send buffers on the target machine
    // in a contiguous array in order by sender rank.
    MPI_Gatherv(send_bufs[machine], send_counts[machine], MPI_UNSIGNED,
		*recv_buf, recv_counts, recv_displs, MPI_UNSIGNED,
		machine, MPI_COMM_WORLD);
  }
}

int main(int argc, char *argv[]) {
  int rank, machines;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &machines);
  // Default input filename.
  const char *fname_in = "./graph.ecg";

  int a = 0;
  while (a < argc) {
    // Output filename argument.
    if (!strcmp("-o", argv[a])) {
      a++;
      if (a < argc) fname_in = argv[a];
    }
    a++;
  }

  // Determine hash function parameters.
  std::random_device rd;
  std::mt19937 prng(rd());

  uint64_t hash_a = 1;
  uint64_t hash_b = 0;
  if (rank == 0) {
    std::uniform_int_distribution<uint64_t> hash_b_dist(0, MERSENNE_61 - 1);
    std::uniform_int_distribution<uint64_t> hash_a_dist(1, MERSENNE_61 - 1);
    hash_b = hash_a_dist(prng);
    hash_a = hash_b_dist(prng);
  }

  MPI_Bcast(&hash_a, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hash_b, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  // Non-uniform gather info.
  int *send_counts = (int *)malloc(machines * sizeof(int));
  int *recv_counts = (int *)malloc(machines * sizeof(int));
  int *recv_displs = (int *)malloc(machines * sizeof(int));
  int *recv_totals = (int *)malloc(machines * sizeof(int));
  uint32_t *recv_buf = NULL;

  // Read input and distribute to appropriate nodes.
  MPI_File mpi_fp_in = MPI_FILE_NULL;
  MPI_Offset total_buf_sz = 0;
  MPI_File_open(MPI_COMM_WORLD, fname_in, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpi_fp_in);
  MPI_File_get_size(mpi_fp_in, &total_buf_sz);
  // Total number of edges = file size (bytes) / 4 (bytes/node) / 2 (nodes/edge)
  size_t total_edges = ((size_t)total_buf_sz) >> 3;
  size_t edges_per_machine = total_edges / machines;
  size_t edges_leftover = total_edges % machines;
  // Evenly distribute remainder of edges if not evenly divisible over
  // the first machines until the remainder is exhausted.
  size_t edges = edges_per_machine + (rank < edges_leftover ? 1 : 0);
  // Total edge buffer size = number of edges * 4 (bytes/node) * 2 (nodes/edge)
  size_t edges_buf_sz = edges << 3;

  uint32_t *edges_buf = (uint32_t *)malloc(edges_buf_sz);
  memset(edges_buf, 0, edges_buf_sz);
  int *send_caps = (int *)malloc(machines * sizeof(int));
  uint32_t **send_bufs = (uint32_t **)malloc(machines * sizeof(uint32_t **));
  memset(send_bufs, 0, machines * sizeof(uint32_t **));
  size_t mpi_fp_in_off = 0;

  // Initialize send buffers for edge exchange.
  for (int machine = 0; machine < machines; machine++) {
    send_bufs[machine] = (uint32_t *)malloc(2 * sizeof(uint32_t));
    send_caps[machine] = 2;
    send_counts[machine] = 0;
    // Increment the file read start offset for this rank.
    if (machine < rank) {
      mpi_fp_in_off += (edges_per_machine + (machine < edges_leftover ? 1 : 0)) << 3;
    }
  }

  // Fill edge input buffers.
  MPI_File_set_view(mpi_fp_in, mpi_fp_in_off, MPI_UNSIGNED, MPI_UNSIGNED, "native", MPI_INFO_NULL);
  MPI_File_read(mpi_fp_in, edges_buf, edges << 1, MPI_UNSIGNED, NULL);
  MPI_File_close(&mpi_fp_in);

  // Fill send buffers for edge exchange.
  for (int i = 0; i < edges; i++) {
    // Get u and v from the file buffer.
    uint32_t u = edges_buf[i*2];
    uint32_t v = edges_buf[i*2+1];
    // Determine the machines housing u and v.
    uint64_t machine_u = ((hash_a * u + hash_b) % MERSENNE_61) % machines;
    uint64_t machine_v = ((hash_a * v + hash_b) % MERSENNE_61) % machines;
    // Resize machine[u]'s buffer if we're out of room.
    if (send_counts[machine_u] == send_caps[machine_u]) {
      send_caps[machine_u] *= 2;
      send_bufs[machine_u] = (uint32_t *)realloc(send_bufs[machine_u],
						 send_caps[machine_u] * sizeof(uint32_t));
    }
    // Add (u,v) to buffer to send to machine[u].
    // By convention, the node housed in the machine is first.
    send_bufs[machine_u][send_counts[machine_u]++] = u;
    send_bufs[machine_u][send_counts[machine_u]++] = v;
    // Resize machine[v]'s buffer if we're out of room.
    if (send_counts[machine_v] == send_caps[machine_v]) {
      send_caps[machine_v] *= 2;
      send_bufs[machine_v] = (uint32_t *)realloc(send_bufs[machine_v],
						 send_caps[machine_v] * sizeof(uint32_t));
    }
    // Add (v,u) to buffer to send to machine[v].
    // By convention, the node housed in the machine is first.
    send_bufs[machine_v][send_counts[machine_v]++] = v;
    send_bufs[machine_v][send_counts[machine_v]++] = u;
  }

  // Clean up so as to not hoard memory.
  free(edges_buf);

  // Exchange edges.
  exchange(rank, machines,
	   send_counts, MPI_UNSIGNED, send_bufs,
	   recv_counts, recv_displs, recv_totals,
	   MPI_UNSIGNED, &recv_buf);

  // Clean up so as to not hoard memory.
  if (send_caps) free( send_caps );
  for (int machine = 0; machine < machines; machine++) {
    if (send_bufs[machine]) free( send_bufs[machine] );
  }

  // Construct local vertex-centric model.
  std::map<uint32_t, vertex_t *> V_machine;
  for (int i = 0; i < recv_totals[rank]; i += 2) {
    uint32_t u = recv_buf[i];
    uint32_t v = recv_buf[i+1];
    // If the node has not been created yet, do so.
    if (V_machine.find(u) == V_machine.end()) {
      vertex_t *vertex_u = (vertex_t *)malloc(sizeof(vertex_t));
      memset(vertex_u, 0, sizeof(vertex_t));
      vertex_u->id = u;
      vertex_u->parent = u;
      vertex_u->group = u;
      vertex_u->group_ct = 1; // population of BFS sub-tree including itself
      // Set state variable. First decrement is upon receipt of broadcast.
      vertex_u->awaiting = 1;
      vertex_u->edge_ct = 0;
      vertex_u->edge_sz = 1;
      vertex_u->neighbors = (uint32_t *)malloc(sizeof(uint32_t));
      V_machine[u] = vertex_u;
    }
    if (u != v) {
      if (V_machine[u]->edge_ct == V_machine[u]->edge_sz) {
	V_machine[u]->edge_sz *= 2;
	V_machine[u]->neighbors = (uint32_t *)realloc(V_machine[u]->neighbors,
						      V_machine[u]->edge_sz * sizeof(uint32_t));
      }
      // Increment edge count and add newest neighbor.
      V_machine[u]->neighbors[V_machine[u]->edge_ct++] = v;
      // Increment state variable for counting down during upcast phase.
      V_machine[u]->awaiting++;
    }
  }

  free( recv_buf );

  // Debug printout of graph in vertex-centric format.
  for (int machine = 0; machine < machines; machine++) {
    if (rank == machine) {
      for (auto &kv : V_machine) {
	std::cout << "[" << machine << "]";
	std::cout << kv.first << ":";
	for (int i = 0; i < kv.second->edge_ct; i++) {
	  std::cout << " " << kv.second->neighbors[i];
	}
	std::cout << std::endl;
      }
    }
  }

  // Execute BFS search for components.


  // Clean up.
  if (send_counts) free( send_counts );
  if (recv_counts) free( recv_counts );
  if (recv_displs) free( recv_displs );
  if (recv_totals) free( recv_totals );

  for (auto &kv : V_machine) {
    if (kv.second) free(kv.second);
    kv.second = NULL;
  }
  V_machine.clear();
  // Tear down MPI.
  MPI_Finalize();

  // Exit.
  return 0;
}
