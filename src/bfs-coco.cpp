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
  uint32_t edge_ct;
  uint32_t edge_sz;
  uint32_t *neighbors;
} vertex_t;

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

  std::uniform_int_distribution<uint64_t> hash_b_dist(0, MERSENNE_61 - 1);
  std::uniform_int_distribution<uint64_t> hash_a_dist(1, MERSENNE_61 - 1);
  uint64_t hash_a = hash_a_dist(prng);
  uint64_t hash_b = hash_b_dist(prng);

  // Read input and distribute to appropriate nodes.
  MPI_File mpi_fp_in = MPI_FILE_NULL;
  MPI_Offset total_buf_sz = 0;
  MPI_File_open(MPI_COMM_WORLD, fname_in, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpi_fp_in);
  MPI_File_get_size(mpi_fp_in, &total_buf_sz);
  // Total number of edges = file size (bytes) / 4 (bytes/node) / 2 (nodes/edge)
  size_t total_edges = ((size_t)total_buf_sz) >> 3;
  size_t edges_per_machine = total_edges / machines;
  size_t edges_leftover = total_edges % machines;
  size_t edges = edges_per_machine + (rank < edges_leftover ? 1 : 0);
  size_t edges_buf_sz = edges << 8;
  size_t edges_buf_off = 0;

  uint32_t *edges_buf = (uint32_t *)malloc(edges_buf_sz);
  memset(edges_buf, 0, edges_buf_sz);
  int *send_sizes = (int *)malloc(machines * sizeof(int));
  int *send_caps = (int *)malloc(machines * sizeof(int));
  uint32_t **send_bufs = (uint32_t **)malloc(machines * sizeof(uint32_t **));
  memset(edges_buf, 0, machines * sizeof(uint32_t **));
  size_t mpi_fp_in_off = 0;
  for (int machine = 0; machine < machines; machine++) {
    send_bufs[machine] = (uint32_t *)malloc(2 * sizeof(uint32_t));
    send_caps[machine] = 2;
    send_sizes[machine] = 0;
    // Increment the file read start offset for this rank.
    if (machine < rank) {
      mpi_fp_in_off += (edges_per_machine + (machine < edges_leftover ? 1 : 0)) << 8;
    }
  }

  // Fill edge input buffers.
  MPI_File_set_view(mpi_fp_in, mpi_fp_in_off, MPI_UNSIGNED, MPI_UNSIGNED, "native", MPI_INFO_NULL);
  MPI_File_read(mpi_fp_in, edges_buf, edges_buf_sz, MPI_UNSIGNED, NULL);
  MPI_File_close(&mpi_fp_in);

  // Fill communication buffers.
  for (int i = 0; i < edges; i++) {
    // Get u and v from the file buffer.
    uint32_t u = edges_buf[edges_buf_off++];
    uint32_t v = edges_buf[edges_buf_off++];
    // Determine the machines housing u and v.
    uint64_t machine_u = ((hash_a * u + hash_b) % MERSENNE_61) % machines;
    uint64_t machine_v = ((hash_a * v + hash_b) % MERSENNE_61) % machines;
    // Resize machine[u]'s buffer if we're out of room.
    if (send_sizes[machine_u] == send_caps[machine_u]) {
      send_caps[machine_u] *= 2;
      send_bufs[machine_u] = (uint32_t *)realloc(send_bufs[machine_u],
						 send_caps[machine_u] * sizeof(uint32_t));
    }
    // Add (u,v) to buffer to send to machine[u].
    // By convention, the node housed in the machine is first.
    send_bufs[machine_u][send_sizes[machine_u]++] = u;
    send_bufs[machine_u][send_sizes[machine_u]++] = v;
    // Resize machine[v]'s buffer if we're out of room.
    if (send_sizes[machine_v] == send_caps[machine_v]) {
      send_caps[machine_v] *= 2;
      send_bufs[machine_v] = (uint32_t *)realloc(send_bufs[machine_v],
						 send_caps[machine_v] * sizeof(uint32_t));
    }
    // Add (u,v) to buffer to send to machine[v].
    // By convention, the node housed in the machine is first.
    send_bufs[machine_v][send_sizes[machine_v]++] = v;
    send_bufs[machine_v][send_sizes[machine_v]++] = u;
  }

  free(edges_buf);

  int *recv_counts = (int *)malloc(machines * sizeof(int));
  int *recv_displs = (int *)malloc(machines * sizeof(int));
  int total_recv_counts;
  memset(recv_counts, 0, sizeof(int) * machines);
  memset(recv_displs, 0, sizeof(int) * machines);
  // Disseminate edges belonging to each machine.
  for (int machine = 0; machine < machines; machine++) {
    // Gather edge counts for receiving machine
    // read from disk by sending machines.
    MPI_Gather(&send_sizes[machine], 1, MPI_INT,
	       recv_counts, 1, MPI_INT,
	       machine, MPI_COMM_WORLD);
    MPI_Reduce(&send_sizes[machine], &total_recv_counts,
	       1, MPI_INT, MPI_SUM,
	       machine, MPI_COMM_WORLD);
    if (rank == machine) {
      for (int i = 1; i < machines; i++) {
	recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
      }
      edges_buf = (uint32_t *)malloc(total_recv_counts * sizeof(uint32_t));
      memset(edges_buf, 0, total_recv_counts * sizeof(uint32_t));
    }
    MPI_Gatherv(send_bufs[machine], send_sizes[machine], MPI_UNSIGNED,
		edges_buf, recv_counts, recv_displs, MPI_UNSIGNED,
		machine, MPI_COMM_WORLD);
  }

  if (send_sizes) free( send_sizes );
  if (send_caps) free( send_caps );
  for (int machine = 0; machine < machines; machine++) {
    if (send_bufs[machine]) free( send_bufs[machine] );
  }
  if (send_bufs) free( send_bufs );
  if (recv_counts) free( recv_counts );
  if (recv_displs) free( recv_displs );

  // Construct local vertex-centric model.

  std::map<uint32_t, vertex_t *> V_machine;
  for (int i = 0; i < total_recv_counts; i += 2) {
    uint32_t u = edges_buf[i];
    uint32_t v = edges_buf[i+1];
    if (V_machine.find(u) == V_machine.end()) {
      vertex_t *vertex_u = (vertex_t *)malloc(sizeof(vertex_t));
      memset(vertex_u, 0, sizeof(vertex_t));
      vertex_u->id = u;
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
      V_machine[u]->neighbors[V_machine[u]->edge_ct++] = v;
    }
  }

  free( edges_buf );

  for (int machine = 0; machine < machines; machine++) {
    if (rank == machine) {
      for (auto &kv : V_machine) {
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

  // Tear down MPI.
  MPI_Finalize();

  // Exit.
  return 0;
}
