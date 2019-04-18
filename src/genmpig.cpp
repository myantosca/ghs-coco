#include <cstdio>
#include <cstdint>
#include <cstring>
#include <random>
#include <mpi.h>

int main(int argc, char *argv[]) {
  int rank, machines, edges_per_block = 1;
  double p_edge = 0;
  uint32_t n = 1;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &machines);

  // Default output filename.
  const char *fname_out = "./graph.ecg";
  int a = 0;
  while (a < argc) {
    // Edge probability argument.
    if (!strcmp("-p", argv[a])) {
      a++;
      if ((a < argc) && ((sscanf(argv[a], "%lf", &p_edge) != 1) || (p_edge < 0) || (p_edge > 1))) {
	if (rank == 0) {
	  fprintf(stderr, "Invalid edge probability: %s\n", argv[a]);
	}
	return -1;
      }
    }
    // Population size argument.
    if (!strcmp("-n", argv[a])) {
      a++;
      if ((a < argc) && (sscanf(argv[a], "%u", &n) != 1)) {
	if (rank == 0) {
	  fprintf(stderr, "Invalid population size: %s\n", argv[a]);
	}
	return -1;
      }
    }
    // Output filename argument.
    if (!strcmp("-o", argv[a])) {
      a++;
      if (a < argc) fname_out = argv[a];
    }
    // Measure of edges per output file block per machine.
    if (!strcmp("-b", argv[a])) {
      a++;
      if ((a < argc) && (sscanf(argv[a], "%u", &edges_per_block) != 1)) {
	if (rank == 0) {
	  fprintf(stderr, "Invalid output block size: %s\n", argv[a]);
	}
	return -1;
      }
    }
    a++;
  }

  MPI_File mpi_fp_out = MPI_FILE_NULL;
  MPI_Request dump_req = MPI_REQUEST_NULL;
  MPI_Status dump_status;

  size_t edge_wr_buf_sz = 2 * edges_per_block;
  uint32_t *edge_wr_buf = new uint32_t[edge_wr_buf_sz * sizeof(uint32_t)];
  size_t edge_wr_buf_off = 0;
  memset(edge_wr_buf, 0x0, edge_wr_buf_sz * sizeof(uint32_t));

  std::random_device rd;
  std::mt19937 prng(rd());
  std::uniform_real_distribution<double> edge_dist(0.0, 1.0);

  // Open the output file for writing.
  MPI_File_open(MPI_COMM_WORLD, fname_out, MPI_MODE_CREATE | MPI_MODE_RDWR,
		MPI_INFO_NULL, &mpi_fp_out);
  MPI_File_set_view(mpi_fp_out, 0, MPI_UNSIGNED, MPI_UNSIGNED, "native", MPI_INFO_NULL);
  // Iterate over each node and calculate whether to
  // add edges based on the given edge probability.
  for (uint32_t u = 1; u <= n; u++) {
    // Add an incident "edge" from u to itself to ensure
    // that u is included in the graph even in the (possibly)
    // rare event that it is a singleton component.
    // But only do this on one machine. And not the same machine.
    if (rank == (u % machines)) {
      edge_wr_buf[edge_wr_buf_off++] = u;
      edge_wr_buf[edge_wr_buf_off++] = u;
    }
    // Guard against buffer overrun by dumping.
    if (edge_wr_buf_off == edge_wr_buf_sz) {
      // Wait for previous file write if still ongoing.
      MPI_Wait(&dump_req, &dump_status);
      // Because of the random nature of edge assignment, it's easier to
      // let MPI do the file synchronization via a shared file pointer
      // instead of trying to predict a priori when buffers will fill.
      // This should exhibit a somewhat pipelined behavior, especially
      // with larger buffers that are multiples of OS block sizes.
      MPI_File_iwrite_shared(mpi_fp_out, edge_wr_buf, edge_wr_buf_off, MPI_UNSIGNED, &dump_req);
      // Reset buffer offset to avoid overflow.
      edge_wr_buf_off = 0;
    }
    // Iterate over potential neighbors. Neighbors v
    // with labels less than u are ignored since
    // their incident edges to u would already have
    // been recorded if they exist.
    for (uint32_t v = u + 1; v <= n; v+= machines) {
      if (v + rank <= n) {
	// Sample the Bernoulli random variable E_uv:
	//   1, if sample <= p
	//   0  if sample > p
	if (edge_dist(prng) <= p_edge) {
	  // Add (u,v + rank) to the machine's output buffer.
	  edge_wr_buf[edge_wr_buf_off++] = u;
	  edge_wr_buf[edge_wr_buf_off++] = v + rank;
	}
	// Guard against buffer overrun by dumping.
	if (edge_wr_buf_off == edge_wr_buf_sz) {
	  MPI_Wait(&dump_req, &dump_status);
	  MPI_File_iwrite_shared(mpi_fp_out, edge_wr_buf, edge_wr_buf_off, MPI_UNSIGNED, &dump_req);
	  edge_wr_buf_off = 0;
	}
      }
    }
  }

  // Dump remainder of partially filled buffer.
  if (edge_wr_buf_off) {
    MPI_Wait(&dump_req, &dump_status);
    MPI_File_iwrite_shared(mpi_fp_out, edge_wr_buf, edge_wr_buf_off, MPI_UNSIGNED, &dump_req);
  }

  // Wait for last write to finish before exiting.
  MPI_Wait(&dump_req, &dump_status);
  MPI_File_close(&mpi_fp_out);

  // Clean up.
  delete[] edge_wr_buf;
  MPI_Finalize();

  // Exit.
  return 0;
}
