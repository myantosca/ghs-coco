#include <cstdio>
#include <cstdint>
#include <cstring>
#include <sys/time.h>
#include <random>
#include <mpi.h>

int main(int argc, char *argv[]) {
  int rank, machines, edges_per_block = 1;
  double p_edge = 0;
  uint32_t n = 1;
  bool monte_carlo = false;

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
	  fprintf(stderr, "Invalid output block size (edges): %s\n", argv[a]);
	}
	return -1;
      }
    }
    // Monte Carlo edge generation (faster, but may not be completely correct)
    if (!strcmp("-m", argv[a])) {
      monte_carlo = true;
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

  struct drand48_data rand_state;
  struct timeval tv;
  gettimeofday(&tv, NULL);
  srand48_r(tv.tv_sec * 1000000 + tv.tv_usec, &rand_state);
  std::random_device rd;
  std::mt19937 prng(tv.tv_sec * 1000000 + tv.tv_usec);
  std::bernoulli_distribution edge_bern_dist(p_edge);
  // Open the output file for writing.
  MPI_File_open(MPI_COMM_WORLD, fname_out, MPI_MODE_CREATE | MPI_MODE_RDWR,
		MPI_INFO_NULL, &mpi_fp_out);
  MPI_File_set_view(mpi_fp_out, 0, MPI_UNSIGNED, MPI_UNSIGNED, "native", MPI_INFO_NULL);
  // Iterate over each node and calculate whether to
  // add edges based on the given edge probability.
  for (uint32_t i = 0; i < n; i+= machines) {
    uint32_t u = i + rank;
    // Only determine edges for u within the population size.
    // The machines outside the remainder if n is not
    // evenly divisible by k will remain idle.
    // The fixity of u is important for ensuring the
    // equivalence of the binomial distribution over
    // u's neighbors to the Bernoulli distribution of
    // each edge. By fixing u, we can keep PRNG state
    // consistent for the uniform distribution of edges
    // to all other nodes v.
    if (u < n) {
      // Add an incident "edge" from u to itself to ensure
      // that u is included in the graph even in the (possibly)
      // rare event that it is a singleton component.
      // Since each u is unique to some machine for the
      // edge generation, we can do this unconditionally.
      edge_wr_buf[edge_wr_buf_off++] = u;
      edge_wr_buf[edge_wr_buf_off++] = u;

      // Guard against buffer overrun by dumping.
      if (edge_wr_buf_off == edge_wr_buf_sz) {
	// Wait for previous file write if still ongoing.
	dump_req = MPI_REQUEST_NULL;
	// Because of the random nature of edge assignment, it's easier to
	// let MPI do the file synchronization via a shared file pointer
	// instead of trying to predict a priori when buffers will fill.
	// This should exhibit a somewhat pipelined behavior, especially
	// with larger buffers that are multiples of OS block sizes.
	MPI_File_write_shared(mpi_fp_out, edge_wr_buf, edge_wr_buf_off, MPI_UNSIGNED, &dump_status);
	// Reset buffer offset to avoid overflow.
	edge_wr_buf_off = 0;
      }

      if (monte_carlo) {
	// Monte Carlo edge generation (faster, but not necessarily correct).

	// Create the binomial distribution for quickly determining
	// the degree of u without having to do n - u Bernoulli tests.
	// We avoid nodes v < u to avoid duplicates.
	// The probability would have been handled in a previous u.
	// Thus, there is an assumption of direction on each edge, i.e.,
	// ∀ (u,v) . u < v, which preserves the property that the
	// edge probability is not over- or undercounted per edge.
	std::binomial_distribution<uint32_t> edge_dist(n - 1 - u, p_edge);
	// Determine the degree of u. This should have some variance
	// by virtue of the PRNG and the distribution depending on it.
	size_t edges = edge_dist(prng);

	std::uniform_int_distribution<uint32_t> vs(u + 1, n - 1);
	// "Pick" each v randomly from a uniform distribution over [u + 1, n).
	for (uint32_t e = 0; e < edges; e++) {
	  // NB: We are NOT guaranteed to avoid duplicates here
	  // because the PRNG could yield the same v more than once.
	  uint32_t v = vs(prng);
	  // Add (u,v) to the machine's output buffer.
	  edge_wr_buf[edge_wr_buf_off++] = u;
	  edge_wr_buf[edge_wr_buf_off++] = v;

	  // Guard against buffer overrun by dumping.
	  if (edge_wr_buf_off == edge_wr_buf_sz) {
	    MPI_File_write_shared(mpi_fp_out, edge_wr_buf, edge_wr_buf_off, MPI_UNSIGNED, &dump_status);
	  edge_wr_buf_off = 0;
	  }
	}
      }
      else {
	// Las-Vegas edge generation (slower, but guaranteed correct).
	for (uint32_t v = u + 1; v < n; v++) {
	  if (edge_bern_dist(prng)) {
	    // Add (u,v) to the machine's output buffer.
	    edge_wr_buf[edge_wr_buf_off++] = u;
	    edge_wr_buf[edge_wr_buf_off++] = v;

	    // Guard against buffer overrun by dumping.
	    if (edge_wr_buf_off == edge_wr_buf_sz) {
	      MPI_File_write_shared(mpi_fp_out, edge_wr_buf, edge_wr_buf_off, MPI_UNSIGNED, &dump_status);
	      edge_wr_buf_off = 0;
	    }
	  }
	}
      }
    }
  }

  // Dump remainder of partially filled buffer.
  if (edge_wr_buf_off) {
    //MPI_Wait(&dump_req, &dump_status);
    MPI_File_write_shared(mpi_fp_out, edge_wr_buf, edge_wr_buf_off, MPI_UNSIGNED, &dump_status);
    edge_wr_buf_off = 0;
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
