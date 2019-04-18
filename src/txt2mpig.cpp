#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

int main(int argc, char *argv[]) {
  int rank, machines, edges_per_block = 1;
  // Initialize MPI.
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &machines);

  // Default input filename.
  const char *fname_in = "./graph.txt";
  // Default output filename.
  const char *fname_out = "./graph.ecg";
  int a = 0;
  while (a < argc) {
    // Input filename argument.
    if (!strcmp("-i", argv[a])) {
      a++;
      if (a < argc) fname_in = argv[a];
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

  FILE* fp_in = NULL;
  MPI_File mpi_fp_out;
  size_t mpi_fp_out_off = 0;

  size_t edge_wr_buf_sz = 2 * edges_per_block;
  size_t edge_rd_buf_sz = edge_wr_buf_sz * machines;
  uint32_t *edge_rd_buf = NULL;
  uint32_t *edge_wr_buf = new uint32_t[edge_wr_buf_sz];
  memset(edge_wr_buf, 0x0, edge_wr_buf_sz * sizeof(uint32_t));
  int *edge_wr_buf_counts = NULL;
  int *edge_wr_buf_displs = NULL;
  int eof = 0;
  if (rank == 0) {
    fp_in = fopen(fname_in, "r");
    eof = feof(fp_in);
    edge_rd_buf = new uint32_t[edge_rd_buf_sz];
    memset(edge_rd_buf, 0x0, edge_rd_buf_sz * sizeof(uint32_t));
    edge_wr_buf_counts = new int[machines];
    edge_wr_buf_displs = new int[machines];
    for (int machine = 0; machine < machines; machine++) {
      edge_wr_buf_displs[machine] = machine * edge_wr_buf_sz;
      edge_wr_buf_counts[machine] = (machine + 1) * edge_wr_buf_sz;
    }
  }

  MPI_File_open(MPI_COMM_WORLD, fname_out, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &mpi_fp_out);
  MPI_File_set_view(mpi_fp_out, 0, MPI_UNSIGNED, MPI_UNSIGNED, "native", MPI_INFO_NULL);
  // Read buffer offset.
  size_t edge_rd_buf_off = 0;

  // Temporary variable for line of text.
  char *line = NULL;
  // Input text line size.
  size_t line_sz = 0;

  MPI_Request dump_req = MPI_REQUEST_NULL;
  // Iterate over input text file.
  while (!eof) {
    if (rank == 0) {
      uint32_t u = 0;
      uint32_t v = 0;
      // Read next line from input file.
      getline(&line, &line_sz, fp_in);
      // Check for two unsigned numbers separated by whitespace
      // representing the endpoints of the edge in the retrieved
      // line of text.
      int found = sscanf(line, "%u %u", &u, &v);
      if (found == 2) {
	edge_rd_buf[edge_rd_buf_off++] = u;
	edge_rd_buf[edge_rd_buf_off++] = v;
      }
      eof = feof(fp_in);
    }

    // Update all other processes on read buffer offset.
    MPI_Bcast(&edge_rd_buf_off, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    // Update all other processes on EOF status.
    MPI_Bcast(&eof, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Check whether read buffer is full or EOF.
    if ((edge_rd_buf_off == edge_rd_buf_sz) || eof) {
      // Wait on file dump before overwriting the write buffer.
      MPI_Wait(&dump_req, NULL);
      // Reset dump request to avoid segfault.
      dump_req = MPI_REQUEST_NULL;
      // Divide the input among all the processes.
      MPI_Scatter(edge_rd_buf, edge_wr_buf_sz, MPI_UNSIGNED,
      		  edge_wr_buf, edge_wr_buf_sz, MPI_UNSIGNED,
      		  0, MPI_COMM_WORLD);
      // Number of processes participating in the actual write.
      int participants = edge_rd_buf_off / edge_wr_buf_sz;
      // Buffer size to be written by a given process.
      int edge_wr_buf_count =
	(rank < participants                     // rank well below the offset cutoff
	 ? edge_wr_buf_sz                        //   Write away!
	 : (rank == participants                 // rank on the cutoff border
	    ? edge_rd_buf_off % edge_wr_buf_sz   //   If remainder exists, write that.
	    : 0));                               //   Otherwise, nothing.
      // Only write if necessary.
      if (edge_wr_buf_count > 0) {
	MPI_File_iwrite_at(mpi_fp_out, mpi_fp_out_off + rank * edge_wr_buf_sz,
			   edge_wr_buf, edge_wr_buf_count, MPI_UNSIGNED, &dump_req);
      }

      // Increment absolute write offset.
      mpi_fp_out_off += edge_rd_buf_off;
      // Reset read buffer offset.
      edge_rd_buf_off = 0;
    }

  }

  // Wait on file dump, just in case.
  MPI_Wait(&dump_req, NULL);

  // Clean up.
  if (rank == 0) {
    fclose(fp_in);
    delete[] edge_rd_buf;
    delete[] edge_wr_buf_counts;
    delete[] edge_wr_buf_displs;
  }

  delete[] edge_wr_buf;
  MPI_File_close(&mpi_fp_out);

  // Tear down MPI.
  MPI_Finalize();

  // Exit.
  return 0;
}
