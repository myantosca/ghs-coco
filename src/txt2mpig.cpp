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
  //MPI_File_set_view(mpi_fp_out, 0, MPI_UNSIGNED, MPI_UNSIGNED, "native", MPI_INFO_NULL);
  size_t edge_rd_buf_off = 0;
  size_t mpi_fp_out_off = 0;
  int c = 0;
  // Iterate over input text file.
  while (!eof) {
    if (rank == 0) {
      uint32_t u = 0;
      uint32_t v = 0;
      int found = fscanf(fp_in, "%u %u", &u, &v);
      if (found == 2) {
	// Line contains exactly two unsigned numbers
	// representing the endpoints of the edge.
	edge_rd_buf[edge_rd_buf_off++] = u;
	edge_rd_buf[edge_rd_buf_off++] = v;
	//fprintf(stderr, "%d: (%u, %u)\n", c++, u, v);
      }
      else {
	fscanf(fp_in, "%*s\n");
      }
      if (edge_rd_buf_off == edge_rd_buf_sz) {
	MPI_Bcast(&edge_rd_buf_off, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }
      if (eof = feof(fp_in)) {
	MPI_Bcast(&eof, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int last_machine = edge_rd_buf_off / (machines * edge_wr_buf_sz) + 1;
	for (int machine = last_machine + 1; machine < machines; machine++) {
	  edge_wr_buf_displs[machine] = 0;
	  edge_wr_buf_counts[machine] = 0;
	}
	if (last_machine < machines) {
	  edge_wr_buf_counts[last_machine] = edge_rd_buf_off % (machines * edge_wr_buf_sz);
	}
      }
    }

    // Check whether read buffer is full or EOF.
    if ((edge_rd_buf_off == edge_rd_buf_sz) || eof) {
      int rcvd = 0;
      MPI_Scatter(edge_wr_buf_counts, 1, MPI_UNSIGNED,
		  &rcvd, 1, MPI_UNSIGNED,
		  0, MPI_COMM_WORLD);
      // Dump buffer to disk.
      MPI_Scatterv(edge_rd_buf, edge_wr_buf_counts, edge_wr_buf_displs, MPI_UNSIGNED,
		   edge_wr_buf, edge_wr_buf_sz, MPI_UNSIGNED,
		   0, MPI_COMM_WORLD);
      if (rcvd > 0) {
	// @FIXME: mpirun -np > 1 hangs.
	MPI_File_set_view(mpi_fp_out, (mpi_fp_out_off + rank * edge_wr_buf_sz) * sizeof(uint32_t),
			  MPI_UNSIGNED, MPI_UNSIGNED, "native", MPI_INFO_NULL);
	MPI_Request dump_req;
	MPI_File_iwrite_at(mpi_fp_out, 0, edge_wr_buf, rcvd, MPI_UNSIGNED, &dump_req);
	MPI_Wait(&dump_req, NULL);
      }

      mpi_fp_out_off += edge_rd_buf_off;
      edge_rd_buf_off = 0;
    }

  }

  // Close file handles.
  if (rank == 0) {
    fclose(fp_in);
  }
  MPI_File_close(&mpi_fp_out);
  // Tear down MPI.
  MPI_Finalize();
  return 0;
}
