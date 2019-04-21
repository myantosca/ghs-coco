#include <cstdio>
#include <cstdint>
#include <cstring>
#include <sys/time.h>
#include <random>
#include <map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <mpi.h>

// https://en.wikipedia.org/wiki/Mersenne_prime#List_of_known_Mersenne_primes
#define MERSENNE_61 2305843009213693951
#define MACHINE_HASH(__VERTEX__) ((hash_a * __VERTEX__ + hash_b) % MERSENNE_61) % machines

typedef enum {
  UNGROUPED = 0,
  BROADCAST = 1,
  PENDING = 2,
  FINISHED = 3
} vertex_state_t;


typedef struct vertex
{
  uint32_t id;
  uint32_t parent;
  uint32_t group;
  uint32_t group_ct;
  vertex_state_t state;
  std::unordered_set<uint32_t> neighbors;
} vertex_t;

typedef struct ucast_msg
{
  uint32_t parent;
  uint32_t child;
  uint32_t group_ct;
} ucast_msg_t;

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
  int verbosity = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &machines);
  // Default input filename.
  const char *fname_in = "./graph.ecg";

  int a = 0;
  while (a < argc) {
    // Input filename argument.
    if (!strcmp("-i", argv[a])) {
      a++;
      if (a < argc) fname_in = argv[a];
    }
    if (!strcmp("-v", argv[a])) {
      verbosity = 1;
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

  uint32_t *realloc_temp;
  // Fill send buffers for edge exchange.
  for (int i = 0; i < edges; i++) {
    // Get u and v from the file buffer.
    uint32_t u = edges_buf[i*2];
    uint32_t v = edges_buf[i*2+1];
    // Determine the machines housing u and v.
    uint64_t machine_u = MACHINE_HASH(u);
    uint64_t machine_v = MACHINE_HASH(v);
    // Resize machine[u]'s buffer if we're out of room.
    if (send_counts[machine_u] == send_caps[machine_u]) {
      send_caps[machine_u] *= 2;
      realloc_temp = (uint32_t *)realloc(send_bufs[machine_u],
					 send_caps[machine_u] * sizeof(uint32_t));
      if (realloc_temp) {
	send_bufs[machine_u] = realloc_temp;
      }
      else {
	exit(errno);
      }
    }
    // Add (u,v) to buffer to send to machine[u].
    // By convention, the node housed in the machine is first.
    send_bufs[machine_u][send_counts[machine_u]++] = u;
    send_bufs[machine_u][send_counts[machine_u]++] = v;
    // Resize machine[v]'s buffer if we're out of room.
    if (send_counts[machine_v] == send_caps[machine_v]) {
      send_caps[machine_v] *= 2;
      realloc_temp = (uint32_t *)realloc(send_bufs[machine_v],
					 send_caps[machine_v] * sizeof(uint32_t));
      if (realloc_temp) {
	send_bufs[machine_v] = realloc_temp;
      }
      else {
	exit(errno);
      }
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
      vertex_u->state = UNGROUPED;
      vertex_u->neighbors = std::unordered_set<uint32_t>();
      vertex_u->neighbors.clear();
      V_machine[u] = vertex_u;
    }
    if (u != v) {
      V_machine[u]->neighbors.insert(v);
    }
  }

  free( recv_buf );

  if (verbosity == 1) {
    // Debug printout of graph in vertex-centric format.
    for (int machine = 0; machine < machines; machine++) {
      if (rank == machine) {
	for (auto &kv : V_machine) {
	  std::cout << "[" << machine << "]";
	  std::cout << kv.first << ":";
	  for (auto &neighbor : kv.second->neighbors) {
	    std::cout << " " << neighbor;
	  }
	  std::cout << std::endl;
	}
      }
    }
  }

  // Execute BFS search for forest.
  uint32_t forest_sz = 2 * sizeof(uint32_t);;
  uint32_t trees = 0;
  uint32_t *forest = (uint32_t *)malloc(forest_sz);
  memset(forest, 0, forest_sz);
  int global_done = 0;
  uint32_t finished = 0;
  int local_done = V_machine.size() == finished;
  std::map <int, std::unordered_set<uint32_t>> bcast_msgs;
  std::map <int, std::vector<ucast_msg_t>> ucast_msgs;
  std::unordered_set<uint32_t> to_delete;
  for (int machine = 0; machine < machines; machine++) {
    bcast_msgs[machine] = std::unordered_set<uint32_t>();
    ucast_msgs[machine] = std::vector<ucast_msg_t>();
  }
  MPI_Allreduce(&local_done, &global_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
  // Continue reducing vertices to forest
  // until all internal vertices are marked finished.
  while (!global_done) {
    // Out of room for forest!
    uint32_t trees_off = trees << 3;
    uint32_t bfs_root;
    if (trees_off == forest_sz) {
      // Double the buffer size.
      forest_sz = forest_sz << 1;
      realloc_temp = (uint32_t *)realloc(forest, forest_sz);
      if (realloc_temp) {
	forest = realloc_temp;
      }
      else {
	exit(errno);
      }
    }
    // // The first element must be ungrouped. Otherwise,
    // // it would have already been removed from the map.
    // bfs_root = V_machine.empty() ? (1 << 31) : V_machine.begin()->first;
    bfs_root = (1 << 31);
    for (auto &kv : V_machine) {
      if (kv.second->state == UNGROUPED) {
	bfs_root = kv.second->id;
	break;
      }
    }
    forest[trees*2 + 1] = 0;
    // Elect the minimum vertex id as the BFS tree root.
    MPI_Allreduce(&bfs_root, &forest[trees*2],
		  1, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);
    // if (bfs_root == 1 << 31) break;
    bfs_root = forest[trees*2];
    int bfs_root_machine = MACHINE_HASH(bfs_root);
    // Set the root node to broadcast state.
    if (rank == bfs_root_machine) {
      vertex_t *r = V_machine[bfs_root];
      r->state = BROADCAST;
    }
    uint32_t tree_done = 0;
    // Continue flooding until root has received
    // population subtotals from all its children.
    while (!forest[trees*2+1]) {
      /*******************
       * Broadcast phase *
       *******************/
      for (auto &kv : V_machine) {
	vertex_t *u = kv.second;
	// If scheduled to broadcast...
	if (u->state == BROADCAST) {
	  u->state = PENDING;
	  // Broadcast group to children.
	  to_delete.clear();
	  for (auto &child : u->neighbors) {
	    uint32_t child_machine = MACHINE_HASH(child);
	    // Local delivery
	    if (rank == child_machine) {
	      vertex_t *v = V_machine[child];
	      if (v->state == UNGROUPED) {
		// If ungrouped, assign parent.
		v->parent = u->id;
		// Assign group label.
		v->group = u->group;
		// Update vertex state.
		v->state = BROADCAST;
	      }
	      else {
		to_delete.insert(child);
	      }
	      // Remove attempted parent regardless of prospective child state.
	      v->neighbors.erase(u->id);
	    }
	    // Remote delivery
	    else {
	      // Add sender to machine-specific broadcast pre-buffer.
	      bcast_msgs[child_machine].insert(u->id);
	    }
	  }
	  for(auto &child : to_delete) {
	    u->neighbors.erase(child);
	  }
	  to_delete.clear();
	}
      }

      // Gather machine-specific broadcast buffers for each machine in turn.
      for (auto &machine_msg : bcast_msgs) {
	send_counts[machine_msg.first] = machine_msg.second.size();
	send_bufs[machine_msg.first] = (uint32_t *)malloc(send_counts[machine_msg.first] * sizeof(uint32_t));
	int i = 0;
	for (auto &v : machine_msg.second) {
	  send_bufs[machine_msg.first][i++] = v;
	}
	machine_msg.second.clear();
      }

      // Exchange broadcast messages between all machines.
      exchange(rank, machines,
	       send_counts, MPI_UNSIGNED, send_bufs,
	       recv_counts, recv_displs, recv_totals,
	       MPI_UNSIGNED, &recv_buf);

      // Perform local receipt of remote flood.
      for (int i = 0; i < recv_totals[rank]; i++) {
	uint32_t id_u = recv_buf[i];
	uint32_t machine_u = MACHINE_HASH(id_u);
	// Check each node for connection to each remote flooding parent.
	for (auto &kv : V_machine) {
	  vertex_t *v = kv.second;
	  // If connected...
	  if (v->neighbors.find(id_u) != v->neighbors.end()) {
	    // Child is ungrouped. Meet your parent!
	    if (v->state == UNGROUPED) {
	      v->parent = id_u;
	      v->state = BROADCAST;
	      v->group = bfs_root;
	    }
	    // Child already has a parent. Upcast to remove dead link.
	    else {
	      // No local check is done here since the parent must be remote.
	      ucast_msg_t msg;
	      msg.parent = id_u;
	      msg.child = v->id;
	      msg.group_ct = 0;
	      ucast_msgs[machine_u].push_back(msg);
	    }
	    // Remove attempted parent regardless of prospective child state.
	    v->neighbors.erase(id_u);
	  }
	}
      }
      if (recv_buf) free(recv_buf);

      /****************
       * Upcast phase *
       ****************/
      for (auto &kv : V_machine) {
	vertex_t *u = kv.second;
	// If no longer waiting on children...
	if (u->state == PENDING && u->neighbors.size() == 0) {
	  u->state == FINISHED;
	  // Upcast subtree population to parent.
	  if (u->id == bfs_root) {
	    // Update the tree count.
	    forest[trees*2+1] = u->group_ct;
	  }
	  else {
	    uint32_t parent_machine = MACHINE_HASH(u->parent);
	    // Local delivery.
	    if (rank == parent_machine) {
	      // Subsume group population into parent node.
	      V_machine[u->parent]->group_ct += u->group_ct;
	      // Decrement parent's awaiting counter.
	      V_machine[u->parent]->neighbors.erase(u->id);
	    }
	    // Remote delivery.
	    else {
	      ucast_msg_t msg;
	      msg.parent = u->parent;
	      msg.child = u->id;
	      msg.group_ct = u->group_ct;
	      ucast_msgs[parent_machine].push_back(msg);
	    }
	  }
	  to_delete.insert(kv.first);
	}
      }

      // Gather machine-specific upcast buffers for each machine in turn.
      for (auto &machine_msgs : ucast_msgs) {
	send_counts[machine_msgs.first] = machine_msgs.second.size() * 3;
	send_bufs[machine_msgs.first] = (uint32_t *)malloc(send_counts[machine_msgs.first] * sizeof(uint32_t));
	int i = 0;
	for (auto &msg: machine_msgs.second) {
	  send_bufs[machine_msgs.first][i++] = msg.parent;
	  send_bufs[machine_msgs.first][i++] = msg.child;
	  send_bufs[machine_msgs.first][i++] = msg.group_ct;
	}
	machine_msgs.second.clear();
      }

      // Exchange broadcast messages between all machines.
      exchange(rank, machines,
	       send_counts, MPI_UNSIGNED, send_bufs,
	       recv_counts, recv_displs, recv_totals,
	       MPI_UNSIGNED, &recv_buf);

      // Perform local receipt of remote upcast.
      for (int i = 0; i < recv_totals[rank]; i+=3) {
	uint32_t parent = recv_buf[i];
	uint32_t child = recv_buf[i+1];
	uint32_t group_ct = recv_buf[i+2];
	V_machine[parent]->neighbors.erase(child);
	V_machine[parent]->group_ct += group_ct;
      }
      if (recv_buf) free(recv_buf);

      // Erase moribund vertices that have sent their
      // respective upcast messages.
      for (auto &v : to_delete) {
	if (V_machine[v]) free(V_machine[v]);
	V_machine.erase(v);
      }
      to_delete.clear();
      // Reduce termination condition.
      MPI_Allreduce(&forest[trees*2 + 1], &forest[trees*2 + 1], 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
    }
    trees++;
    local_done = V_machine.empty();
    MPI_Allreduce(&local_done, &global_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
  }

  if (rank == 0) {
    for (int tree = 0; tree < trees; tree+=2) {
      std::cout << forest[tree] << "," << forest[tree+1] << std::endl;
    }
  }

  // Clean up.
  bcast_msgs.clear();
  ucast_msgs.clear();

  if (send_counts) free( send_counts );
  if (recv_counts) free( recv_counts );
  if (recv_displs) free( recv_displs );
  if (recv_totals) free( recv_totals );

  for (auto &kv : V_machine) {
    if (kv.second) free(kv.second);
    kv.second = NULL;
  }
  V_machine.clear();
  free(forest);
  // Tear down MPI.
  MPI_Finalize();

  // Exit.
  return 0;
}
