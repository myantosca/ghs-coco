#include <cstdio>
#include <cstdint>
#include <cstring>
#include <sys/time.h>
#include <random>
#include <map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <chrono>
#include <mpi.h>

using namespace std::chrono;

// https://en.wikipedia.org/wiki/Mersenne_prime#List_of_known_Mersenne_primes
#define MERSENNE_61 2305843009213693951
#define MACHINE_HASH(__VERTEX__) ((hash_a * __VERTEX__ + hash_b) % MERSENNE_61) % machines

typedef enum {
  UNGROUPED = 0,
  BROADCAST = 1,
  PENDING = 2,
  FINISHED = 3
} vertex_state_t;


typedef struct exchange_info
{
  int rank;
  int machines;
  int *send_counts;
  int *send_caps;
  uint32_t **send_bufs;
  int *recv_counts;
  int *recv_displs;
  int recv_caps;
  uint32_t *recv_buf;
} exchange_info_t;

exchange_info_t *exchange_info_new(int rank, int machines) {
  exchange_info_t *info = (exchange_info_t *)malloc(sizeof(exchange_info_t));
  info->rank = rank;
  info->machines = machines;
  info->send_counts = (int *)malloc(info->machines * sizeof(int));
  memset(info->send_counts, 0, info->machines * sizeof(int));
  info->send_caps = (int *)malloc(info->machines * sizeof(int));
  memset(info->send_caps, 0, info->machines * sizeof(int));
  info->send_bufs = (uint32_t**)malloc(info->machines * sizeof(uint32_t *));
  memset(info->send_bufs, 0, info->machines * sizeof(uint32_t *));
  info->recv_counts = (int *)malloc(info->machines * sizeof(int));
  memset(info->recv_counts, 0, info->machines * sizeof(int));
  info->recv_displs = (int *)malloc(info->machines * sizeof(int));
  memset(info->recv_displs, 0, info->machines * sizeof(int));
  info->recv_caps = 0;
  info->recv_buf = NULL;
  return info;
}

void exchange_info_free(exchange_info_t *info) {
  if (info) {
    if (info->send_counts) { free(info->send_counts); }
    if (info->send_caps)   { free(info->send_caps); }
    if (info->send_bufs)   { free(info->send_bufs); }
    if (info->recv_counts) { free(info->recv_counts); }
    if (info->recv_displs) { free(info->recv_displs); }
    if (info->recv_buf)    { free(info->recv_buf); }
    free(info);
  }
}

void exchange_info_send_buf_insert(exchange_info_t *info, int machine, uint32_t *vals, int count) {
  if (info->send_counts[machine] + count > info->send_caps[machine]) {
    info->send_caps[machine] += count * 2;
    uint32_t *tmp = (uint32_t *)realloc(info->send_bufs[machine], info->send_caps[machine] * sizeof(uint32_t));
    if (tmp) {
      info->send_bufs[machine] = tmp;
    }
    else {
      exit(errno);
    }
  }
  memcpy(&(info->send_bufs[machine][info->send_counts[machine]]), vals, count * sizeof(uint32_t));
  info->send_counts[machine] += count;
}

void exchange_info_send_buf_resize(exchange_info_t *info, int machine, int caps) {
  if (caps > info->send_caps[machine]) {
    uint32_t *tmp = (uint32_t *)realloc(info->send_bufs[machine], caps * sizeof(uint32_t));
    if (tmp) {
      info->send_bufs[machine] = tmp;
      info->send_caps[machine] = caps;
    }
    else {
      exit(errno);
    }
  }
}

void exchange_info_rewind(exchange_info_t *info) {
  for (int machine = 0; machine < info->machines; machine++) {
    info->send_counts[machine] = 0;
  }
}

typedef struct vertex
{
  uint32_t id;
  uint32_t parent;
  uint32_t group;
  uint32_t group_ct;
  vertex_state_t state;
  uint32_t awaiting;
  std::unordered_set<uint32_t> neighbors;
  std::unordered_set<uint32_t> children;
} vertex_t;

typedef struct ucast_msg
{
  uint32_t parent;
  uint32_t child;
  uint32_t group_ct;
} ucast_msg_t;

void exchange(exchange_info_t *info, uint64_t *messages) {
  memset(info->recv_counts, 0, sizeof(int) * info->machines);
  memset(info->recv_displs, 0, sizeof(int) * info->machines);
  // Gather for each machine in turn from every other machine's
  // targeted send buffers.
  int received = 0;
  for (int machine = 0; machine < info->machines; machine++) {
    // Determine the receive sub-buffer element counts.
    MPI_Gather(&info->send_counts[machine], 1, MPI_INT,
               info->recv_counts, 1, MPI_INT,
               machine, MPI_COMM_WORLD);
    *messages += info->machines;
    // Only one machine is target of the gather operation at a time.
    if (info->rank == machine) {
      for (int i = 1; i < info->machines; i++) {
        info->recv_displs[i] = info->recv_displs[i-1] + info->recv_counts[i-1];
      }
      // Calculate the total receive buffer size so each machine knows its extent.
      info->recv_caps = info->recv_counts[info->machines - 1] + info->recv_displs[info->machines - 1];
      received = info->recv_caps;
      // Allocate the receive buffer according to the calculated size.
      info->recv_buf = (uint32_t *)malloc(info->recv_caps * sizeof(uint32_t));
      // Zero out the receive buffer as a matter of safety.
      memset(info->recv_buf, 0, info->recv_caps * sizeof(uint32_t));
    }
    // Make sure everyone knows the total messages in this gather for accounting.
    MPI_Bcast(&received, 1, MPI_INT, machine, MPI_COMM_WORLD);
    *messages += info->machines;
    // Gather all the targeted send buffers on the target machine
    // in a contiguous array in order by sender rank.
    MPI_Gatherv(info->send_bufs[machine], info->send_counts[machine], MPI_UNSIGNED,
                info->recv_buf, info->recv_counts, info->recv_displs, MPI_UNSIGNED,
                machine, MPI_COMM_WORLD);
    *messages += received;
  }
}

int main(int argc, char *argv[]) {
  int rank, machines;
  uint32_t local_vertices = 0, global_vertices = 0;
  uint64_t messages = 0, rounds = 0;
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

  // Capture starting time for preprocessing.
  time_point<system_clock, nanoseconds> tp_prep_a = system_clock::now();

  // Determine hash function parameters.
  std::random_device rd;
  std::mt19937 prng(rd());

  uint64_t hash_a = 1;
  uint64_t hash_b = 0;
  if (rank == 0) {
    std::uniform_int_distribution<uint64_t> hash_a_dist(1, MERSENNE_61 - 1);
    std::uniform_int_distribution<uint64_t> hash_b_dist(0, MERSENNE_61 - 1);
    hash_a = hash_a_dist(prng);
    hash_b = hash_b_dist(prng);
  }

  // Distribute hash parameters to all machines.
  MPI_Bcast(&hash_a, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hash_b, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

  exchange_info_t *edges_xinfo = exchange_info_new(rank, machines);
  exchange_info_t *bcast_xinfo = exchange_info_new(rank, machines);
  exchange_info_t *ucast_xinfo = exchange_info_new(rank, machines);

  std::vector<std::unordered_set<uint32_t>> bcast_msgs;
  std::map<uint32_t, ucast_msg_t> ucast_msgs;
  std::vector<vertex_t *>bcast_vertices;
  std::vector<vertex_t *>ucast_vertices;
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
  size_t mpi_fp_in_off = 0;

  // Initialize send buffers for edge exchange.
  for (int machine = 0; machine < machines; machine++) {
    exchange_info_send_buf_resize(edges_xinfo, machine, edges >> 1);
    edges_xinfo->send_counts[machine] = 0;
    // Increment the file read start offset for this rank.
    if (machine < rank) {
      mpi_fp_in_off += (edges_per_machine + (machine < edges_leftover ? 1 : 0)) << 3;
    }
    bcast_msgs.push_back(std::unordered_set<uint32_t>());
  }

  // Fill edge input buffers.
  MPI_File_set_view(mpi_fp_in, mpi_fp_in_off, MPI_UNSIGNED, MPI_UNSIGNED, "native", MPI_INFO_NULL);
  MPI_File_read(mpi_fp_in, edges_buf, edges << 1, MPI_UNSIGNED, NULL);
  MPI_File_close(&mpi_fp_in);

  uint32_t *realloc_temp;
  // Fill send buffers for edge exchange.
  for (int i = 0; i < edges; i++) {
    // Get u and v from the file buffer.
    uint32_t uv[2];
    uv[0] = edges_buf[i*2];
    uv[1] = edges_buf[i*2+1];
    uint32_t vu[2];
    vu[0] = uv[1];
    vu[1] = uv[0];
    // Determine the machines housing u and v.
    uint64_t machine_u = MACHINE_HASH(uv[0]);
    uint64_t machine_v = MACHINE_HASH(vu[0]);
    // Add (u,v) to buffer to send to machine[u].
    // By convention, the node housed in the machine is first.
    exchange_info_send_buf_insert(edges_xinfo, machine_u, uv, 2);

    // Add (v,u) to buffer to send to machine[v].
    // By convention, the node housed in the machine is first.
    exchange_info_send_buf_insert(edges_xinfo, machine_v, vu, 2);
  }

  // Clean up so as to not hoard memory.
  free(edges_buf);

  // Exchange edges.
  exchange(edges_xinfo, &messages);

  // Construct local vertex-centric model.
  std::map<uint32_t, vertex_t *> V_in;
  std::vector<vertex_t *> V_out;
  std::map<uint32_t, std::unordered_set<uint32_t>> E_incoming;
  for (int i = 0; i < edges_xinfo->recv_caps; i += 2) {
    uint32_t u = edges_xinfo->recv_buf[i];
    uint32_t v = edges_xinfo->recv_buf[i+1];
    vertex_t *vertex_u = V_in[u];
    // If the node has not been created yet, do so.
    if (!vertex_u) {
      vertex_u = (vertex_t *)malloc(sizeof(vertex_t));
      memset(vertex_u, 0, sizeof(vertex_t));
      vertex_u->id = u;
      vertex_u->parent = u;
      vertex_u->group = u;
      vertex_u->group_ct = 1; // population of BFS sub-tree including itself
      vertex_u->state = UNGROUPED;
      vertex_u->awaiting = 0;
      vertex_u->neighbors = std::unordered_set<uint32_t>();
      vertex_u->neighbors.clear();
      vertex_u->children = std::unordered_set<uint32_t>();
      vertex_u->children.clear();
      V_in[u] = vertex_u;
    }
    if (u != v) {
      vertex_u->neighbors.insert(v);
      if (E_incoming.find(v) == E_incoming.end()) { E_incoming[v] = std::unordered_set<uint32_t>(); }
      E_incoming[v].insert(u);
    }
  }

  // Initialize awaiting counter here to avoid double counting
  // edges if graph storage has duplicates.
  for (auto &kv : V_in) {
    kv.second->awaiting = kv.second->neighbors.size();
  }

  // Clean up so as not hoard memory.
  exchange_info_free(edges_xinfo);

  // Get local stats.
  local_vertices = V_in.size();

  // Capture finishing time for preprocessing.
  time_point<system_clock, nanoseconds> tp_prep_b = system_clock::now();

  if (verbosity == 1) {
    // Debug printout of graph in vertex-centric format.
    for (int machine = 0; machine < machines; machine++) {
      if (rank == machine) {
        for (auto &kv : V_in) {
          std::cerr << "[" << machine << "]";
          std::cerr << kv.first << ":";
          for (auto &neighbor : kv.second->neighbors) {
            std::cerr << " " << neighbor;
          }
          std::cerr << std::endl;
        }
      }
    }
  }

  // Capture starting time for connected component search.
  time_point<system_clock, nanoseconds> tp_coco_a = system_clock::now();

  // Execute BFS search for forest.
  uint32_t forest_sz = 2 * sizeof(uint32_t);;
  uint32_t trees = 0;
  uint32_t tree_popl_min = 1 << 31;
  uint32_t tree_popl_max = 0;
  double tree_popl_av1 = 0;
  double tree_popl_av0 = 0;
  double tree_popl_var = 0;
  double tree_popl_ess = 0;
  uint32_t *forest = (uint32_t *)malloc(forest_sz);
  memset(forest, 0, forest_sz);
  int global_done = 0;
  int local_done = V_in.empty();
  std::unordered_set<uint32_t> finished;

  MPI_Allreduce(&local_done, &global_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
  messages += machines << 1;
  // Continue reducing vertices to forest until all internal vertices
  // have been moved from the input map to the output map.
  uint32_t trees_off;
  uint32_t bfs_root;
  while (!global_done) {
    trees_off = trees << 3;
    // Out of room for forest!
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

    // The first element must be ungrouped. Otherwise,
    // it would have already been removed from the map.
    bfs_root = V_in.empty() ? (1 << 31) : V_in.begin()->second->id;
    forest[trees*2 + 1] = 0;
    // Elect the minimum vertex id as the BFS tree root.
    MPI_Allreduce(&bfs_root, &forest[trees*2],
                  1, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);
    messages += machines << 1;
    bfs_root = forest[trees*2];
    int bfs_root_machine = MACHINE_HASH(bfs_root);
    // Set the root node to broadcast state.
    if (rank == bfs_root_machine) {
      vertex_t *r = V_in[bfs_root];
      r->state = BROADCAST;
      bcast_vertices.push_back(r);
    }
    uint32_t tree_popl = 0;
    // Continue flooding until root has received
    // population subtotals from all its children.
    while (!tree_popl) {
      /*******************
       * Broadcast phase *
       *******************/
      for (vertex *u : bcast_vertices) {
        // If scheduled to broadcast...
        if (u->state == BROADCAST) {
          u->state = PENDING;
          for (auto &id_v : u->neighbors) {
            uint32_t machine_v = MACHINE_HASH(id_v);
	    bcast_msgs[machine_v].insert(u->id);
          }
	  // Be sure to add singleton nodes to the upcast queue.
	  if (u->awaiting == 0) { ucast_vertices.push_back(u); }
        }
      }

      for (int machine = 0; machine < machines; machine++) {
	for (auto &u : bcast_msgs[machine]) {
	  uint32_t id_u = u;
	  exchange_info_send_buf_insert(bcast_xinfo, machine, &id_u, 1);
	}
      }

      // Exchange broadcast messages between all machines.
      exchange(bcast_xinfo, &messages);

      // Reset, rewind.
      exchange_info_rewind(bcast_xinfo);
      for (auto &msg : bcast_msgs) {
	msg.clear();
      }
      bcast_vertices.clear();

      // Broadcast receipt.
      for (int i = 0; i < bcast_xinfo->recv_caps; i++) {
        uint32_t id_u = bcast_xinfo->recv_buf[i];
        uint32_t machine_u = MACHINE_HASH(id_u);
	// Use reverse lookup to only hit nodes that will have received something.
	for (auto &id_v : E_incoming[id_u]) {
	  vertex_t *v = V_in[id_v];
          // If connected...
          if (v->neighbors.find(id_u) != v->neighbors.end()) {
            // Child is ungrouped. Meet your parent!
            if (v->state == UNGROUPED) {
              v->parent = id_u;
              v->state = BROADCAST;
              v->group = bfs_root;
	      bcast_vertices.push_back(v);
            }
            // Child already has a parent. Upcast to remove dead link.
            else {
              ucast_msg_t msg;
              msg.parent = id_u;
              msg.child = v->id;
              msg.group_ct = 0;
	      exchange_info_send_buf_insert(ucast_xinfo, machine_u, (uint32_t *)&msg, 3);
            }
          }
        }
      }

      /****************
       * Upcast phase *
       ****************/
      for (vertex_t *u : ucast_vertices) {
        // If no longer waiting on children...
	u->state = FINISHED;
	// Upcast subtree population to parent.
	if (u->id == bfs_root) {
	  // Update the tree count.
	  forest[trees*2+1] = u->group_ct;
	}
	else {
	  uint32_t parent_machine = MACHINE_HASH(u->parent);
	  ucast_msg_t msg;
	  msg.parent = u->parent;
	  msg.child = u->id;
	  msg.group_ct = u->group_ct;
	  exchange_info_send_buf_insert(ucast_xinfo, parent_machine, (uint32_t *)&msg, 3);
	}
	V_out.push_back(u);
	finished.insert(u->id);
      }

      ucast_vertices.clear();

      // Exchange upcast messages between all machines.
      exchange(ucast_xinfo, &messages);
      exchange_info_rewind(ucast_xinfo);

      // Perform local receipt of remote upcast.
      for (int i = 0; i < ucast_xinfo->recv_caps; i+=3) {
        uint32_t parent = ucast_xinfo->recv_buf[i];
        uint32_t child = ucast_xinfo->recv_buf[i+1];
        uint32_t group_ct = ucast_xinfo->recv_buf[i+2];
        vertex_t *p = V_in[parent];
        p->group_ct += group_ct;
        if (group_ct > 0) { p->children.insert(child); }
	p->awaiting--;
	// Add nodes with nothing to await to the upcast queue.
	if (p->awaiting == 0) { ucast_vertices.push_back(p); }
      }

      // Reduce termination condition. Implicit synchronization point.
      MPI_Allreduce(&forest[trees*2 + 1], &tree_popl, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
      messages += machines << 1;
      forest[trees*2 + 1] = tree_popl;
      // Move finished vertices to the output map.
      for (auto &v : finished) {
        V_in.erase(v);
      }
      finished.clear();
      rounds++;
    }

    // Component population statistics.
    trees++;
    tree_popl_min = (tree_popl < tree_popl_min) ? tree_popl : tree_popl_min;
    tree_popl_max = (tree_popl > tree_popl_max) ? tree_popl : tree_popl_max;
    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
    tree_popl_av0 = tree_popl_av1;
    tree_popl_av1 += (tree_popl - tree_popl_av0) / trees;
    tree_popl_ess += (tree_popl - tree_popl_av0) * (tree_popl - tree_popl_av1);
    tree_popl_var = tree_popl_ess / trees;

    // Global termination condition.
    local_done = V_in.empty();
    MPI_Allreduce(&local_done, &global_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    messages += machines << 1;
  }

  // Capture finishing time for connected component search.
  time_point<system_clock, nanoseconds> tp_coco_b = system_clock::now();
  int64_t T[2] = { (tp_prep_b - tp_prep_a).count(), (tp_coco_b - tp_coco_a).count() };
  int64_t T0[2] = { 0, 0 };
  MPI_Reduce(&T, &T0, 2, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_vertices, &global_vertices, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  // Print out connected component labels with component population counts.
  if (rank == 0) {
    std::cout << "n,k,M,r,Tp,Tc,f,cl,cm,cu,cv" << std::endl;
    std::cout << global_vertices << ",";
    std::cout << machines << ",";
    std::cout << messages << ",";
    std::cout << rounds << ",";
    std::cout << T0[0] << ",";
    std::cout << T0[1] << ",";
    std::cout << trees << ",";
    std::cout << tree_popl_min << ",";
    std::cout << tree_popl_av1 << ",";
    std::cout << tree_popl_max << ",";
    std::cout << tree_popl_var << std::endl;

    std::cout << "iC,nC" << std::endl;
    for (int tree = 0; tree < trees; tree++) {
      std::cout << forest[tree*2] << "," << forest[tree*2+1] << std::endl;
    }
  }

  // Clean up.
  if (bcast_xinfo) exchange_info_free(bcast_xinfo);
  if (ucast_xinfo) exchange_info_free(ucast_xinfo);

  for (auto &kv : V_in) {
    if (kv.second) free(kv.second);
    kv.second = NULL;
  }
  V_in.clear();

  for (auto &v : V_out) {
    if (v) free(v);
  }
  V_out.clear();

  for (auto &e : E_incoming) {
    e.second.clear();
  }
  E_incoming.clear();

  free(forest);
  // Tear down MPI.
  MPI_Finalize();

  // Exit.
  return 0;
}
