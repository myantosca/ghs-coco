#include <cstdio>
#include <cstdint>
#include <cstring>
#include <sys/time.h>
#include <random>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <chrono>
#include <mpi.h>
#include <queue>
#include <cassert>
using namespace std::chrono;

// https://en.wikipedia.org/wiki/Mersenne_prime#List_of_known_Mersenne_primes
#define MERSENNE_61 2305843009213693951
#define MACHINE_HASH(__VERTEX__) ((hash_a * __VERTEX__ + hash_b) % MERSENNE_61) % machines

typedef enum {
  IDLE = 0,
  FIND_SEND = 1,
  FIND_WAIT = 2,
  FIND_TEST = 3,
  FIND_RPLY = 4,
  MWOE_SEND = 5
} vertex_state_t;

typedef enum {
  PING = 0,
  PONG = 1,
  FIND = 2,
  FOUND = 3,
  MWOE = 4,
  JOIN = 5,
  NEW = 6
} msg_type_t;

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


typedef struct edge
{
  uint32_t u;
  uint32_t v;
} edge_t;

typedef struct vertex
{
  uint32_t id;
  uint32_t parent;
  uint32_t group;
  uint32_t group_ct;
  vertex_state_t state;
  uint32_t awaiting;
  edge_t mwoe;
  std::set<uint32_t> inactive_neighbors;
  std::set<uint32_t> active_neighbors;
  std::unordered_set<uint32_t> children;
} vertex_t;

typedef struct quad_msg
{
  uint32_t typ;
  uint32_t dst;
  uint32_t a;
  uint32_t b;
} quad_msg_t;

typedef struct census_msg {
  uint32_t id;
  uint32_t ct;
} census_msg_t;

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

void exchange_one(exchange_info_t *info, uint64_t *messages, int machine) {
  int received = 0;
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

void exchange_all(exchange_info_t *info, uint64_t *messages) {
  memset(info->recv_counts, 0, sizeof(int) * info->machines);
  memset(info->recv_displs, 0, sizeof(int) * info->machines);
  // Gather for each machine in turn from every other machine's
  // targeted send buffers.
  for (int machine = 0; machine < info->machines; machine++) {
    exchange_one(info, messages, machine);
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
  exchange_info_t *ghs_xinfo = exchange_info_new(rank, machines);

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
  exchange_all(edges_xinfo, &messages);

  // Construct local vertex-centric model.
  std::map<uint32_t, vertex_t *> V_r;
  std::map<uint32_t, std::unordered_set<uint32_t>> E_incoming;
  for (int i = 0; i < edges_xinfo->recv_caps; i += 2) {
    uint32_t id_u = edges_xinfo->recv_buf[i];
    uint32_t id_v = edges_xinfo->recv_buf[i+1];
    vertex_t *u = V_r[id_u];
    // If the node has not been created yet, do so.
    if (!u) {
      u = (vertex_t *)malloc(sizeof(vertex_t));
      memset(u, 0, sizeof(vertex_t));
      u->id = id_u;
      u->parent = id_u;
      u->group = id_u;
      u->group_ct = 1; // population of GHS sub-tree including itself
      u->state = IDLE;
      u->awaiting = 0;
      u->active_neighbors = std::set<uint32_t>();
      u->inactive_neighbors = std::set<uint32_t>();
      u->children = std::unordered_set<uint32_t>();
      V_r[id_u] = u;
    }
    if (id_u != id_v) {
      u->active_neighbors.insert(id_v);
      if (E_incoming.find(id_v) == E_incoming.end()) { E_incoming[id_v] = std::unordered_set<uint32_t>(); }
      E_incoming[id_v].insert(id_u);
    }
  }

  // Clean up so as not hoard memory.
  exchange_info_free(edges_xinfo);

  // Get local stats.
  local_vertices = V_r.size();

  // Capture finishing time for preprocessing.
  time_point<system_clock, nanoseconds> tp_prep_b = system_clock::now();

  if (verbosity == 1) {
    // Debug printout of graph in vertex-centric format.
    for (int machine = 0; machine < machines; machine++) {
      if (rank == machine) {
        for (auto &kv : V_r) {
          std::cerr << "[" << machine << "]";
          std::cerr << kv.first << ":";
          for (auto &neighbor : kv.second->active_neighbors) {
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
  std::map<uint32_t, vertex_t *> T_r;
  std::map<uint32_t, vertex_t *> T_0;
  std::queue<vertex_t *> S_r;
  // Separate singletons from connected vertices.
  for (auto &kv : V_r) {
    if (kv.second->inactive_neighbors.empty() && kv.second->active_neighbors.empty()) {
      T_0[kv.first] = kv.second;
    }
    else {
      T_r[kv.first] = kv.second;
    }
  }

  // Initialize termination conditions.
  int global_done = 0;
  int local_done = (T_r.size() == 0);

  MPI_Allreduce(&local_done, &global_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
  messages += machines << 1;
  // Continue until all edges have been rendered inactive.

  //std::queue<vertex_t *> M_r;

  while (!global_done) {
    // Find phase

    // All roots of connected components must initiate the MWOE search.
    for (auto &kv : T_r) {
      // Vertices find local MWOE and then query children to avoid code duplication.
      kv.second->state = FIND_TEST;
      // Defaults for MWOE = link to self (special meaning, i.e., no outgoing edge)
      kv.second->mwoe.u = kv.second->id;
      kv.second->mwoe.v = kv.second->id;
      S_r.push(kv.second);
    }

    uint32_t w_L = V_r.size() - T_0.size();
    int local_link_done = (w_L == 0);
    int global_link_done = 0;
    MPI_Allreduce(&local_link_done, &global_link_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    while (!global_link_done) {
      int m_v;
      // Link send
      std::cerr << "------ SEND ------" << std::endl;
      while (!S_r.empty()) {
	quad_msg_t req;
	vertex_t *u = S_r.front();
	S_r.pop();
	if (u->state == FIND_SEND) {
	  u->awaiting = u->children.size();
	  if (u->awaiting > 0) {
	    u->state = FIND_WAIT;
	    for (auto &id_v : u->children) {
	      req.typ = FIND;
	      req.dst = id_v;
	      req.a = u->id;
	      req.b = u->group;
	      m_v = MACHINE_HASH(req.dst);
	      exchange_info_send_buf_insert(ghs_xinfo, m_v, (uint32_t *)&req, 4);
	      std::cerr << "FIND " << req.dst << " " << req.a << " " << req.b << std::endl;
	    }
	  }
	  else {
	    u->state = FIND_RPLY;
	    S_r.push(u);
	  }
	}
	else if (u->state == FIND_TEST) {
	  if (u->active_neighbors.size() > 0) {
	    // Edges left to check
	    req.typ = PING;
	    req.dst = *u->active_neighbors.begin();
	    req.a = u->id;
	    req.b = u->group;
	    m_v = MACHINE_HASH(req.dst);
	    exchange_info_send_buf_insert(ghs_xinfo, m_v, (uint32_t *)&req, 4);
	    std::cerr << "PING " << req.dst << " " << req.a << " " << req.b << std::endl;
	  }
	  else {
	    // No more incident edges to check. MWOE is what it is.
	    u->state = FIND_SEND;
	    S_r.push(u);
	  }
	}
	else if (u->state == FIND_RPLY) {
	  if (u->parent == u->id) {
	    // At the root. Time to announce the MWOE.
	    u->state = MWOE_SEND;
	    S_r.push(u);
	  }
	  else {
	    u->state = IDLE;
	    req.typ = FOUND;
	    req.dst = u->parent;
	    req.a = u->mwoe.u;
	    req.b = u->mwoe.v;
	    m_v = MACHINE_HASH(req.dst);
	    exchange_info_send_buf_insert(ghs_xinfo, m_v, (uint32_t *)&req, 4);
	    std::cerr << "FOUND " << req.dst << " " << req.a << " " << req.b << std::endl;
	  }
	}
	else if (u->state == MWOE_SEND) {
	  u->state = IDLE;
	  w_L--;
	  for (auto &id_v : u->children) {
	    req.typ = MWOE;
	    req.dst = id_v;
	    req.a = u->mwoe.u;
	    req.b = u->mwoe.v;
	    m_v = MACHINE_HASH(req.dst);
	    exchange_info_send_buf_insert(ghs_xinfo, m_v, (uint32_t *)&req, 4);
	    std::cerr << "MWOE " << req.dst << " " << req.a << " " << req.b << std::endl;
	  }
	}
      }
      // Exchange messages between all machines.
      exchange_all(ghs_xinfo, &messages);
      // Reset, rewind.
      exchange_info_rewind(ghs_xinfo);

      std::cerr << "------ RECV ------" << std::endl;
      // Link receive
      for (int i = 0; i < ghs_xinfo->recv_caps >> 2; i++) {
	quad_msg_t req = ((quad_msg_t *)ghs_xinfo->recv_buf)[i];
	quad_msg_t rsp;
	vertex_t *v = V_r[req.dst];

	if (req.typ == FIND) {
	  assert(v->state == IDLE);
	  v->awaiting = v->children.size();
	  std::cerr << "FIND " << req.dst << " " << req.a << " " << req.b << std::endl;
	  v->state = FIND_SEND;
	  S_r.push(v);
	}
	else if (req.typ == PING) {
	  std::cerr << "PING " << req.dst << " " << req.a << " " << req.b << std::endl;
	  rsp.typ = PONG;
	  rsp.dst = req.a;
	  rsp.a = v->id;
	  rsp.b = v->group;
	  int m_u = MACHINE_HASH(rsp.dst);
	  // No state change. Just respond.
	  exchange_info_send_buf_insert(ghs_xinfo, m_u, (uint32_t *)&rsp, 4);
	}
	else if (req.typ == PONG) {
	  assert(v->state == FIND_TEST);
	  std::cerr << "PONG " << req.dst << " " << req.a << " " << req.b << std::endl;
	  if (req.b == v->group) {
	    // Keep testing.
	    v->active_neighbors.erase(req.a);
	    v->inactive_neighbors.insert(req.a);
	  }
	  else {
	    // Found a minimum edge.
	    v->mwoe.u = v->id;
	    v->mwoe.v = req.a;
	    v->state = FIND_SEND;
	  }
	  // Regardless of the PONG, take action next cycle.
	  S_r.push(v);
	}
	else if (req.typ == FOUND) {
	  assert(v->state == FIND_WAIT);
	  std::cerr << "FOUND " << req.dst << " " << req.a << " " << req.b << "(" << v->awaiting << ")" << std::endl;
	  v->awaiting--;
	  if (req.a != req.b) {
	    if (v->mwoe.u != v->mwoe.v) {
	      edge_t e_msg = (req.a < req.b
			      ? (edge_t){ req.a, req.b}
			      : (edge_t){ req.b, req.a });
	      edge_t e_vtx = (v->mwoe.u < v->mwoe.v
			      ? (edge_t){ v->mwoe.u, v->mwoe.v }
			      : (edge_t){ v->mwoe.v, v->mwoe.u });
	      if ((e_msg.u < e_vtx.u) ||
		  ((e_msg.u == e_vtx.u) && (e_msg.v < e_vtx.v))) {
		v->mwoe.u = req.a;
		v->mwoe.v = req.b;
	      }
	    }
	    else {
	      v->mwoe.u = req.a;
	      v->mwoe.v = req.b;
	    }
	  }
	  if (v->awaiting == 0) {
	    v->state = FIND_RPLY;
	    S_r.push(v);
	  }
	}
	else if (req.typ == MWOE) {
	  assert(v->state == IDLE);
	  std::cerr << "MWOE " << req.dst << " " << req.a << " " << req.b << std::endl;

	  if (v->id == req.a) {
	    if (v->id != v->group) {
	      // Fragment re-roots to vertex with minimum edge.
	      T_r[v->id] = v;
	      T_r.erase(v->group);
	    }
	  }
	  v->mwoe.u = req.a;
	  v->mwoe.v = req.b;
	  v->state = MWOE_SEND;

	  S_r.push(v);
	}
	else {
	  std::cerr << req.typ << " " << req.dst << " " << req.a << " " << req.b << std::endl;
	}
      }
      // Link phase termination condition.
      local_link_done = (w_L == 0);
      MPI_Allreduce(&local_link_done, &global_link_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    }

    // Merge Phase
    std::cerr << "****** MERGE *******" << std::endl;
    uint32_t joins = 0;
    // Join send
    for (auto &kv : T_r) {
      vertex_t *u = kv.second;
      if (u->active_neighbors.size() > 0) {
	quad_msg_t msg;
	msg.typ = JOIN;
	msg.dst = u->mwoe.v;
	msg.a = u->id;
	msg.b = u->group;
	int m_v = MACHINE_HASH(msg.dst);
	exchange_info_send_buf_insert(ghs_xinfo, m_v, (uint32_t *)&msg, 4);
	joins++;
      }
    }

    if (joins) {
      exchange_all(ghs_xinfo, &messages);
      exchange_info_rewind(ghs_xinfo);
      // Join receive
      for (int i = 0; i < ghs_xinfo->recv_caps >> 2; i++) {
	quad_msg_t req = ((quad_msg_t *)ghs_xinfo->recv_buf)[i];
	vertex_t *v = V_r[req.dst];
	std::cerr << "JOIN " << req.dst << " " << req.a << " " << req.b << std::endl;
	if (req.a == v->mwoe.v) {
	  if (v->id > req.a) {
	    if (v->parent != v->id) { v->children.insert(v->parent); }
	    std::cerr << "***** " << v->id << " *****" << std::endl;
	    v->parent = v->id;
	    v->group = v->id;
	    S_r.push(v);
	    v->children.insert(req.a);
	  }
	}
	else {
	  v->children.insert(req.a);
	  v->inactive_neighbors.insert(req.a);
	  v->active_neighbors.erase(req.a);
	}
      }
    }

    int local_merge_done = S_r.empty();
    int global_merge_done = 0;
    MPI_Allreduce(&local_merge_done, &global_merge_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    while (!global_merge_done) {
      while (!S_r.empty()) {
	vertex_t *u = S_r.front();
	S_r.pop();
	for (auto &id_v : u->children) {
	  int m_v = MACHINE_HASH(id_v);
	  quad_msg_t msg;
	  msg.typ = NEW;
	  msg.dst = id_v;
	  msg.a = u->id;
	  msg.b = u->group;
	  exchange_info_send_buf_insert(ghs_xinfo, m_v, (uint32_t *)&msg, 4);
	}
      }
      exchange_all(ghs_xinfo, &messages);
      exchange_info_rewind(ghs_xinfo);
      for (int i = 0; i < ghs_xinfo->recv_caps >> 2; i++) {
	quad_msg_t req = ((quad_msg_t *)ghs_xinfo->recv_buf)[i];
	vertex_t *v = V_r[req.dst];
	std::cerr << "NEW " << req.dst << " " << req.a << " " << req.b << std::endl;
	v->group = req.b;
	if ((v->parent != v->id) && (v->parent != req.a)) { v->children.insert(v->parent); }
	v->parent = req.a;
	v->children.erase(req.a);
	v->inactive_neighbors.insert(req.a);
	v->active_neighbors.erase(req.a);
	// Defaults for MWOE = link to self (special meaning, i.e., no outgoing edge)
	v->mwoe.u = v->id;
	v->mwoe.v = v->id;
	T_r.erase(v->id);
	S_r.push(v);
      }
      local_merge_done = S_r.empty();
      MPI_Allreduce(&local_merge_done, &global_merge_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    }

    local_done = (joins == 0);
    MPI_Allreduce(&local_done, &global_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

    messages += machines << 1;
  }

  std::cerr << "merging" << std::endl;
  // Merge all local component lists into global component list.
  for (auto &kv : T_r) {
    kv.second->awaiting = kv.second->children.size();
    S_r.push(kv.second);
  }

  uint32_t w_G = S_r.size();
  int local_census_done = (w_G == 0);
  int global_census_done = 0;
  MPI_Allreduce(&local_census_done, &global_census_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
  std::cerr << "CENSUS" << std::endl;
  while (!global_census_done) {
    while(! S_r.empty()) {
      vertex_t *u = S_r.front();
      S_r.pop();

      if (u->awaiting == 0) {
	if (u->group == u->id) {
	  w_G--;
	}
	else {
	  quad_msg_t msg;
	  msg.typ = PONG;
	  msg.dst = u->parent;
	  msg.a = u->id;
	  msg.b = u->group_ct;
	  int m_v = MACHINE_HASH(u->parent);
	  exchange_info_send_buf_insert(ghs_xinfo, m_v, (uint32_t *)&msg, 4);
	}
      }
      else {
	for (auto &id_v : u->children) {
	  quad_msg_t msg;
	  msg.typ = PING;
	  msg.dst = id_v;
	  msg.a = u->id;
	  msg.b = 0;
	  int m_v = MACHINE_HASH(id_v);
	  exchange_info_send_buf_insert(ghs_xinfo, m_v, (uint32_t *)&msg, 4);
	}
      }
    }

    exchange_all(ghs_xinfo, &messages);
    exchange_info_rewind(ghs_xinfo);
    for (int i = 0; i < ghs_xinfo->recv_caps >> 2; i++) {
      quad_msg_t req = ((quad_msg_t *)ghs_xinfo->recv_buf)[i];
      vertex_t *v = V_r[req.dst];
      if (req.typ == PING) {
	v->awaiting = v->children.size();
	S_r.push(v);
      }
      else if (req.typ == PONG) {
	v->awaiting--;
	v->group_ct += req.b;
	if (v->awaiting == 0) { S_r.push(v); }
      }
    }
    local_census_done = (w_G == 0);
    MPI_Allreduce(&local_census_done, &global_census_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
  }

  for (auto &kv : T_r) {
    vertex_t *u = kv.second;
    census_msg_t msg;
    msg.id = u->id;
    msg.ct = u->group_ct;
    exchange_info_send_buf_insert(ghs_xinfo, 0, (uint32_t *)&msg, 2);
  }
  for (auto &kv : T_0) {
    vertex_t *u = kv.second;
    census_msg_t msg;
    msg.id = u->id;
    msg.ct = u->group_ct;
    exchange_info_send_buf_insert(ghs_xinfo, 0, (uint32_t *)&msg, 2);
  }
  exchange_one(ghs_xinfo, &messages, 0);

  // Capture finishing time for connected component search.
  time_point<system_clock, nanoseconds> tp_coco_b = system_clock::now();
  int64_t T[2] = { (tp_prep_b - tp_prep_a).count(), (tp_coco_b - tp_coco_a).count() };
  int64_t T0[2] = { 0, 0 };
  MPI_Reduce(&T, &T0, 2, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_vertices, &global_vertices, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  // Print out connected component labels with component population counts.
  if (rank == 0) {
    trees = ghs_xinfo->recv_caps >> 1;
    forest = ghs_xinfo->recv_buf;
    for (int tree = 0; tree < trees; tree++) {
      uint32_t tree_popl = forest[tree*2+1];
      if (tree_popl < tree_popl_min) { tree_popl_min = tree_popl; }
      if (tree_popl > tree_popl_max) { tree_popl_max = tree_popl; }
      tree_popl_av0 = tree_popl_av1;
      tree_popl_av1 += (tree_popl - tree_popl_av0) / (tree + 1);
      tree_popl_ess += (tree_popl - tree_popl_av0) * (tree_popl - tree_popl_av1);
      tree_popl_var = tree_popl_ess / (tree + 1);
    }
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
  if (ghs_xinfo) exchange_info_free(ghs_xinfo);

  for (auto &kv : V_r) {
    if (kv.second) free(kv.second);
    kv.second = NULL;
  }
  V_r.clear();

  T_r.clear();

  for (auto &e : E_incoming) {
    e.second.clear();
  }
  E_incoming.clear();

  // Tear down MPI.
  MPI_Finalize();

  // Exit.
  return 0;
}
