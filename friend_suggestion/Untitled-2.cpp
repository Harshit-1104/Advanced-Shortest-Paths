#include <cstdio>
#include <cassert>
#include <vector>
#include <queue>
#include <limits>
#include <utility>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <iostream>
#include <algorithm>

using namespace std;
typedef long long ll;

const ll INFINIT = numeric_limits<ll>::max() / 4;

struct Node {
    ll id;
    ll dist;
    vector<pair<ll, ll>> adj;
    bool processed;
    ll prev;

    bool wp = false;
};

class Heap {
  ll maxSize;
  ll size;
  vector<Node> H;

public:
  Heap(ll maxSize) {
    this->size = 0;
    this->maxSize = maxSize;
    this->H.resize(maxSize);
  } 

  void SiftUp(ll i) {
    //int dist = H[parent(i)].dist;  
    while (i > 0 && H[parent(i)].dist > H[i].dist) {
      swap(H[parent(i)], H[i]);
      i = parent(i);
    }
  }

  void SiftDown(ll i) {
    ll maxIdx = i;
    ll left = leftChild(i), right = rightChild(i);

    if (left <= size-1 && H[left].dist < H[maxIdx].dist) {
      maxIdx = left;
    }
    
    if (right <= size-1 && H[right].dist < H[maxIdx].dist) {
      maxIdx = right;
    }

    if (i != maxIdx) {
      swap(H[i], H[maxIdx]);
      SiftDown(maxIdx);
    }
  }

  void Insert(Node node) {
    maxSize = maxSize +1;

/*    if (size == maxSize) {
      exit(0);
    }*/

    size = size+1;
    H[size-1] = node;
    SiftUp(size-1);
  }

  Node ExtractMin() {
    Node node = H[0];
    H[0] = H[size -1];
    size = size-1;
    SiftDown(0);
    return node;
  }

  bool empty() {
    if (size == 0)
      return true;
    else
      return false;
  }

  ll parent(ll i) {
    return floor((i-1)/2);
  }

  ll leftChild(ll i) {
    return (2*i +1);
  }

  ll rightChild(ll i) {
    return (2*i +2);
  }
};

// Distances can grow out of int type

// Vector of two priority queues - for fwd and bwd searches.
// Each priority queue stores the closest unprocessed node in its head.
// Syntax to create a min heap priority queue.

class Bidijkstra {
    // Number of nodes
    ll n_;
    // Graph adj_[0] and cost_[0] correspond to the intital graph,
    // adj_[1] and cost_[1] correspond to the reversed graph.
    // Graphs are stored as vectors of adjacency lists corresponding
    // to nodes.
    // Adjacency list itself is stored in adj_, and the corresponding 
    // edge costs are stored in cost_.
    vector<Node> graph_, graph_rev_;

    // distance_[0] stores distances for the forward search,
    // and distance_[1] stores distances for the backward search.
    vector<Node> proc_fwd_, proc_bwd_;

    // Stores all the nodes visited either by forward or backward search.
    vector<Node> workset_;

public:
    Bidijkstra(ll n, vector<Node> graph, vector<Node> graph_rev)
        : n_(n), graph_(graph), graph_rev_(graph_rev)
    { workset_.reserve(n); } // Preserves this memory from the start so that not much time is wasted everytime while increasing vector size

    // Initialize the data structues before new query,
    // clear the changes made by the previous query.
    void clear() {
        for (ll i = 0; i < workset_.size(); i++) {
            Node v = workset_[i];
            graph_[v.id].wp = graph_rev_[v.id].wp = false;
            graph_[v.id].dist = graph_rev_[v.id].dist = INFINIT;
            graph_[v.id].processed = graph_rev_[v.id].processed = false;
            graph_[v.id].prev = graph_rev_[v.id].prev = -1; 
        }
        workset_.clear();
    }

    // Processes visit of either fwd or bwd search to node v trying to
    // relax the current distance by dist.
    ll Algo(ll s, ll t) {
        clear();
        graph_[s].dist = graph_rev_[t].dist = 0;
        Heap fwd(graph_.size()*1000), bwd(graph_rev_.size()*1000);
        for (ll i = 0 ; i < n_; i++) {
            fwd.Insert(graph_[i]);
            bwd.Insert(graph_rev_[i]);
        }

        while(true) {
            Node tmp = fwd.ExtractMin();
            Process_fwd(&tmp, fwd);

            Node *ptr = &graph_[tmp.id];
            if (ptr->wp == false) {
              ptr->wp = true;
              workset_.push_back(*ptr);
            }

            if (graph_rev_[tmp.id].processed == true) {
                return ShortestPath(s, t); 
            }

            tmp = bwd.ExtractMin();
            Process_bwd(&tmp, bwd);
            
            Node *ptr_ = &graph_rev_[tmp.id];
            if (ptr_->wp == false) {
              ptr_->wp = true;
              workset_.push_back(*ptr_);
            }

            if (graph_[tmp.id].processed == true) {
                return ShortestPath(s, t);
            }
        }
    }

    ll ShortestPath(ll s, ll t) {
        ll distance = INFINIT;
        Node u_best;

        for (ll i = 0; i < proc_fwd_.size() + proc_bwd_.size() +2; i++) {
            Node u;

            if (!proc_fwd_.empty()) {
                u = proc_fwd_.back();
                proc_fwd_.pop_back();
            } else {
                u = proc_bwd_.back();
                proc_bwd_.pop_back();
            }

            if (graph_[u.id].dist + graph_rev_[u.id].dist < distance) {
                u_best = u;
                distance = graph_rev_[u.id].dist + graph_[u.id].dist;
            }
        }

        if (distance >= INFINIT) {
            return -1;
        }
        return distance;
    }

    void Process_fwd(Node *u, Heap &fwd) {

        for (ll i = 0; i < u->adj.size(); i++) {
            pair<ll, ll> tmp = u->adj[i];
            if (graph_[tmp.first].dist > u->dist + tmp.second) {
                graph_[tmp.first].dist = u->dist + tmp.second;
                
                if (graph_[tmp.first].wp == false) {
                  graph_[tmp.first].wp = true;
                  workset_.push_back(graph_[tmp.first]);
                }

                graph_[tmp.first].prev = u->id;
                fwd.Insert(graph_[tmp.first]);
            }
        }
        graph_[u->id].processed = true;
        proc_fwd_.push_back(graph_[u->id]);
    }

    void Process_bwd(Node *u, Heap &bwd) {
        for (ll i = 0; i < u->adj.size(); i++) {
            pair<ll, ll> tmp = u->adj[i];
            if (graph_rev_[tmp.first].dist > u->dist + u->tmp.second) {
                graph_rev_[tmp.first].dist = u->dist + tmp.second;
                
                if (graph_rev_[tmp.first].wp == false) {
                  graph_rev_[tmp.first].wp = true;
                  workset_.push_back(graph_rev_[tmp.first]);
                }

                graph_rev_[tmp.first].prev = u->id;
                bwd.Insert(graph_rev_[tmp.first]);
            }
        }
        graph_rev_[u->id].processed = true;
        proc_bwd_.push_back(graph_rev_[u->id]);
    }
};

int main() {
    ll n, m;
    scanf("%lld%lld", &n, &m);
    vector<Node> graph(n), graph_rev(n);

    for (int i = 0; i < n; i++) {
        graph[i].id = graph_rev[i].id = i;
        graph[i].dist = graph_rev[i].dist = INFINIT;
        graph[i].processed = graph_rev[i].processed = false;
        graph[i].prev = graph_rev[i].prev = -1;
    }

    for (int i = 0; i < m; i++) {
        int u, v, c;
        scanf("%d%d%d", &u, &v, &c);
        graph[u-1].adj.push_back(make_pair(v-1, c));
        graph_rev[v-1].adj.push_back(make_pair(u-1, c));
    }

    Bidijkstra bidij(n, graph, graph_rev);

    int t;
    scanf("%d", &t);
 
    for (int i = 0; i < t; i++) {
      ll u, v;
      scanf("%lld%lld", &u, &v);
      ll dist1 = bidij.Algo(u-1, v-1);
      printf("%lld\n", dist1);
    }
      // Testing zone
}