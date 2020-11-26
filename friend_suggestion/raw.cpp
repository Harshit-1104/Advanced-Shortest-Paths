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
    ll dist = INFINIT;
    vector<pair<ll, ll>> adj;
    bool processed = false;
    bool wp = false;
};
/*
class Heap {
  ll maxSize;
  ll size;
  vector<Node> H;
  vector<Node> *graph__;

public:
  Heap(ll maxSize, vector<Node> *graph__) {
    this->size = 0;
    this->maxSize = maxSize;
    this->H.resize(maxSize);
    this->graph__ = graph__;
  } 

  void SiftUp(ll i) {
    while (i > 0 && H[parent(i)].dist > H[i].dist) {
      swap(H[parent(i)], H[i]);
      graph__->at(H[i].id).heap_idx = i;    // Updating the heap_idx of parent
      i = parent(i);
    }

    graph__->at(H[i].id).heap_idx = i;    // Updating the original node id
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
      graph__->at(H[i].id).heap_idx = i;
      SiftDown(maxIdx);
    }
      graph__->at(H[maxIdx].id).heap_idx = maxIdx;   
  }

  void Insert(Node node) {
    maxSize = maxSize +1;

    if (size == maxSize) {
      assert(0);
    }

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

  void ChangePriority(ll id, ll dist) {
    Node old = H[id];
    H[id].dist = dist;

    if (dist < old.dist) {
      SiftUp(id);
    }
    else {
      SiftDown(id);
    }
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
};*/
// Distances can grow out of int type

// Vector of two priority queues - for fwd and bwd searches.
// Each priority queue stores the closest unprocessed node in its head.
// Syntax to create a min heap priority queue.

// STL Queue

struct CompareHeight { 
	bool inline operator()(Node* const& n1, Node* const& n2) 
	{ 
		// return "true" if "p1" is ordered 
		// before "p2", for example: 
		return n1->dist > n2->dist; 
	} 
}; 

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
    vector<Node*> proc_;

    // Stores all the nodes visited either by forward or backward search.
    vector<Node*> workset_;

public:
    Bidijkstra(ll n, vector<Node> graph, vector<Node> graph_rev)
        : n_(n), graph_(graph), graph_rev_(graph_rev)
    { workset_.reserve(2*n);
      proc_.reserve(2*n); 
    } // Preserves this memory from the start so that not much time is wasted everytime while increasing vector size

    // Initialize the data structues before new query,
    // clear the changes made by the previous query.
    void clear() {
        Node* v;
        for (ll i = 0; i < workset_.size(); i++) {
            v = workset_[i];
            graph_[v->id].wp = graph_rev_[v->id].wp = graph_[v->id].processed = graph_rev_[v->id].processed = false;
            graph_[v->id].dist = graph_rev_[v->id].dist = INFINIT;
        }
        workset_.clear();
        proc_.clear();
    }

    // Processes visit of either fwd or bwd search to node v trying to
    // relax the current distance by dist.
    ll Algo(ll s, ll t) {
        clear();
        graph_[s].dist = graph_rev_[t].dist = 0;
        graph_[s].id = s;
        graph_rev_[t].id = t;

	      priority_queue<Node*, vector<Node*>, CompareHeight> fwd, bwd;

        fwd.emplace(&graph_[s]);
        bwd.emplace(&graph_rev_[t]);
/*
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
        }*/

        
        while(true) {
            Node *tmp = fwd.top();
            fwd.pop();

            Process_fwd(tmp, &fwd);

            if (tmp->wp == false) {
              tmp->wp = true;
              workset_.emplace_back(tmp);
            }

            if (graph_rev_[tmp->id].processed == true) {
                return ShortestPath(s, t); 
            }

            tmp = bwd.top();
            bwd.pop();
            Process_bwd(tmp, &bwd);
            
            if (tmp->wp == false) {
              tmp->wp = true;
              workset_.emplace_back(tmp);
            }

            if (graph_[tmp->id].processed == true) {
                return ShortestPath(s, t);
            }

            if (fwd.empty() || bwd.empty()) {
              return -1;
            }
        }
    }
    
    void Process_fwd(Node *u, priority_queue<Node*, vector<Node*>, CompareHeight> *fwd) {

        for (ll i = 0; i < u->adj.size(); i++) {
            pair<ll, ll> *tmp = &u->adj[i];
            if (graph_[tmp->first].dist > u->dist + tmp->second) {
                graph_[tmp->first].dist = u->dist + tmp->second;
                
                graph_[tmp->first].id = tmp->first;

                if (graph_[tmp->first].wp == false) {
                  graph_[tmp->first].wp = true;
                  workset_.emplace_back(&graph_[tmp->first]);
                }

                //graph_[tmp.first].prev = u->id;
                fwd->emplace(&graph_[tmp->first]);
            }
        }
        graph_[u->id].processed = true;
        proc_.emplace_back(&graph_[u->id]);
    }

    void Process_bwd(Node *u, priority_queue<Node*, vector<Node*>, CompareHeight> *bwd) {
        for (ll i = 0; i < u->adj.size(); i++) {
            pair<ll, ll>* tmp = &u->adj[i];
            if (graph_rev_[tmp->first].dist > u->dist + tmp->second) {
                graph_rev_[tmp->first].dist = u->dist + tmp->second;
                
                graph_rev_[tmp->first].id = tmp->first;

                if (graph_rev_[tmp->first].wp == false) {
                  graph_rev_[tmp->first].wp = true;
                  workset_.emplace_back(&graph_rev_[tmp->first]);
                }

                //graph_rev_[tmp.first].prev = u->id;
                bwd->emplace(&graph_rev_[tmp->first]);
            }
        }
        graph_rev_[u->id].processed = true;
        proc_.emplace_back(&graph_rev_[u->id]);
    }

    ll ShortestPath(ll s, ll t) {
        ll distance = INFINIT;
        Node* u_best;

        for (ll i = 0; i < proc_.size(); i++) {
            Node* u = proc_[i];

            if (graph_[u->id].dist + graph_rev_[u->id].dist < distance) {
                u_best = u;
                distance = graph_rev_[u->id].dist + graph_[u->id].dist;
            }
        }

        if (distance >= INFINIT) {
            return -1;
        }
        return distance;
    }

    

/*    
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
                fwd.ChangePriority(graph_[tmp.first].heap_idx, graph_[tmp.first].dist);
            }
        }
        graph_[u->id].processed = true;
        proc_fwd_.push_back(graph_[u->id]);
    }

    void Process_bwd(Node *u, Heap &bwd) {
        for (ll i = 0; i < u->adj.size(); i++) {
            pair<ll, ll> tmp = u->adj[i];
            if (graph_rev_[tmp.first].dist > u->dist + tmp.second) {
                graph_rev_[tmp.first].dist = u->dist + tmp.second;
                
                if (graph_rev_[tmp.first].wp == false) {
                  graph_rev_[tmp.first].wp = true;
                  workset_.push_back(graph_rev_[tmp.first]);
                }

                graph_rev_[tmp.first].prev = u->id;
                bwd.ChangePriority(graph_rev_[tmp.first].heap_idx, graph_rev_[tmp.first].dist);
            }
        }
        graph_rev_[u->id].processed = true;
        proc_bwd_.push_back(graph_rev_[u->id]);
    }
*/
};

/*
ll dijkstra(vector<Node> graph, int s, int t) {
  graph[s].dist = 0; // Source
  priority_queue<Node, vector<Node>, CompareHeight> heap; // The keys will be distance values
  heap.push(graph[s]);

  if (t == s) {
    return 0;
  }
  
  while (!heap.empty()) {
    int u = (heap.top()).id;
    heap.pop();

    if (graph[u].processed == true) {
      continue;
    }

    for (int i = 0; i < graph[u].adj.size(); i++) {
      if (graph[graph[u].adj[i].first].dist > graph[u].dist + graph[u].adj[i].second) {
        graph[graph[u].adj[i].first].dist = graph[u].dist + graph[u].adj[i].second;
        graph[graph[u].adj[i].first].prev = graph[u].id;
        heap.push(graph[graph[u].adj[i].first]);
      }
    }
  }
  
  return (graph[t].dist != INFINIT) ? (graph[t].dist) : (-1);
}
*/
    //graph[u].processed = true;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);

    ll n, m;
    scanf("%lld%lld", &n, &m);
    vector<Node> graph(n), graph_rev(n);

    for (int i = 0; i < m; i++) {
        int u, v, c;
        scanf("%d%d%d", &u, &v, &c);
        graph[u-1].adj.emplace_back(make_pair(v-1, c));
        graph_rev[v-1].adj.emplace_back(make_pair(u-1, c));
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
  /*    
    for (int i = 1; i <= n; ++i) {
      for (int j = i; j <= n; j++) {
        ll u = i, v = j;
        //scanf("%d%d", &u, &v);
        ll dist1 = dijkstra(graph, u-1, v-1);
        ll dist2 = bidij.Algo(u-1, v-1);
        printf("%lld : %lld -> %lld \t", u-1, v-1, dist2);
        if (dist1 != dist2) {
          exit(0);
        }

        dist1 = dijkstra(graph, v-1, u-1);
        dist2 = bidij.Algo(v-1, u-1);
        printf("%lld : %lld -> %lld \t", v-1, u-1, dist2);
        if (dist1 != dist2) {
          exit(0);
        }
      }
    }*/
}