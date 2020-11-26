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
//#include <Windows.h>
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
            workset_.push_back(graph_[tmp.id]);

            if (graph_rev_[tmp.id].processed == true) {
                return ShortestPath(s, t); 
            }

            tmp = bwd.ExtractMin();
            Process_bwd(&tmp, bwd);
            workset_.push_back(graph_rev_[tmp.id]);

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
            cout << "-1" << endl;
            return -1;
        }

        cout << distance << " :: ";

        vector<int> Path;
        Node last = u_best;

        while (last.id != graph_[s].id && last.prev != -1) {
          Path.push_back(last.id);
          last = graph_[last.prev];
        }
        Path.push_back(last.id);

        reverse(begin(Path), end(Path));

        last = u_best;
        while (last.id != graph_rev_[t].id && graph_rev_[last.id].prev != -1) {
          last = graph_rev_[graph_rev_[last.id].prev];
          Path.push_back(last.id);
        }

        for (int i = 0; i < Path.size(); i++) {
          cout << Path[i] << " ";
        }
        cout <<"\n";

        return distance;
    }

    void Process_fwd(Node *u, Heap &fwd) {
        int size = u->adj.size();
        for (ll i = 0; i < u->adj.size(); i++) {
            pair<ll, ll> tmp = u->adj[i];
            ll dist = graph_[u->adj[i].first].dist;
            if (graph_[u->adj[i].first].dist > u->dist + u->adj[i].second) {
                graph_[u->adj[i].first].dist = u->dist + u->adj[i].second;
                
                workset_.push_back(graph_[u->adj[i].first]);
                graph_[u->adj[i].first].prev = u->id;
                fwd.Insert(graph_[u->adj[i].first]);
            }
        }
        graph_[u->id].processed = true;
        proc_fwd_.push_back(graph_[u->id]);
    }

    void Process_bwd(Node *u, Heap &bwd) {
        ll size = u->adj.size();
        for (ll i = 0; i < u->adj.size(); i++) {
            pair<ll, ll> tmp = u->adj[i];
            ll dist = graph_rev_[u->adj[i].first].dist;
            if (graph_rev_[u->adj[i].first].dist > u->dist + u->adj[i].second) {
                graph_rev_[u->adj[i].first].dist = u->dist + u->adj[i].second;
                
                workset_.push_back(graph_rev_[u->adj[i].first]);
                graph_rev_[u->adj[i].first].prev = u->id;
                bwd.Insert(graph_rev_[u->adj[i].first]);
            }
        }
        graph_rev_[u->id].processed = true;
        proc_bwd_.push_back(graph_rev_[u->id]);
    }
};


ll dijkstra(vector<Node> graph, int s, int t) {
  graph[s].dist = 0; // Source
  Heap heap(graph.size() * 1000); // The keys will be distance values
  heap.Insert(graph[s]);
  
  while (!heap.empty()) {
    int u = (heap.ExtractMin()).id;
    if (graph[u].processed == true) {
      continue;
    }

    for (int i = 0; i < graph[u].adj.size(); i++) {
      if (graph[graph[u].adj[i].first].dist > graph[u].dist + graph[u].adj[i].second) {
        graph[graph[u].adj[i].first].dist = graph[u].dist + graph[u].adj[i].second;
        graph[graph[u].adj[i].first].prev = graph[u].id;
        heap.Insert(graph[graph[u].adj[i].first]);
      }
    }

    //graph[u].processed = true;
  }

  return (graph[t].dist != INFINIT) ? (graph[t].dist) : (-1);
}

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
        graph[u].adj.push_back(make_pair(v, c));
       graph_rev[v].adj.push_back(make_pair(u, c));
    }

    Bidijkstra bidij(n, graph, graph_rev);

    ll t;
//    scanf("%d", &t);
    for (int i = 1; i <= n; ++i) {
      for (int j = i; j <= n; j++) {
        ll u = i, v = j;
        //scanf("%d%d", &u, &v);
        ll dist1 = dijkstra(graph, u-1, v-1);
        printf("%lld : %lld -> %lld \t", u-1, v-1, dist1);
        ll dist2 = bidij.Algo(u-1, v-1);
        if (dist1 != dist2) {
          exit(0);
        }

        dist1 = dijkstra(graph, v-1, u-1);
        printf("%lld : %lld -> %lld \t", v-1, u-1, dist1);
        dist2 = bidij.Algo(v-1, u-1);
        if (dist1 != dist2) {
          exit(0);
        }
      }
    }

      // Testing zone
/*
  ll later = 0;
  while (true) {
    Sleep(1000);
    srand(time(0));
    ll n = rand() % 20 + 2;              // Number of vertices
    if (n <= 3) {
      n += later;
    }

    Sleep(1000);
    srand(time(0));
    ll m = (n * (n-1));       // Number of edges
    std::cout << n << " " << m << "\n";
    vector<Node> graph(n), graph_rev(n);

    for (ll i = 0; i < n; i++) {
        graph[i].id = graph_rev[i].id = i;
        graph[i].dist = graph_rev[i].dist = INFINIT;
        graph[i].processed = graph_rev[i].processed = false;
    }
    int ctr2 = 0;

    for (ll i = 0; i < n; i++) {
      for (ll j = i+1; j < n; j++) {
        if (ctr2 == 60) { 
          Sleep(1000);
          srand(time(0));
        }
        ll c = (rand() % 100 + 50) % 20 + rand() % 10;
        graph[i].adj.push_back(make_pair(j, c));
        graph_rev[j].adj.push_back(make_pair(i, c));

        cout << i << " " << j << " " << c << "\n";
        
        if (ctr2 == 60) { 
          Sleep(1000);
          srand(time(0));
          ctr2 = 0;
        }
        c = (rand() % 100 + 20) % 18 + rand() % 7;
        graph[j].adj.push_back(make_pair(i, c));
        graph_rev[i].adj.push_back(make_pair(j, c));

        cout << j << " " << i << " " << c << "\n";
        ctr2++;
      }
    }

    Bidijkstra bidij(n, graph, graph_rev);
    srand(time(0));
    int ctr = 1;
    later = rand() % 5;

    for (int i = 1; i <= n; ++i) {
      for (int j = i+1; j <= n; j++) {
        ll u = i, v = j;
        //scanf("%d%d", &u, &v);
        ll dist1 = dijkstra(graph, u-1, v-1);
        printf("%d) %lld : %lld -> %lld \t", ctr, u-1, v-1, dist1);
        ll dist2 = bidij.Algo(u-1, v-1);
        if (dist1 != dist2) {
          exit(0);
        }
        ctr++;

        dist1 = dijkstra(graph, v-1, u-1);
        printf("%d) %lld : %lld -> %lld \t", ctr, v-1, u-1, dist1);
        dist2 = bidij.Algo(v-1, u-1);
        if (dist1 != dist2) {
          exit(0);
        }
        ctr++;
      }
    }

    cout << endl << endl;
  }*/
}