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
    { 
      workset_.reserve(2*n);
      proc_.reserve(2*n); 
    } 
    // Preserves this memory from the start so that not much time is wasted everytime while increasing vector size

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
};


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
}