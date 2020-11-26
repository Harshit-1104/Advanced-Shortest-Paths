#include <cstdio>
#include <cassert>
#include <vector>
#include <queue>
#include <limits>
#include <utility>
#include <math.h>
#include <iostream>
#include <fstream>
//#include <Windows.h>
//#include <time.h>
//#include <limits.h>
#include <algorithm>

using namespace std;
typedef long long ll;

const ll INFINI = numeric_limits<ll>::max() / 4;

struct Node {
  ll id;
  ll dist_act = INFINI;
  ll dist_pi = INFINI;
  vector<pair<ll, ll>> adj;
  bool processed = false;
  bool wp = false;

  pair<ll, ll> cord;
};

struct CompareHeight { 
	bool inline operator()(Node* const& n1, Node* const& n2) 
	{ 
		// return "true" if "p1" is ordered 
		// before "p2", for example: 
		return n1->dist_pi > n2->dist_pi; 
	} 
}; 

class AStar {
    // See the descriptions of these fields in the starter for friend_suggestion
    ll n_;
    ll s_, t_;
    vector<Node> graph_, graph_rev_;

    vector<Node*> proc_;
    vector<Node*> workset_;
    
public:
    AStar(ll n, vector<Node> graph, vector<Node> graph_rev)
        : n_(n), graph_(graph), graph_rev_(graph_rev) { 
      workset_.reserve(2*n); 
      proc_.reserve(2*n);
    }

    inline ll dist(pair<ll, ll> s, pair<ll, ll> t) {
      return (sqrt(pow(s.first-t.first, 2) + pow(s.second-t.second, 2)));
    }

    // See the description of this method in the starter for friend_suggestion
    void clear() {
        Node* v;
        for (ll i = 0; i < workset_.size(); ++i) {
            v = workset_[i];
            graph_[v->id].wp = graph_rev_[v->id].wp = graph_[v->id].processed = graph_rev_[v->id].processed = false;
            graph_[v->id].dist_act = graph_rev_[v->id].dist_act = INFINI;
            graph_[v->id].dist_pi = graph_rev_[v->id].dist_pi = INFINI;
        }
        workset_.clear();
        proc_.clear();
    }

    // Processes visit of either fwd or bwd search to node v trying to
    // relax the current distance by dist.
    ll Algo(ll s, ll t) {
        clear();
        s_ = s;
        t_ = t;

        graph_[s].dist_act = graph_rev_[t].dist_act = 0;
        graph_[s].dist_pi = graph_rev_[t].dist_pi = 0;
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
              tmp->wp = graph_rev_[tmp->id].wp = true;       // Add for the other graph too
              workset_.emplace_back(tmp);
            }

            if (graph_rev_[tmp->id].processed == true) {
                return ShortestPath(s, t); 
            }

            tmp = bwd.top();
            bwd.pop();

            Process_bwd(tmp, &bwd);
            
            if (tmp->wp == false) {
              tmp->wp = graph_[tmp->id].wp = true;
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
      
        vector<ll> potential = potential_func_(u->adj);
        for (ll i = 0; i < u->adj.size(); ++i) {
            pair<ll, ll> *node = &u->adj[i];
            if (graph_[node->first].dist_pi > u->dist_act + node->second + potential[i]) {
                graph_[node->first].dist_act = u->dist_act + node->second;
                graph_[node->first].dist_pi = graph_[node->first].dist_act + potential[i];
              
                if (graph_[node->first].wp == false) {
                  graph_[node->first].wp = graph_rev_[node->first].wp = true;
                  workset_.emplace_back(&graph_[node->first]);
                }

                fwd->emplace(&graph_[node->first]);
            }
        }
        graph_[u->id].processed = true;
        proc_.emplace_back(&graph_[u->id]);
    }

    void Process_bwd(Node *u, priority_queue<Node*, vector<Node*>, CompareHeight> *bwd) {
        
        vector<ll> potential = potential_func_(u->adj, -1);
        for (ll i = 0; i < u->adj.size(); ++i) {
            pair<ll, ll>* node = &u->adj[i];
            if (graph_rev_[node->first].dist_pi > u->dist_act + node->second + potential[i]) {
                graph_rev_[node->first].dist_act = u->dist_act + node->second;
                graph_rev_[node->first].dist_pi = graph_rev_[node->first].dist_act + potential[i];
                
                if (graph_rev_[node->first].wp == false) {
                  graph_rev_[node->first].wp = graph_[node->first].wp = true;
                  workset_.emplace_back(&graph_rev_[node->first]);
                }

                bwd->emplace(&graph_rev_[node->first]);
            }
        }
        graph_rev_[u->id].processed = true;
        proc_.emplace_back(&graph_rev_[u->id]);
    }

    vector<ll> potential_func_(vector<pair<ll, ll>> adj, int x = 1) {
      vector<ll> tmp;
      ll size = adj.size();
      tmp.reserve(2*size);

      for (int i = 0; i < size; ++i) {
        pair<ll, ll> s = graph_[adj[i].first].cord;
        tmp.emplace_back((-dist(s, graph_[t_].cord) + dist(s, graph_rev_[s_].cord)) * x / 2);     // Coordinates remains the same in both the graphs
      }
      return tmp;
    }

    ll ShortestPath(ll s, ll t) {
        ll distance = INFINI;
        Node* u_best;

        for (ll i = 0; i < proc_.size(); i++) {
            Node* u = proc_[i];
       
            if (graph_[u->id].dist_act + graph_rev_[u->id].dist_act < distance) {
                u_best = u;
                distance = graph_rev_[u->id].dist_act + graph_[u->id].dist_act;
            }
        }

        if (distance >= INFINI) {
            return -1;
        }
        return distance;
    }

/*    ll distance(vector<Node> graph, ll s, ll t) {
    s_ = s;
    t_ = t;  

    graph[s].dist_act = 0; // Source
    graph[s].dist_pi = 0;
    priority_queue<Node*, vector<Node*>, CompareHeight> Q; // The keys will be distance values
    Q.emplace(&graph[s]);

    while (!Q.empty()) {
      Node* u = Q.top();
      Q.pop();

      //cout << "hello";

      if (graph[u->id].processed == true) {
        continue;
      }

      for (int i = 0; i < graph[u->id].adj.size(); i++) {
        pair<ll, ll>* Node = &u->adj[i];
        double potential = potential_func_(graph[Node->first].cord);
        if (graph[Node->first].dist_pi > u->dist_act + Node->second + potential) {
          graph[Node->first].dist_act = u->dist_act + Node->second;
          graph[Node->first].dist_pi = graph[Node->first].dist_act + potential;

          //graph[Node->first].prev = graph[u].id;
          Q.emplace(&graph[Node->first]);
        }
      }

      graph[u->id].processed = true;
    }

    return (graph[t].dist_act != INFINI) ? (graph[t].dist_act) : (-1);
    }

    
    ll Algo(ll s, ll t) {
        clear();
        s_ = s;
        t_ = t;

        graph_[s].dist_act = graph_rev_[t].dist_act = 0;
        graph_[s].dist_pi = graph_rev_[t].dist_pi = 0;

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
            if (graph_[tmp->first].dist_act > u->dist_act + tmp->second) {
                graph_[tmp->first].dist_act = u->dist_act + tmp->second;
                graph_[tmp->first].dist_pi = graph_[tmp->first].dist_act;

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
            if (graph_rev_[tmp->first].dist_act > u->dist_act + tmp->second) {
                graph_rev_[tmp->first].dist_act = u->dist_act + tmp->second;
                graph_rev_[tmp->first].dist_pi = graph_rev_[tmp->first].dist_act;
                
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
        ll distance = INFINI;
        Node* u_best;

        for (ll i = 0; i < proc_.size(); i++) {
            Node* u = proc_[i];

            if (graph_[u->id].dist_act + graph_rev_[u->id].dist_act < distance) {
                u_best = u;
                distance = graph_rev_[u->id].dist_act + graph_[u->id].dist_act;
            }
        }

        if (distance >= INFINI) {
            return -1;
        }
        return distance;
    }*/
};
/*
class Alt {
  // See the descriptions of these fields in the starter for friend_suggestion
  ll n_;
  ll s_, t_;
  vector<Node> graph_, graph_rev_;

  vector<Node*> proc_;
  vector<Node*> workset_;

public:

    Alt(ll n, vector<Node> graph, vector<Node> graph_rev)
        : n_(n), graph_(graph), graph_rev_(graph_rev) { 
      workset_.reserve(2*n); 
      proc_.reserve(2*n);
    }
  
    void clear() {
        Node* v;
        for (ll i = 0; i < workset_.size(); ++i) {
            v = workset_[i];
            graph_[v->id].wp = graph_rev_[v->id].wp = graph_[v->id].processed = graph_rev_[v->id].processed = false;
            graph_[v->id].dist_act = graph_rev_[v->id].dist_act = INFINI;
            graph_[v->id].dist_pi = graph_rev_[v->id].dist_pi = INFINI;
        }
        workset_.clear();
        proc_.clear();
    }

  ll Algo(ll s, ll t) {
        clear();
        s_ = s;
        t_ = t;

        graph_[s].dist_act = graph_rev_[t].dist_act = 0;
        graph_[s].dist_pi = graph_rev_[t].dist_pi = 0;

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
            if (graph_[tmp->first].dist_act > u->dist_act + tmp->second) {
                graph_[tmp->first].dist_act = u->dist_act + tmp->second;
                graph_[tmp->first].dist_pi = graph_[tmp->first].dist_act;

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
            if (graph_rev_[tmp->first].dist_act > u->dist_act + tmp->second) {
                graph_rev_[tmp->first].dist_act = u->dist_act + tmp->second;
                graph_rev_[tmp->first].dist_pi = graph_rev_[tmp->first].dist_act;
                
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
        ll distance = INFINI;
        Node* u_best;

        for (ll i = 0; i < proc_.size(); i++) {
            Node* u = proc_[i];

            if (graph_[u->id].dist_act + graph_rev_[u->id].dist_act < distance) {
                u_best = u;
                distance = graph_rev_[u->id].dist_act + graph_[u->id].dist_act;
            }
        }

        if (distance >= INFINI) {
            return -1;
        }
        return distance;
    }
};*/

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);
//    ifstream fin;
    
    ll n, m;
    scanf("%lld%lld", &n, &m);
    vector<Node> graph(n), graph_rev(n);
    for (int i=0;i<n;++i){
        ll a, b;
        scanf("%lld%lld", &a, &b);
        graph[i].cord.first = graph_rev[i].cord.first = a;
        graph[i].cord.second = graph_rev[i].cord.second = b;
        graph[i].id = graph_rev[i].id = i;
    }

    for (int i=0; i<m; ++i) {
        ll u, v, c;
        scanf("%lld%lld%lld", &u, &v, &c);
        graph[u-1].adj.emplace_back(make_pair(v-1, c));
        graph_rev[v-1].adj.emplace_back(make_pair(u-1, c));
    }

    AStar astar(n, graph, graph_rev);
//    Alt alt(n, graph, graph_rev);
//    fin.open("sample.txt");

    int t;
    ll num;
/*
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i != j) {
          cout << i << " : " << j << " => " << graph[i].adj[j].second << "\t";
          cout << j << " : " << i << " => " << graph_rev[j].adj[i].second << "\n";
        }
      }
    }*/

    scanf("%d", &t);
    for (int i=0; i<t; ++i) {
        ll u, v;
        scanf("%lld%lld", &u, &v);
        ll dist = astar.Algo(u-1, v-1);
        printf("%lld\n", dist);

/*        fin >> num;
        if (num != dist) {
          assert(0);
        } */
    }


/*    
  ll later = 0;
  while (true) {
    Sleep(1000);
    srand(time(0));
    ll n = rand() % 20 + 3;              // Number of vertices
    if (n <= 3) {
      n += later;
    }

    Sleep(1000);
    srand(time(0));
    ll m = (n * (n-2));       // Number of edges
    std::cout << n << " " << m << "\n";
    vector<Node> graph(n), graph_rev(n);
    int ctr2 = 0;
    ll xlim = 100, ylim = 100;

    for (ll i = 0; i < n; i++) {
        graph[i].id = graph_rev[i].id = i;
        graph[i].dist_act = graph_rev[i].dist_act = INFINI;
        graph[i].processed = graph_rev[i].processed = false;

        Sleep(1000);
        srand(time(0));
        ll x = rand() % xlim;

        Sleep(1000);
        ll y = (rand() % ylim + x) % ylim;

        graph[i].cord.first = graph_rev[i].cord.first = x;
        graph[i].cord.second = graph_rev[i].cord.second = y;

        cout << "Coordinates of " << i << "th node : ("<< x << "," << y << ")\n";
    }

    for (ll i = 0; i < n; i++) {
      for (ll j = i+1; j < n; j++) {
        ll dista = dist(graph[i].cord, graph[j].cord);

        ll c = (rand() % 100 + 50) % (dista +1) + dista;
        graph[i].adj.push_back(make_pair(j, c));
        graph_rev[j].adj.push_back(make_pair(i, c));

        //cout << "Coordinates of " << i << " node: (" << graph[i].cord.first << "," << graph[i].cord.second << ")\t" << j << " node: (" << graph[j].cord.first << "," << graph[j].cord.second << ")\n";
        cout << i << " " << j << " " << c << "\n";
        
        c = (rand() % 100 + 20) % (dista +1) + dista;
        graph[j].adj.push_back(make_pair(i, c));
        graph_rev[i].adj.push_back(make_pair(j, c));

        cout << j << " " << i << " " << c << "\n";
      }
    }

    AStar astar(n, graph, graph_rev);
    Alt alt(n, graph, graph_rev);
    srand(time(0));
    int ctr = 1;
    later = rand() % 5;

    for (int i = 1; i <= n; ++i) {
      for (int j = i+1; j <= n; j++) {
        ll u = i, v = j;
        //scanf("%d%d", &u, &v);
        ll dist1 = astar.Algo(u-1, v-1);
        printf("%d) %lld : %lld -> %lld \n", ctr, u-1, v-1, dist1);
        ll dist2 = alt.Algo(u-1, v-1);
        if (dist1 != dist2) {
          assert(0);
        }
        ctr++;

        dist1 = astar.Algo(v-1, u-1);
        printf("%d) %lld : %lld -> %lld \n", ctr, v-1, u-1, dist1);
        dist2 = alt.Algo(v-1, u-1);
        if (dist1 != dist2) {
          assert(0);
        }
        ctr++;
      }
    }
    cout << endl << endl;
  }*/
}
