#include <iostream>
#include "timer.cpp"
#include <vector>
#include <iomanip>
#include <Windows.h>
#include <algorithm>

using namespace std;
typedef long long ll;

struct Node {
    ll id;
    ll dist;
    vector<pair<ll, ll>> adj;
    bool processed;
    ll prev;

    bool wp = false;
    ll heap_idx = -1;
};

int main() {
  Node a, b;
  a.id = b.id = a.dist = b.dist = a.processed = b.processed = a.prev = b.prev = 1;
  int a1 = 20, b1 = 4000000;

  bool l = true;


  Timer t;
  if (a1) {
    a1 = 2;
  }
  Sleep(1000);
  cout << setprecision(8) << t.elapsed();
}
