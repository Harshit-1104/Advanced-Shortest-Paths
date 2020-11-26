#include <iostream>

using namespace std;

int main() {
  int a[13][13];

  for (int i = 0; i < 13; i++) {
    a[i][i] = 0;
  }
  int x, y, w;

  while(cin >> x >> y >> w) {
    a[x][y] = w;
  }

  for (int i = 0; i< 13; i++) {
    for (int j = 0; j< 13; j++) {
      cout << a[i][j] << ", ";
    }
    cout <<endl;
  }
}