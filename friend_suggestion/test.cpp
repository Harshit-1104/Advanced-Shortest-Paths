#include <iostream>
#include "timer.cpp"

int main()
{
    Timer t;

    cout << "Time elapsed: " << t.elapsed() << " seconds\n";

    return 0;
}