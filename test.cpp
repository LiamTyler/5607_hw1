#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    srand(time(NULL));

    auto f = [](int c) {
        double x = rand() / (float) RAND_MAX * 510 - 255;
        return c + x;
    };
    cout << f(0) << endl;

    return 0;
}
