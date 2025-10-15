#include <iostream>
#include <memory>
#include <vector>
using namespace std;


int main() {
        vector<shared_ptr<int>> v1;
        int a = 2;
        int* b = &a;
        v1.push_back(shared_ptr<int> (new int(5)));
        shared_ptr<int> c (new int(3));
        cout << *v1[0];
        *v1[0] = 4;
        cout << a;
        cout << *v1[0];
        v1[0] = c;
        cout << *v1[0];
}