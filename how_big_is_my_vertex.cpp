#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <memory>

struct Vertex {
        int state;
        std::vector<std::shared_ptr<Vertex>> edges;
};

struct Test {int a; int c; int d; std::vector<int> b;};


int main()
{
        std::cout<< "test";
        std::vector<int> testvec = {1, 2, 3};
        for (int i = 0; i < testvec.size(); i++) {
                std::cout << testvec[i] << "\n";
        }
        std::cout << std::endl;
        std::cout << sizeof(int);
        std::shared_ptr<Vertex> v1(new Vertex);
        std::shared_ptr<Vertex> v2(new Vertex);
        std::shared_ptr<Vertex> v3(new Vertex);
        std::shared_ptr<Vertex> v4(new Vertex);
        std::shared_ptr<Vertex> v5(new Vertex);
        v1->edges.push_back(v2);
        v1->edges.push_back(v3);
        v1->edges.push_back(v4);
        v1->edges.push_back(v5);

        std::cout << std::endl << sizeof(*v1);
        std::cout << std::endl << v1;
        std::cout << std::endl << sizeof(Vertex);
        std::cout << std::endl << sizeof(std::shared_ptr<Vertex>);
        std::cout << std::endl << sizeof(std::vector<std::shared_ptr<Vertex>>);
        std::cout << std::endl << sizeof(Test);
        std::cout << std::endl << sizeof(std::vector<int>);
}