#include <stdio.h>
#include <iostream>
#include <memory>
#include <vector>
#include <stdlib.h>


struct Vertex {
        char state;
        std::vector<Vertex *> edges;
        Vertex();
        Vertex(char start);
};


Vertex::Vertex()
{
        state = 'N';
}


Vertex::Vertex(char start)
{
        state = start;
}


std::ostream& operator<<(std::ostream& os, const Vertex& vOut)
{
        os << "State: " << vOut.state << std::endl;
        os << "Edges:\n";
        for (int i = 0; i < vOut.edges.size(); i++) {
                os << "Edge " << i << ": " << vOut.edges[i] << std::endl;
        }
        return os;
}


int main()
{
        std::vector<std::vector<Vertex *>> lattice(4);
        
        for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                        lattice[i].push_back(new Vertex('S'));
                }
        }

        for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                        std::cout << "Vertex (" << i << ", " << j << "):\n";
                        std::cout << *lattice[i][j];
                }
        }
        
        for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                        delete lattice[i][j];
                }
        }
}