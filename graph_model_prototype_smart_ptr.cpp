#include <stdio.h>
#include <iostream>
#include <memory>
#include <vector>
#include <stdlib.h>


struct Vertex {
        char state;
        std::vector<std::weak_ptr<Vertex>> edges;
        Vertex();
        Vertex(char start);
};


Vertex::Vertex()
{
        state = 'S';
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
                os << "Edge " << i << ": " << vOut.edges[i].lock() << " count: "
                   << vOut.edges[i].use_count() << std::endl;
        }
        return os;
}

int main()
{
        std::vector<std::vector<std::shared_ptr<Vertex>>> lattice(4);
        
        // Allocate all the shared pointers to vertices
        for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                        std::shared_ptr<Vertex> vtemp(new Vertex('S')); 
                        lattice[i].push_back(vtemp);
                }
        }

        // Add edges to vertices to make a square lattice
        for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 4; j++) {
                        std::weak_ptr<Vertex> neighbourTemp (lattice[i+1][j]);
                        lattice[i][j]->edges.push_back(neighbourTemp);
                }
        }

        for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 3; j++) {
                        std::weak_ptr<Vertex> neighbourTemp (lattice[i][j+1]);
                        lattice[i][j]->edges.push_back(neighbourTemp);
                }
        }

        for (int i = 1; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                        std::weak_ptr<Vertex> neighbourTemp (lattice[i-1][j]);
                        lattice[i][j]->edges.push_back(neighbourTemp);
                }
        }

        for (int i = 0; i < 4; i++) {
                for (int j = 1; j < 4; j++) {
                        std::weak_ptr<Vertex> neighbourTemp (lattice[i][j-1]);
                        lattice[i][j]->edges.push_back(neighbourTemp);
                }
        }


        // Output all vertex info
        for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                        std::cout << "Vertex (" << i << ", " << j << "): "
                                  << lattice[i][j] << std::endl;
                        std::cout << *lattice[i][j];
                }
        }

        /* leftover deletion code from raw pointer implementation
        for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                        delete lattice[i][j];
                }
        }
        */
}