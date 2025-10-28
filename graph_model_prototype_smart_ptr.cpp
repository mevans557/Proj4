#include <stdio.h>
#include <iostream>
#include <memory>
#include <vector>
#include <stdlib.h>

#define WIDTH 10
#define HEIGHT 10

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


bool connected(const std::shared_ptr<Vertex>& v1,
               const std::shared_ptr<Vertex>& v2)
{
        for (int i = 0; i < v1->edges.size(); i++) {
                if (v1->edges[i].lock() == v2) return true;
        }
        return false;
}

// TODO make this nicer
void draw_graph(std::ostream& os,
                const std::vector<std::vector<std::shared_ptr<Vertex>>>& graph)
{
        for (int i = 0; i < graph.size()-1; i++) {
                std::string bufferline = "";
                for (int j = 0; j < graph[i].size()-1; j++) {
                        os << graph[i][j]->state;
                        if (connected(graph[i][j], graph[i][j+1])) os << "-";
                        else os << " ";
                        if (connected(graph[i][j], graph[i+1][j])) bufferline += "|";
                        else bufferline += " ";
                        bufferline += " ";
                }
                os << graph[i][graph[i].size()-1]->state << std::endl;
                if (connected(graph[i][graph[i].size()-1], graph[i+1][graph[i].size()-1])) bufferline += "|";
                else bufferline += " ";
                os << bufferline << std::endl;
        }
        for (int j = 0; j < graph[graph.size()-1].size()-1; j++) {
                os << graph[graph.size()-1][j]->state;
                if (connected(graph[graph.size()-1][j], graph[graph.size()-1][j+1])) os << "-";
                else os << " ";
        }
        os << graph[graph.size()-1][graph[graph.size()-1].size()-1]->state;
}


int main()
{
        std::vector<std::vector<std::shared_ptr<Vertex>>> lattice(HEIGHT);
        
        // Allocate all the shared pointers to vertices
        for (int i = 0; i < HEIGHT; i++) {
                for (int j = 0; j < WIDTH; j++) {
                        std::shared_ptr<Vertex> vtemp(new Vertex('S')); 
                        lattice[i].push_back(vtemp);
                }
        }

        // Add edges to vertices to make a square lattice
        for (int i = 0; i < HEIGHT-1; i++) {
                for (int j = 0; j < WIDTH; j++) {
                        std::weak_ptr<Vertex> neighbourTemp (lattice[i+1][j]);
                        lattice[i][j]->edges.push_back(neighbourTemp);
                }
        }

        for (int i = 0; i < HEIGHT; i++) {
                for (int j = 0; j < WIDTH-1; j++) {
                        std::weak_ptr<Vertex> neighbourTemp (lattice[i][j+1]);
                        lattice[i][j]->edges.push_back(neighbourTemp);
                }
        }

        for (int i = 1; i < HEIGHT; i++) {
                for (int j = 0; j < WIDTH; j++) {
                        std::weak_ptr<Vertex> neighbourTemp (lattice[i-1][j]);
                        lattice[i][j]->edges.push_back(neighbourTemp);
                }
        }

        for (int i = 0; i < HEIGHT; i++) {
                for (int j = 1; j < WIDTH; j++) {
                        std::weak_ptr<Vertex> neighbourTemp (lattice[i][j-1]);
                        lattice[i][j]->edges.push_back(neighbourTemp);
                }
        }


        // Output all vertex info
        for (int i = 0; i < HEIGHT; i++) {
                for (int j = 0; j < WIDTH; j++) {
                        std::cout << "Vertex (" << i << ", " << j << "): "
                                  << lattice[i][j] << std::endl;
                        std::cout << *lattice[i][j];
                }
        }

        draw_graph(std::cout, lattice);


        /* leftover deletion code from raw pointer implementation
        for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                        delete lattice[i][j];
                }
        }
        */
}