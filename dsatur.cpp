#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <filesystem>

using namespace std;

// !! this file is adapted from an example provided by prof. Willem-Jan van Hoeve

// Graph structure holding number of vertices and an adjacency list.
struct Graph {
    int numVertices;
    vector< set<int> > adj;  // Note the space between > and > to support older standards.
};

// Function to read a DIMACS file and build the graph.
// It expects a line starting with "p" to indicate the number of vertices and edges,
// and lines starting with "e" that list the edges. Vertices in DIMACS files are 1-indexed.
Graph readDIMACS(const string &filename) {
    ifstream infile(filename.c_str());
    if (!infile) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }

    string line;
    int n = 0, m = 0;
    Graph g;

    while(getline(infile, line)) {
        if (line.empty() || line[0] == 'c') // skip comments and empty lines
            continue;
        if (line[0] == 'p') {
            // Expected format: p edge n m
            string dummy, edge;
            istringstream iss(line);
            iss >> dummy >> edge >> n >> m;
            g.numVertices = n;
            g.adj.resize(n);  // initialize the adjacency list for n vertices
        } else if (line[0] == 'e') {
            int u, v;
            char c;
            istringstream iss(line);
            iss >> c >> u >> v;
            // Convert from 1-indexed to 0-indexed
            u--; v--;
            if(u < g.numVertices && v < g.numVertices){
                g.adj[u].insert(v);
                g.adj[v].insert(u);
            }
        }
    }
    return g;
}

// DSATUR algorithm implementation
// The algorithm selects at each step an uncolored vertex with the highest saturation
// (number of distinct colors among its neighbors) and, as a tie-breaker, with the highest degree.
void dsatur(string filename) {
    Graph g = readDIMACS(filename);
    int n = g.numVertices;

    // color[i] holds the assigned color for vertex i (-1 means uncolored).
    vector<int> color(n, -1);
    // neighborColors[i] holds the set of colors used by the neighbors of vertex i.
    vector< set<int> > neighborColors(n);
    // colored[i] indicates whether vertex i has been colored.
    vector<bool> colored(n, false);

    vector<int> order(n);

    // Start DSATUR by selecting the vertex with maximum degree for the first assignment.
    int start = 0;
    int maxDegree = -1;
    for (int i = 0; i < n; i++) {
        int deg = g.adj[i].size();
        if (deg > maxDegree) {
            maxDegree = deg;
            start = i;
        }
    }

    int count = 0;
    order[count] = start;
    count++;

    // Color the starting vertex with color 0.
    color[start] = 0;
    colored[start] = true;
    // Update neighbor saturation using iterators.
    for (set<int>::iterator it = g.adj[start].begin(); it != g.adj[start].end(); ++it) {
        neighborColors[*it].insert(0);
    }

    int numColored = 1;

    // Continue until all vertices are colored.
    while (numColored < n) {
        int chosen = -1;
        int maxSat = -1;
        int maxDeg = -1;

        // Select the uncolored vertex with the highest saturation.
        // Tiebreak: vertex with the highest degree.
        for (int i = 0; i < n; i++) {
            if (!colored[i]) {
                int sat = neighborColors[i].size();
                int deg = g.adj[i].size();
                if (sat > maxSat || (sat == maxSat && deg > maxDeg)) {
                    maxSat = sat;
                    maxDeg = deg;
                    chosen = i;
                }
            }
        }

        // Assign the smallest available color to the chosen vertex.
        int c = 0;
        while (neighborColors[chosen].find(c) != neighborColors[chosen].end()) {
            c++;
        }
        color[chosen] = c;
        colored[chosen] = true;
        numColored++;

        order[count] = chosen;
        count++;

        // Update the saturation sets for the uncolored neighbors using iterators.
        for (set<int>::iterator it = g.adj[chosen].begin(); it != g.adj[chosen].end(); ++it) {
            if (!colored[*it]) {
                neighborColors[*it].insert(c);
            }
        }
    }

    // Determine the total number of colors used.
    int numColors = 0;
    for (int i = 0; i < n; i++) {
        if (color[i] > numColors)
            numColors = color[i];
    }

    cout << "Number of colors used: " << numColors + 1 << endl;


    vector<int> reverse_order(n);
    for (int i = 0; i < n; i++) {
        reverse_order[order[i]] = i;
    }
    // cout<<endl;
    // for (auto i = reverse_order.begin(); i != reverse_order.end(); ++i)
    //     std::cout << *i << ' ';
    // cout<<endl;

    vector< set<int> > new_neighbours(n);
    for (int i = 0; i < n; i++) {
        for (int k : g.adj[order[i]]) {
            new_neighbours[i].insert(reverse_order[k]);
        }
    }

    size_t pos = filename.find_last_of("/\\");
    string log_file = (pos == string::npos) ? filename : filename.substr(pos + 1);
    pos = log_file.find_last_of(".");
    log_file = (pos == string::npos) ? log_file : log_file.substr(0, pos);

    log_file = "./data/generated/toy/" + log_file + ".py";              // !! change target directory here
    ofstream outFile(log_file);
    cout <<log_file << endl;

    if (!outFile) {
        cerr << "cant find" << log_file << endl;
        return;
    }

    outFile << "N = " << n << "\nE = { ";
    bool first = true;

    for (int i = 0; i < n; i++) {
        for (int j : new_neighbours[i]) {
            if (i < j) {
                if (!first) outFile << ",";
                outFile << "(" << i << "," << j << ")";
                first = false;
            }
        }
    }
    outFile << " }\n";
    outFile << "K = " << numColors+1 << endl;
}

int main() {
    // change target directory in line 173 for where all processed DSATUR .py files will end up

    // path to directory with all .clq files
    string path = "./data/generated/toy";
    for (const auto &entry : std::filesystem::directory_iterator(path)) {
        dsatur(entry.path().string());
    }
    return 0;
}
