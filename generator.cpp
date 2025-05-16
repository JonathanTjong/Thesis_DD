#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <random>

// pick random set of elements from set (Fisher-Yates shuffle)
std::vector<int> random_selection(std::vector<int> edges, int density_nr) {
    std::vector<int> shuffled;
    int count = 0;
    while (count < density_nr || !edges.empty()) {
        int k = rand() % edges.size();
        shuffled.push_back(edges[k]);
        edges.erase(edges.begin() + k);
        count++;
    }
    return shuffled;
}


int main() {

    std::random_device rd;  //seed
    // std::mt19937 gen(rd());
    std::mt19937 gen(0);        //TODO change seed back
    std::uniform_real_distribution<> real_d(0, 1);

    std::vector<int> sizes = {100};    // graph size(s)
     std::vector<double> densities = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};    // graph densities
     int redundancy = 25; //nr of times same instance type is generated (but randomized)
     std::string suite = "suite6_n75";    // target folder name, need create empty one beforehand

    for (int r = 0; r < redundancy; r++) {
        for (int size : sizes) {
            // vector of all edges (n(n-1)/2)
            int nr_e = (size * (size-1)) / 2;
            std::vector<int> total_e(nr_e);
            for(int i=0; i<nr_e; ++i)
                total_e[i] = i;
            for (double density : densities) {

                std::vector<int> selected = random_selection(total_e, density * nr_e);

                // need to know exact number of edges at start of writing
                std::vector<int> edges(size*size, 0);
                int idx = 0;
                if (density > 0) {
                    for (int i=0; i<density * nr_e; i++) {
                        edges[selected[i]] = 1;
                    }
                }

                std::ofstream myfile;
                std::stringstream stream;
                stream << std::fixed << std::setprecision(3) << density;
                std::string density_string = stream.str();
                std::string title = "./data/generated/" + suite + "/gen_n" +  std::to_string(size) + "_d0_" + density_string.erase(0, 2) + "r" + std::to_string(r) + ".clq";
                myfile.open(title);

                std::cout << "size: " + std::to_string(size) + " density: " + std::to_string(density) + " r: " + std::to_string(r) + " open: " + std::to_string(myfile.is_open()) << std::endl; // print if could create/open location

                //write 1st line
                int n_edges = std::reduce(edges.begin(), edges.end());
                myfile << "p edge " << size << " " << n_edges << std::endl;
                int count = 0;
                for (int i = 1; i < size; i++) {
                    for (int j = i+1; j < size+1; j++) {
                        // myfile << i << " " << j << std::endl;
                        if (edges[count] == 1) {
                            myfile << "e " <<  i << " " << j << std::endl;
                        }
                        count++;
                    }
                }
                myfile.close();
            }
        }
    }
}
