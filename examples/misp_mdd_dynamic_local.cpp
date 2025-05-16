#include <queue>

#include "codd.hpp"
#include "math.h"

struct MISPMDD {
   GNSet sel; // set of eligible vertices
   std::vector<std::pair<int, int>> deg; // list of pairs <vertex_label, local degree>
   int transitions; // nr. of transitions from the root

   MISPMDD(GNSet sel, std::vector<std::pair<int, int>> data, int transitions) : sel(std::move(sel)), transitions(transitions), deg(std::move(data)) {}

   friend std::ostream& operator<<(std::ostream& os, const MISPMDD& m) {
      return os << "<" << m.sel << ',' << m.transitions << ">";
   }
};

template<> struct std::equal_to<MISPMDD> {
   bool operator()(const MISPMDD& s1,const MISPMDD& s2) const {
      return s1.sel == s2.sel && s1.transitions == s2.transitions;
   }
};

template<> struct std::hash<MISPMDD> {
   std::size_t operator()(const MISPMDD& v) const noexcept {
      return std::rotl( std::hash<GNSet>{}(v.sel),32 ^ std::hash<int>{}(v.transitions));
   }
};

struct GE {
   int a,b;
   friend bool operator==(const GE& e1,const GE& e2) {
      return e1.a == e2.a && e1.b == e2.b;
   }
   friend std::ostream& operator<<(std::ostream& os,const GE& e) {
      return os  << e.a << "-->" << e.b;
   }
};


struct Instance {
   int nv;
   int ne;
   std::vector<int> weights;
   std::vector<GE> edges;
   FArray<GNSet> adj;
   Instance() : adj(0) {}
   Instance(Instance&& i) : nv(i.nv),ne(i.ne),edges(std::move(i.edges)) {}
   GNSet vertices() {
      return setFrom(std::views::iota(0,nv));
   }
   auto getEdges() const noexcept { return edges;}
   void convert() {
      adj = FArray<GNSet>(nv+1);
      for(const auto& e : edges) {
         adj[e.a].insert(e.b);
         adj[e.b].insert(e.a);
      }
   }
};

Instance readFile(const char* fName)
{
   Instance i;
   using namespace std;
   ifstream f(fName);
   while (!f.eof()) {
      char c;
      f >> c;
      if (f.eof()) break;
      switch(c) {
         case 'n': {
            int index;
            int weight;
            f >> index >> weight;
            i.weights.push_back(weight);
         }
         case 'c': {
            std::string line;
            std::getline(f,line);
         }break;
         case 'p': {
            string w;
            f >> w >> i.nv >> i.ne;
         }break;
         case 'e': {
            GE edge;
            f >> edge.a >> edge.b;
            edge.a--,edge.b--;      // make it zero-based
            assert(edge.a >=0);
            assert(edge.b >=0);
            i.edges.push_back(edge);
         }break;
      }
   }
   f.close();
   i.convert();
   return i;
}

int main(int argc,char* argv[])
{
   // using STL containers for the graph
   if (argc < 2) {
      std::cout << "usage: coloring <fname> <width>\n";
      exit(1);
   }
   const char* fName = argv[1];
   const int w = argc==3 ? atoi(argv[2]) : 64;
   auto instance = readFile(fName);

   const GNSet ns = instance.vertices();
   const int top = ns.size();
   const std::vector<GE> es = instance.getEdges();
   std::cout << "TOP=" << top << "\n";
   std::vector<int> end = {};

   Bounds bnds([&es](const std::vector<int>& inc)  {
      bool ok = true;
      for(const auto& e : es) {
         bool v1In = (std::find(inc.begin(), inc.end(), e.a) != inc.end()) ? true : false;
         bool v2In = (std::find(inc.begin(), inc.end(), e.b) != inc.end()) ? true : false;
         if (v1In && v2In) {
            std::cout << e << " BOTH ep in inc: " << inc << "\n";
            assert(false);
         }
         ok &= !(v1In && v2In);
      }
      std::cout << "CHECKER is " << ok << "\n";
   });
   const auto labels = ns | GNSet { top };     // using a plain set for the labels
   std::vector<int> weight(ns.size()+1);
   // if no weights, set all to 1
   if (instance.weights.empty()) {
      weight[ns.size()] = 0;
      for(auto v : ns) weight[v] = 1;
   }
   else {weight = instance.weights; weight.push_back(0);}
   int beam_n = argc>=4 ? atoi(argv[3]) : 1;

   auto neighbors = instance.adj;
   std::vector<std::pair<int, int>> pq_empty;
   auto empty_state = MISPMDD { GNSet{}, pq_empty, -1};

   const auto myInit = [top,&neighbors]() {   // The root state
      GNSet U = {};
      std::vector<std::pair<int, int>> pq_init;

      for(auto i : std::views::iota(0,top)) {
         U.insert(i);
         pq_init.emplace_back(i, neighbors[i].size());
      }
      return MISPMDD{U, pq_init, 0};
   };
   const auto myTarget = [&empty_state]() {    // The sink state
      return empty_state;
   };
   const auto lgf = [top, beam_n](const MISPMDD& s, DDContext)  {
      if (s.sel.size() <= beam_n) {    // can replace with positives in s.deg <= beam_n
         GNSet u = s.sel;
         return u;
      }
      GNSet out = {top};
      auto deg = s.deg;
      for (int i = 0; i < beam_n; i++) {     // go over deg list, pick max deg vertex, set its deg to -1, repeat beam_n times
         std::pair<int, int> max_pair = std::make_pair(0, -1);
         int max_idx = 0;
         for (int j = 0; j < deg.size(); j++) {
            if (deg[j].second > max_pair.second) {
               max_pair = deg[j];
               max_idx = j;
            }
         }
         out.insert(max_pair.first);
         deg[max_idx].second = -1;
      }
      return out;
   };
   auto myStf = [top, &neighbors, &empty_state, beam_n](const MISPMDD& s,const int label) -> std::optional<MISPMDD> {
      if (s.transitions == -1) {
         if (s.sel.empty()) {
            return empty_state;
         }
         else {
            std::cout << "top error" << "\n";
            std::cout << s.sel << "\n";
            return std::nullopt;
         }
      }else {
         if (label == top) {     // extra transition

            if (s.sel.size() <= beam_n) return empty_state;

            GNSet out = s.sel;
            GNSet removed = {};
            auto deg = s.deg;

            // remove beam
            for (int i = 0; i < beam_n; i++) {
               std::pair<int, int> max_pair = std::make_pair(0, -1);
               int max_idx = 0;
               for (int j = 0; j < deg.size(); j++) {
                  if (deg[j].second > max_pair.second) {
                     max_pair = deg[j];
                     max_idx = j;
                  }
               }
               out.remove(max_pair.first);
               deg[max_idx].second = -1;
               removed.insert(max_pair.first);
            }

            // remove encountered vertices from deg
            for (auto it = deg.begin(); it != deg.end(); ) {
               if (it->second == -1) {
                  it = deg.erase(it);
               } else {
                  it++;
               }
            }

            GNSet adj_nb = {};
            for (int r : removed) {
               adj_nb = adj_nb | (neighbors[r] & out);
            }

            // change degrees neighbours of those deleted
            if (!adj_nb.empty()) {
               for (auto& elem : deg) {
                  if (adj_nb.contains(elem.first)) {
                     elem.second -= (neighbors[elem.first] & removed).size();
                  }
               }
            }

            if (out.empty()) return empty_state;
            return MISPMDD { std::move(out), deg, s.transitions+1};
         }

         if (!s.sel.contains(label)) return std::nullopt; // we cannot take n (label==1) if not legal.

         GNSet out = s.sel;
         GNSet removed = {};
         auto deg = s.deg;

         // remove up to label
         for (int i = 0; i < beam_n; i++) {
            std::pair<int, int> max_pair = std::make_pair(0, -1);
            int max_idx = 0;
            for (int j = 0; j < deg.size(); j++) {
               if (deg[j].second > max_pair.second) {
                  max_pair = deg[j];
                  max_idx = j;
               }
            }
            out.remove(max_pair.first);
            deg[max_idx].second = -1;
            if (max_pair.first == label) break;
            removed.insert(max_pair.first);
         }

         // remove neighbors label vertex
         out.diffWith(neighbors[label]);
         removed = removed | (s.sel & neighbors[label]);

         for (auto it = deg.begin(); it != deg.end(); ) {
            if (it->second == -1 || neighbors[label].contains(it->first)) {
               it = deg.erase(it);
            } else {
               it++;
            }
         }

         GNSet adj_nb = {};
         for (int r : removed) {
            adj_nb = adj_nb | (neighbors[r] & out);   // was approx before
         }

         if (!adj_nb.empty()) {
            for (auto& elem : deg) {
               if (adj_nb.contains(elem.first)) {
                  elem.second -= (neighbors[elem.first] & removed).size();
               }
            }
         }

         if (out.empty()) return empty_state;

         return MISPMDD { std::move(out), deg, s.transitions+1};
      }
   };
   const auto scf = [weight, top](const MISPMDD& s,int label) { // cost function
      if (label == top) return 0;

      return weight[label];
   };
   const auto smf = [&neighbors](const MISPMDD& s1,const MISPMDD& s2) -> std::optional<MISPMDD> { // merge function
      if (s1.transitions != s2.transitions)
         return std::nullopt;

      std::map<int, int> merged_map;
      GNSet un = (s1.sel | s2.sel);

      for (int v : un) {
         merged_map[v] = (neighbors[v] & un).size();
      }

      std::vector<std::pair<int, int>> new_deg;
      for (auto [key, deg] : merged_map) {
         new_deg.emplace_back(key, deg);
      }

      return MISPMDD{
         un, new_deg, s1.transitions};
   };
   const auto eqs = [](const MISPMDD& s) -> bool {
      return s.sel.size() == 0;
   };
   // can uncomment first line to disable new local bound
   const auto local = [&weight, &neighbors, &es, top](const MISPMDD& s, LocalContext) -> double {
      // return sum(s.sel,[&weight](auto v) { return weight[v];});

      int maxD = 0;           // unweighted vertices, no weight considered yet
      int edg = 0;
      for (auto& elem : s.deg) {
         if (elem.second > maxD) maxD = elem.second;
         edg += elem.second;
      }
      edg = edg / 2.0;
      if (maxD == 0) return sum(s.sel,[&weight](auto v) { return weight[v];});
      else return s.sel.size() - std::floor(edg/maxD);
      // else return s.sel.size() - std::ceil(edg/maxD);
   };
   const auto sDom = [](const MISPMDD& a,const MISPMDD& b) -> bool {
      return  a.sel == b.sel;
   };

   BAndB engine(DD<MISPMDD,Maximize<double>, // to maximize
                //decltype(myInit),
                decltype(myTarget), // MISPMDD(*)(),
                decltype(lgf),
                decltype(myStf),
                decltype(scf),
                decltype(smf),
                decltype(eqs),
                decltype(local)
                >::makeDD(myInit,myTarget,lgf,myStf,scf,smf,eqs,labels,local),w);
   engine.search(bnds);
   return 0;
}

