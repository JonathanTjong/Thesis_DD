#include <queue>
#include <random>

#include "codd.hpp"
#include "math.h"

struct MISPMDD {
   struct Compare {
      bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b) const {
         return a.first < b.first;
      }
   };

   GNSet sel;
   std::set<std::pair<int, int>, Compare> pq;
   int transitions;

   MISPMDD(GNSet sel, std::vector<std::pair<int, int>> data, int transitions) : sel(std::move(sel)), transitions(transitions) {
      pq.insert(data.begin(), data.end());
   }
   MISPMDD(GNSet sel, std::set<std::pair<int, int>, Compare> data, int transitions) : sel(std::move(sel)), transitions(transitions), pq(data) {}

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

// select index with highest connectivity, ties broken with highest degree
int get_max_connectivity(const std::vector<int> connectivity, const std::vector<int> degrees) {
   int max_idx = 0;
   int max_conn = INT_MIN;
   int max_deg = INT_MIN;

   for (int i = 0; i < connectivity.size(); i++) {
      if (connectivity[i] > max_conn || (connectivity[i] == max_conn && degrees[i] > max_deg)) {
         max_conn = connectivity[i];
         max_deg = degrees[i];
         max_idx = i;
      }
   }
   return max_idx;
}

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

   // ========== static variable orderings =========================================================================

   // << this variable ordering was not used in the end, as it did not perform as well >>
   // min degree heuristic static ordering
   //    degrees: for each v in V, degree/nr of neighbors of v
   //    neighbors: for each v in V, set of neighboring vertices
   auto min_deg_heuristic = [top](std::vector<int> degrees, FArray<GNSet> neighbors) {
      auto deg = degrees;
      auto nb = neighbors;
      std::vector<int> result(top);
      for (int i = 0; i < top; i++) {
         auto it = std::min_element(std::begin(deg), std::end(deg));
         int idx = std::distance(std::begin(deg), it);
         result[i] = idx;
         deg[idx] = INT_MAX;
         for (auto n : nb[idx]) {
            nb[n].remove(idx);
            nb[n] = (nb[n] | nb[idx]).remove(n); // new neighbors of neighbor n = union of nb[i] and old nb[n] w/o n
            deg[n]--;
         }
      }

      return result;
   };

   // The "max degree" ordering in the thesis is called the "simple max deg" in this code
   // This max_deg here is based on the min degree heuristic above
   // << this variable ordering was not used in the end, as it did not perform as well >>
   // max degree static ordering (for experimentation)
   auto max_deg = [top](std::vector<int> degrees, FArray<GNSet> neighbors) {
      auto deg = degrees;
      auto nb = neighbors;
      std::vector<int> result(top);

      // max deg heuristic
      for (int i = 0; i < top; i++) {
         auto it = std::max_element(std::begin(deg), std::end(deg));
         int idx = std::distance(std::begin(deg), it);
         result[i] = idx;
         deg[idx] = -top*top; // MIN_VALUE
         for (auto n : nb[idx]) {
            nb[n].remove(idx);
            nb[n] = (nb[n] | nb[idx]).remove(n);
            deg[n]--;
         }
      }
      return result;
   };

   auto min_width = [top](std::vector<int> degrees, FArray<GNSet> neighbors) {
      auto deg = degrees;
      auto nb = neighbors;
      std::vector<int> result(top);
      for (int i = 0; i < top; i++) {
         auto it = std::min_element(std::begin(deg), std::end(deg));
         int idx = std::distance(std::begin(deg), it);
         result[top-i-1] = idx; // put results from back to front of list
         deg[idx] = INT_MAX;
         for (auto n : nb[idx]) {
            nb[n].remove(idx);
            deg[n]--;
         }
      }

      return result;
   };

   // max connectivity
   //among all unselected vertices, select the one that is connected to the most vertices that have been selected so far.
   //In case of ties, select a vertex with the highest degree
   auto max_connectivity = [top](std::vector<int> degrees, FArray<GNSet> neighbors) {
      // auto deg = degrees;
      std::vector<int> connectivity(top); // keep track of connectivity
      auto nb = neighbors;
      std::vector<int> result(top);

      // select first vertex
      auto it = std::max_element(std::begin(degrees), std::end(degrees));
      int idx = std::distance(std::begin(degrees), it);
      result[0] = idx;
      connectivity[idx] = -top*top; // MIN_VALUE
      for (auto n : nb[idx]) {
         nb[n].remove(idx);
         connectivity[n]++;
      }

      // create set of unselected vertices
      GNSet unselected = {};
      for(auto i : std::views::iota(0,top)) {   //exclude 0, already covered
         unselected.insert(i);
      }
      unselected.remove(idx);
      int count = 1;
      while (!unselected.empty()) {
         int idx = get_max_connectivity(connectivity, degrees);
         result[count] = idx;
         unselected.remove(idx);
         count++;
         connectivity[idx] = -top*top; // MIN_VALUE
         for (auto n : nb[idx]) {
            nb[n].remove(idx);
            connectivity[n]++;
         }
      }

      return result;
   };

   auto max_paths = [top](FArray<GNSet> neighbors) {
      auto nb = neighbors;
      std::vector<int> result(top);

      // create set of unselected vertices
      GNSet unselected = {};
      for(auto i : std::views::iota(0,top)) {
         unselected.insert(i);
      }
      int count = 0;
      while (!unselected.empty()) {
         int last_v = min(unselected); // get an unselected vertex
         result[count] = last_v;
         unselected.remove(last_v);
         count++;
         for (auto n : neighbors[last_v]) {
            nb[n].remove(last_v);
         }

         int first_v = last_v;

         // while path not maximal
         while (true) {
            // try add to last vertex in path
            if (!nb[last_v].empty()) {
               last_v = min(nb[last_v]);
               result[count] = last_v;
               unselected.remove(last_v);
               count++;
               for (auto n : neighbors[last_v]) {
                  nb[n].remove(last_v);
               }
            }
            // try add to first vertex in path
            else if (!nb[first_v].empty()) {
               first_v = min(nb[first_v]);
               result[count] = first_v;
               unselected.remove(first_v);
               count++;
               for (auto n : neighbors[first_v]) {
                  nb[n].remove(first_v);
               }
            }
            // path is maximal
            else break;
         }
      }
      return result;
   };

   auto max_cliques = [top](std::vector<int> degrees, FArray<GNSet> neighbors) {
      auto deg = degrees;
      auto nb = neighbors;
      std::vector<int> result(top);

      // create set of unselected vertices
      GNSet unselected = {};
      for(auto i : std::views::iota(0,top)) {
         unselected.insert(i);
      }

      std::vector<GNSet> cliques;

      while (!unselected.empty()) {
         auto it = std::max_element(std::begin(deg), std::end(deg));
         int idx = std::distance(std::begin(deg), it);
         unselected.remove(idx);
         for (auto n : neighbors[idx]) {
            nb[n].remove(idx);
         }
         GNSet candidates = nb[idx];
         GNSet clique = {idx};
         while (!candidates.empty()) {
            int idx2 = min(candidates);
            if (!unselected.contains(idx2)) {
               candidates.remove(idx2);
               continue;
            }
            // if n can maintain clique
            if ((clique - neighbors[idx2]).empty()) {
               clique.insert(idx2);
               unselected.remove(idx2);
               for (auto n : neighbors[idx2]) {
                  nb[n].remove(n);
               }
            }
            // if n cannot maintain clique
            else {
               candidates.remove(idx2);
            }
         }

         cliques.push_back(clique);
         for (int v : clique) {
            deg[v] = -top*top;
         }
      }

      GNSet prev_clique = {};
      int clique_idx = 0;
      // start with largest clique
      for (int i = 0; i < cliques.size(); i++) {
         if (cliques[i].size() > prev_clique.size()) {
            prev_clique = cliques[i];
            clique_idx = i;
         }
      }
      cliques.erase(cliques.begin() + clique_idx);

      int count = 0;
      for (int v : prev_clique) {
         result[count] = v;
         count++;
      }

      while(!cliques.empty()) {
         GNSet prev_adjacency = {};
         for (int v : prev_clique) {
            prev_adjacency = prev_adjacency | neighbors[v];
         }
         GNSet current_clique = {};
         int max_adjacency = -1;
         int clique_idx2 = 0;
         for (int i = 0; i < cliques.size(); i++) {
            int adjacency = (cliques[i] & prev_adjacency).size();
            if (adjacency > max_adjacency) {
               max_adjacency = adjacency;
               current_clique = cliques[i];
               clique_idx2 = i;
            }
         }

         cliques.erase(cliques.begin() + clique_idx2);
         for (int v : current_clique) {
            result[count] = v;
            count++;
         }
         prev_clique = current_clique;
      }

      return result;
   };

   // in the thesis this is used as the max degree ordering
   // simple max degree static ordering (no cliques)
   auto simple_max_deg = [top](std::vector<int> degrees, FArray<GNSet> neighbors) {
      auto deg = degrees;
      auto nb = neighbors;
      std::vector<int> result(top);

      // max deg heuristic
      for (int i = 0; i < top; i++) {
         auto it = std::max_element(std::begin(deg), std::end(deg));
         int idx = std::distance(std::begin(deg), it);
         result[i] = idx;
         deg[idx] = -top*top; // MIN_VALUE
         for (auto n : nb[idx]) {
            nb[n].remove(idx);
            deg[n]--;
         }
      }
      return result;
   };


   auto random_order = [top]() {
      std::random_device rd;
      std::mt19937 gen(0);

      std::vector<int> result(top);
      std::iota(result.begin(), result.end(), 0);

      std::shuffle (result.begin(), result.end(), gen);
      return result;
   };
   // ============== end of static variable orderings ============================================================

   auto neighbors_prev = instance.adj;
   std::vector<int> degrees(top);
   for (int i = 0; i < top; i++) {
      degrees[i] = neighbors_prev[i].size();
   }

   // CHOOSE WHICH STATIC ORDERING (5th arg)
   // !!WARNING!! solutions are not labelled correctly. Labels are shuffled due to variable ordering

   int ordering_selection = argc>=5 ? atoi(argv[4]) : 0;
   std::vector<int> sorted(top);
   switch (ordering_selection) {
   case 0:
      sorted = max_deg(degrees, neighbors_prev);
      break;
   case 1:
      sorted = min_width(degrees, neighbors_prev);
      break;
   case 2:
      sorted = max_connectivity(degrees, neighbors_prev);
      break;
   case 3:
      sorted = max_paths(neighbors_prev);
      break;
   case 4:
      sorted = max_cliques(degrees, neighbors_prev);
      break;
   case 5:
      sorted = random_order();
      break;
   case 6:
      sorted = min_deg_heuristic(degrees, neighbors_prev);
      break;
   case 7:
      sorted = simple_max_deg(degrees, neighbors_prev);
      break;
   default:
      sorted = max_deg(degrees, neighbors_prev);
   }

   // check for any errors in variable orderings
   std::cout<< sorted << std::endl;
   GNSet sorted_set = {};     //check duplicates
   for (int i : sorted) {
   sorted_set.insert(i);
   }
   std::cout<< sorted_set.size() << std::endl;
   for (int i = 0; i < top; i++) {
      if (!(sorted_set.contains(i))) {
         std::cout<< i<< std::endl;
         throw std::invalid_argument("ERROR ordering not sound");
      }
   }


   // CHOOSE BEAM SIZE (4th arg)
   int beam_n = argc>=4 ? atoi(argv[3]) : 4;

   std::vector<int> reverse_sorted(top); // for CHECKER bnds to work properly, reverse mapping from sorted to indices
   for (int i = 0; i < top; i++) {
      reverse_sorted[sorted[i]] = i;
   }

   // neighbors sorted by static ordering -> enforces the new ordering (along with the CHECKER changing)
   FArray<GNSet> neighbors(neighbors_prev.size());
   for (int i = 0; i < top; i++) {
      neighbors[i] = GNSet{};
      for (auto n : neighbors_prev[sorted[i]]) {
         neighbors[i].insert(reverse_sorted[n]);
      }
   }

   // now uses reverse_sorted to have a mapping to the original labels of the vertices
   Bounds bnds([&es, &reverse_sorted](const std::vector<int>& inc) {
      bool ok = true;
      for(const auto& e : es) {
         // std::find(v.begin(), v.end(), x) != v.end()
         bool v1In = (std::find(inc.begin(), inc.end(), reverse_sorted[e.a]) != inc.end()) ? true : false;
         bool v2In = (std::find(inc.begin(), inc.end(), reverse_sorted[e.b]) != inc.end()) ? true : false;
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

   std::vector<std::pair<int, int>> pq_empty;
   auto empty_state = MISPMDD { GNSet{}, pq_empty, -1};

   const auto myInit = [top,&neighbors]() {   // The root state
      GNSet U = {};
      std::vector<std::pair<int, int>> pq_init;

      for(auto i : std::views::iota(0,top)) {
         U.insert(i);
         pq_init.emplace_back(i, neighbors[i].size());
      }
      return MISPMDD{U, std::move(pq_init), 0};
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
      auto pq = s.pq;
      for (int i = 0; i < beam_n; i++) {
         auto p = pq.extract(pq.begin()).value();
         out.insert(p.first);
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
         if (label == top) {

            if (s.sel.size() <= beam_n) return empty_state;

            GNSet out = s.sel;
            GNSet removed = {};
            auto pq = s.pq;

            // remove beam
            for (int i = 0; i < beam_n; i++) {
               auto p = pq.extract(pq.begin()).value();
               out.remove(p.first);
               removed.insert(p.first);
            }

            GNSet adj_nb = {};
            for (int r : removed) {
               adj_nb = adj_nb | (neighbors[r] & out);
            }

            // change degrees neighbours of those deleted
            if (!adj_nb.empty()) {
               std::vector<std::pair<int, int>> temp(pq.begin(), pq.end());

               for (auto& elem : temp) {
                  if (adj_nb.contains(elem.first)) {
                     elem.second -= (neighbors[elem.first] & removed).size();
                  }
               }

               std::set<std::pair<int, int>, MISPMDD::Compare> new_pq(
                   std::make_move_iterator(temp.begin()),
                   std::make_move_iterator(temp.end())
               );
               pq = std::move(new_pq);
            }

            if (out.empty()) return empty_state;
            return MISPMDD { std::move(out), std::move(pq), s.transitions+1};
         }

         if (!s.sel.contains(label)) return std::nullopt; // we cannot take n (label==1) if not legal.

         GNSet out = s.sel;
         GNSet removed = {};
         auto pq = s.pq;

         // remove up to label
         for (int i = 0; i < beam_n; i++) {
            auto p = pq.extract(pq.begin()).value();
            out.remove(p.first);
            if (p.first == label) break;
            removed.insert(p.first);
         }
         // remove neighbors label vertex
         out.diffWith(neighbors[label]);
         removed = removed | (s.sel & neighbors[label]);

         for (auto it = pq.begin(); it != pq.end(); ) {
            if (neighbors[label].contains(it->first)) {
               it = pq.erase(it);
            } else {
               it++;
            }
         }

         GNSet adj_nb = {};
         for (int r : removed) {
            adj_nb = adj_nb | (neighbors[r] & out);   // was approx before
         }

         if (!adj_nb.empty()) {
               std::vector<std::pair<int, int>> temp(pq.begin(), pq.end());

               for (auto& elem : temp) {
                  if (adj_nb.contains(elem.first)) {
                     elem.second -= (neighbors[elem.first] & removed).size();
                  }
               }

               std::set<std::pair<int, int>, MISPMDD::Compare> new_pq(
                   std::make_move_iterator(temp.begin()),
                   std::make_move_iterator(temp.end())
               );
               pq = std::move(new_pq);
         }

         if (out.empty()) return empty_state;

         return MISPMDD { std::move(out), std::move(pq), s.transitions+1};
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

      std::set<std::pair<int, int>, MISPMDD::Compare> new_pq;
      for (auto [key, deg] : merged_map) {
         new_pq.emplace(key, deg);
      }

      return MISPMDD{
         un, std::vector<std::pair<int, int>>(new_pq.begin(), new_pq.end()), std::min(s1.transitions, s2.transitions)};
   };
   const auto eqs = [](const MISPMDD& s) -> bool {
      // if (s.sel.size() == 0) std::cout << s.l << std::endl;
      return s.sel.size() == 0;
   };
   const auto local = [&weight, &neighbors, &es, top](const MISPMDD& s, LocalContext) -> double {
      // return sum(s.sel,[&weight](auto v) { return weight[v];});

      int maxD = 0;
      int edg = 0;
      for (auto& elem : s.pq) {
         if (elem.second > maxD) maxD = elem.second;
         edg += elem.second;
      }
      edg = edg / 2.0;
      if (maxD == 0) return sum(s.sel,[&weight](auto v) { return weight[v];});
      else return s.sel.size() - std::floor(edg/maxD);
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

