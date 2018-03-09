// Read a mutation probability and gene data source file name from the
// command line.  Initialize a bitvector population from the data
// stream, and calculate the relation distance between all pairs of
// vectors in the population.
//
// Find an undirected graph that spans all the bitvectors and
// minimizes the bit difference between adjacent bitvectors normalized
// to the number of expected mutations.
//
// Orient the graph into a rooted tree by first transforming the
// spanning graph's vector and edge representation into a neighborhood
// representation.  The leaves of the tree are the bitvectors with
// only a single neighbor.  Find the leaves, note their parents, then
// trim them from the graph, and repeat until there is a single
// bitvector left, the "progenitor".
//
// Run 'make test' to build the executables, run the tests, and check
// the results.  This code was developed on MacOSX, but should run on
// any standardly-endowed Unix system with a Makefile tweak or two.
//
// Or just run 'make bitvectors-parents.data' for the "solution".

#include <algorithm>
#include <bitset>
#include <cassert>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;


static void showUsage(ostream &errs, const char *cmd)
{
  errs << endl
       << "Usage: " << cmd << " <prob> <data>" << endl
       << "Where: <prop> is the bitwise probability of mutation " << endl
       << "              as an integer percentage (20 for example). " << endl
       << "       <data> is a file of " << SCALE
       << " bit strings of length " << SCALE << "." << endl
       << "Each line matches the regular expression "
       << "'^[01]{" << SCALE << "}$', " << endl
       << "and there are " << SCALE << " lines in <data>." << endl
       << endl;
}


// A vector of SCALE bits with mutationPercentage, the bitwise
// probability of mutation each generation.
//
// The difference of two BitVectors is the number of bits different
// between lhs and rhs normalized to the expected mutation count.
//
struct BitVector
{
  static int mutationPercentage;
  int index;
  bitset<SCALE> bits;

  friend size_t operator-(const BitVector &lhs, const BitVector &rhs)
  {
    static const size_t expected = (SCALE * mutationPercentage / 100);
    bitset<SCALE> difference(lhs.bits ^ rhs.bits);
    const size_t distance = difference.count();
    if (distance > expected) return distance - expected;
    return expected - distance;
  }

  BitVector(): index(-1), bits() {}
  BitVector(int n, const string &s): index(n), bits(s) {}
};


// A population of BitVectors initialized from the stream s.
// There is a problem with error on line when *this is false.
//
struct Population
{
  vector<BitVector> bitVectors;
  int line;
  string error;

  operator bool() { return line == -1; }

  Population(istream &s): bitVectors(SCALE), line(-1), error()
  {
    string bitvector;
    int n = 0;
    for (; s && n < bitVectors.size(); ++n) {
      if (getline(s, bitvector) && bitvector.size() == SCALE) {
        bitVectors[n] = BitVector(n, bitvector);
      } else {
        break;
      }
    }
    if (n < bitVectors.size()) {
      line = n; error = bitvector;
    }
  }
};


// Relation between two BitVectors, with indexes left and right, and
// the normalized bit distance nbd between them.
//
struct Relation
{
  size_t nbd;
  int left; int right;

  // True if lhs is probably a closer relation than rhs.
  //
  friend bool operator<(const Relation &lhs, const Relation &rhs)
  {
    return lhs.nbd < rhs.nbd;
  }

  Relation(): nbd(0), left(0), right(0) {}
  Relation(size_t d, int p, int c): nbd(d), left(p), right(c) {}
};


// A connected graph of Relations built while constructing a
// SpanningGraph.
//
struct ConnectedGraph
{
  set<int> vertexes;
  vector<Relation> edges;

  // True if this graph contains all of the bitvectors.
  //
  bool full() { return vertexes.size() >= SCALE; }

  // Add e to this.
  //
  void add(const Relation &e) {
    vertexes.insert(e.left);
    vertexes.insert(e.right);
    edges.push_back(e);
  }

  // True if e connects to (has some vertex in common with) this.
  // Otherwise false.
  //
  bool connectsTo(const Relation &e) {
    return vertexes.count(e.left) || vertexes.count(e.right);
  }

  // Merge that with this (when some new Edge connects two subgraphs).
  //
  void mergeWith(const ConnectedGraph &that) {
    vertexes.insert(that.vertexes.begin(), that.vertexes.end());
    edges.insert(edges.end(), that.edges.begin(), that.edges.end());
  }

  ConnectedGraph(): vertexes(), edges() {}
  ConnectedGraph(const Relation &e): vertexes(), edges() {
    add(e);
  }
};


// A undirected graph spanning all bitvectors and minimizing the
// Relation distance between connected bitvectors.
//
struct SpanningGraph
{
  typedef list<ConnectedGraph>::iterator Cp;
  ConnectedGraph result;

  // Merge any two subgraphs in connected joined by r.
  //
  Cp merge(list<ConnectedGraph> &connected, Cp resultP, const Relation &r) {
    for (Cp pC = connected.begin(); pC != connected.end(); ++pC) {
      if (pC->connectsTo(r)) {
        if (resultP == connected.end()) {
          resultP = pC;
        } else {
          resultP->mergeWith(*pC);
          connected.erase(pC);
          break;
        }
      }
    }
    return resultP;
  }

  // Add r to some subgraph.  Start a new subgraph if no connection
  // was already found by merge().  Return true if graph at resultP
  // spans all bitvectors.  Return false if the graph extraction
  // should continue.
  //
  bool add(list<ConnectedGraph> &connected, Cp resultP, const Relation &r) {
    if (resultP == connected.end()) {
      ConnectedGraph cg(r);
      connected.push_back(cg);
    } else {
      resultP->add(r);
      if (resultP->full()) return true;
    }
    return false;
  }

  // Return all the relations in population.
  //
  static vector<Relation> findAll(const Population &population)
  {
    const vector<BitVector> &bitVectors = population.bitVectors;
    const size_t size = bitVectors.size();
    typedef vector<BitVector>::const_iterator BvP;
    vector<Relation> result(size * (size - 1) / 2);
    vector<Relation>::iterator op = result.begin();
    for (BvP lp = bitVectors.begin(); lp != bitVectors.end(); ++lp) {
      for (BvP rp = lp + 1; rp < bitVectors.end(); ++op, ++rp) {
        *op = Relation(*lp - *rp, lp->index, rp->index);
      }
    }
    return result;
  }

  // *this is true when the result graph spans all the bitvectors.
  // Otherwise false.
  //
  operator bool() { return result.vertexes.size() == SCALE; }

  // Find a graph spanning all the bitvectors in relations that
  // minimizes the normalized bit distances between bitvectors.
  //
  SpanningGraph(const Population &population): result()
  {
    typedef vector<Relation>::const_iterator Rp;
    vector<Relation> relations(findAll(population));
    sort(relations.begin(), relations.end());
    list<ConnectedGraph> connected;
    Cp resultP = connected.end();
    for (Rp pR = relations.begin(); pR != relations.end(); ++pR) {
      resultP = merge(connected, resultP, *pR);
      if (add(connected, resultP, *pR)) break;
      resultP = connected.end();
    }
    if (resultP != connected.end()) result = *resultP;
  }
};


// From the undirected connected graph cg, extract a directed tree,
// where result is a 0-based array such that result[n] is the parent
// of child n and result[r] == -1 means that r is the root of the
// tree.
//
struct Genealogy {
  typedef map<int, set<int> > Neighbors;
  typedef Neighbors::iterator Np;
  typedef vector<Np> Leaves;
  typedef Leaves::iterator Pl;
  vector<int> result;
  bool itsOk;

  // Return a map of bitvector indexes to the indexes of their
  // neighbors from the Relations in edges.
  //
  static Neighbors discover(const vector<Relation> &edges) {
    typedef vector<Relation>::const_iterator Rp;
    Neighbors result;
    for (Rp pE = edges.begin(); pE != edges.end(); ++pE) {
      result[pE->left].insert(pE->right);
      result[pE->right].insert(pE->left);
    }
    return result;
  }

  // Remove all leaves from the map neighbors.  Return true if a leaf
  // was removed from neighbors.  Otherwise return false.
  //
  static bool trim(Neighbors &neighbors, Leaves &leaves) {
    bool result = false;
    for (Pl pL = leaves.begin(); pL != leaves.end(); ++pL) {
      neighbors.erase(*pL);
      result = true;
    }
    return result;
  }

  friend ostream &operator<<(ostream &s, const Genealogy &g) {
    typedef vector<int>::const_iterator Pi;
    for (Pi pI = g.result.begin(); pI != g.result.end(); ++pI) {
      s << *pI << endl;
    }
    return s;
  }

  // True if this converged to a rooted tree.
  //
  operator bool() { return itsOk; }

  // Extract neighborhoods from cg.  Find the leaves (vertexes with
  // only one neighbor).  For each leaf, note its parent in result,
  // remove it from its parent's neighborhood, then erase it from
  // neighbors.  Consume neighbors until only one vertex is left.
  //
  Genealogy(const ConnectedGraph &cg): result(SCALE, -1)
  {
    Neighbors neighbors(discover(cg.edges));
    while (neighbors.size() > 1) {
      vector<Np> leaves;
      for (Np pC = neighbors.begin(); pC != neighbors.end(); ++pC) {
        if (pC->second.size() == 1) {
          const size_t child = pC->first;
          const size_t parent = *pC->second.begin();
          result[child] = parent;
          pC->second.erase(pC->second.begin());
          leaves.push_back(pC);
          neighbors.find(parent)->second.erase(child);
        }
      }
      itsOk = trim(neighbors, leaves);
      if (!itsOk) break;
    }
  }
};


int BitVector::mutationPercentage;
static bool initializeMutationPercentage(char *percentageString)
{
  istringstream issMp; issMp.str(percentageString);
  return issMp >> BitVector::mutationPercentage
    && BitVector::mutationPercentage >= 0
    && BitVector::mutationPercentage <= 100;
}

int main(int ac, char *av[])
{
  if (ac == 3) {
    if (initializeMutationPercentage(av[1])) {
      ifstream dataStream(av[2]);
      Population population(dataStream);
      if (population) {
        SpanningGraph graph(population);
        if (graph) {
          Genealogy genealogy(graph.result);
          if (genealogy) {
            cout << genealogy;
            return 0;
          } else {
            cerr << av[0] << ": Error: The genealogy did not converge." << endl;
          }
        } else {
          cerr << av[0] << ": Error: Cannot relate entire population." << endl;
        }
      } else {
        cerr << av[0] << ": Error on line " << population.line
             << ": " << population.error << endl;
      }
    } else {
      cerr << av[0] << ": Error: First argument '" << av[1]
           << "' should be an integer between 0 and 100." << endl;
    }
  }
  showUsage(cerr, av[0]);
  return 1;
}
