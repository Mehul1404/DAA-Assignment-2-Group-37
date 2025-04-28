#include <bits/stdc++.h>
using namespace std;

// Implementation of Algorithm 1: Exact (flow-based binary search) for h-clique densest subgraph
// This is the pure Exact algorithm, not CoreExact or any core-based optimizations.
// Input: edge-list file, clique size h

class Graph {
public:
    int n;
    vector<vector<int>> adj;
    Graph(int _n=0): n(_n), adj(n) {}
    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    static Graph readEdgeList(const string &fname) {
        ifstream in(fname);
        if (!in) { cerr << "Cannot open " << fname << "\n"; exit(1); }
        int u, v, maxid = -1;
        vector<pair<int,int>> edges;
        string line;
        while (getline(in, line)) {
            if (line.empty() || line[0]=='#') continue;
            istringstream iss(line);
            if (!(iss>>u>>v)) continue;
            edges.emplace_back(u,v);
            maxid = max(maxid, max(u,v));
        }
        Graph G(maxid+1);
        for (auto &e: edges) G.addEdge(e.first, e.second);
        return G;
    }
};

// Dinic max-flow implementation
struct Edge { int to; long long cap; int rev; };
struct Dinic {
    int N, s, t;
    vector<vector<Edge>> g;
    vector<int> level, ptr;
    Dinic(int _N, int _s, int _t): N(_N), s(_s), t(_t), g(N), level(N), ptr(N) {}
    void addEdge(int u, int v, long long c) {
        g[u].push_back({v, c, (int)g[v].size()});
        g[v].push_back({u, 0, (int)g[u].size()-1});
    }
    bool bfs() {
        fill(level.begin(), level.end(), -1);
        queue<int>q; q.push(s); level[s] = 0;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (auto &e: g[u]) if (e.cap>0 && level[e.to]<0) {
                level[e.to] = level[u] + 1;
                q.push(e.to);
            }
        }
        return level[t]>=0;
    }
    long long dfs(int u, long long f) {
        if (!f || u==t) return f;
        for (int &cid = ptr[u]; cid < (int)g[u].size(); ++cid) {
            auto &e = g[u][cid];
            if (e.cap>0 && level[e.to]==level[u]+1) {
                long long pushed = dfs(e.to, min(f, e.cap));
                if (pushed) {
                    e.cap -= pushed;
                    g[e.to][e.rev].cap += pushed;
                    return pushed;
                }
            }
        }
        return 0;
    }
    long long maxflow() {
        long long flow = 0;
        while (bfs()) {
            fill(ptr.begin(), ptr.end(), 0);
            while (long long pushed = dfs(s, LLONG_MAX))
                flow += pushed;
        }
        return flow;
    }
};

// Enumerate all (h-1)-cliques in G, return list of vertex lists
void enumCliques(const Graph &G, int h1, vector<vector<int>> &out) {
    int n = G.n;
    vector<int> clique;
    function<void(vector<int>&)> dfs = [&](vector<int> &cand) {
        if ((int)clique.size() == h1) {
            out.push_back(clique);
            return;
        }
        while (!cand.empty()) {
            int v = cand.back(); cand.pop_back();
            vector<int> nxt;
            for (int u: cand) {
                if (find(G.adj[v].begin(), G.adj[v].end(), u) != G.adj[v].end())
                    nxt.push_back(u);
            }
            clique.push_back(v);
            dfs(nxt);
            clique.pop_back();
        }
    };
    for (int v = 0; v < n; ++v) {
        vector<int> cand;
        for (int u: G.adj[v]) if (u > v) cand.push_back(u);
        clique = {v};
        dfs(cand);
    }
}

int main(int argc, char*argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <edge-file> <h>\n";
        return 1;
    }

    // Input
    Graph G = Graph::readEdgeList(argv[1]);
    int h = stoi(argv[2]);
    if (h < 2) {
        cerr << "h must be >= 2\n";
        return 1;
    }
    int n = G.n;

    // Precompute all (h-1)-cliques
    vector<vector<int>> hm1;
    enumCliques(G, h-1, hm1);
    int H = hm1.size();

    // Binary search on density alpha
    double l = 0, u = 0;
    // upper bound = max clique-degree
    vector<long long> cliqueDeg(n,0);
    for (int i = 0; i < H; ++i) {
        for (int v : hm1[i]) cliqueDeg[v]++;
    }
    u = *max_element(cliqueDeg.begin(), cliqueDeg.end());

    pair<vector<int>, double> best = {{}, 0.0};
    double eps = 1.0 / (n*(n-1));
    while (u - l > eps) {
        double alpha = (l + u) / 2;

        // build flow network
        int N = 1 + n + H + 1;
        int s = 0, t = N-1;
        Dinic din(N, s, t);
        // s->v with cap cliqueDeg[v]
        for (int v = 0; v < n; ++v) din.addEdge(s, 1+v, cliqueDeg[v]);
        // v->t with cap alpha*h
        for (int v = 0; v < n; ++v) din.addEdge(1+v, 1+n+H, (long long)(alpha * h));
        // (h-1)-clique nodes: index 1+n ... 1+n+H-1
        for (int i = 0; i < H; ++i) {
            int cid = 1 + n + i;
            for (int v: hm1[i]) din.addEdge(cid, 1+v, LLONG_MAX);
        }
        for (int v = 0; v < n; ++v) {
            for (int i = 0; i < H; ++i) {
                // check if v completes hm1[i]
                bool inClique = false;
                for (int u2: hm1[i]) if (u2 == v) { inClique = true; break; }
                if (inClique) continue;
                // v adjacent to all in hm1[i]?
                bool ok = true;
                for (int u2: hm1[i]) {
                    if (find(G.adj[v].begin(), G.adj[v].end(), u2) == G.adj[v].end()) { ok = false; break; }
                }
                if (ok) din.addEdge(1+v, 1+n+i, 1);
            }
        }
        // compute flow and min-cut
        long long flow = din.maxflow();
        // total cap from source = sum cliqueDeg[v]
        long long total = accumulate(cliqueDeg.begin(), cliqueDeg.end(), 0LL);

        if (flow < total) {
            // feasible: extract S side
            vector<char> vis(N,0);
            queue<int>q; vis[s]=1; q.push(s);
            while(!q.empty()){
                int x=q.front(); q.pop();
                for(auto &e: din.g[x]) if (e.cap>0 && !vis[e.to]) {
                    vis[e.to]=1; q.push(e.to);
                }
            }
            vector<int> sub;
            for (int v = 0; v < n; ++v) if (vis[1+v]) sub.push_back(v);
            best = {sub, alpha};
            l = alpha;
        } else {
            u = alpha;
        }
    }

                // Output the vertex list of the densest subgraph
     cout << "Densest subgraph size = " << best.first.size() <<endl<< " Vertices:";
    for (int v : best.first) cout << " " << v;
    cout <<endl;

    int sz = best.first.size();
    // Compute edge-density of the returned subgraph
    unordered_set<int> inSub(best.first.begin(), best.first.end());
    long long edgeCount = 0;
    for (int v : best.first) {
        for (int u : G.adj[v]) {
            if (inSub.count(u) && u > v) {
                edgeCount++;
            }
        }
    }
    double edgeDensity = (double)edgeCount / sz;
    cout << "Edge-density = " << edgeDensity << endl;

    // Compute h-clique-density of the returned subgraph
    // Build induced subgraph S with remapped vertices
    Graph S(sz);
    unordered_map<int,int> rmap;
    for (int i = 0; i < sz; ++i) rmap[best.first[i]] = i;
    for (int v : best.first) {
        for (int u : G.adj[v]) {
            if (inSub.count(u) && rmap[u] > rmap[v]) {
                S.addEdge(rmap[v], rmap[u]);
            }
        }
    }
    vector<vector<int>> subCliques;
    enumCliques(S, h, subCliques);
    long long cCount = subCliques.size();
    double cliqueDensity = (double)cCount / sz;
    cout << h << "-clique-density = " << cliqueDensity
         << " (" << cCount << " cliques over " << sz << " vertices) "<<endl;

    return 0;
}
