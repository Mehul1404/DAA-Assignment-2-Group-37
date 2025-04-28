#include <bits/stdc++.h>
using namespace std;

// CoreExact: exact algorithm for h-clique densest subgraph (CDS)
// Implements Algorithm 3 (k-clique-core decomposition) and Algorithm 4 (CoreExact) with Algorithm 1 flow-based search

class Graph {
public:
    int n;
    vector<vector<int>> adj;
    Graph(int _n=0): n(_n), adj(n) {}
    void addEdge(int u, int v) {
        if (u<0||v<0||u>=n||v>=n) return;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    static Graph readEdgeList(const string &fname) {
        ifstream in(fname);
        if (!in) { cerr<<"Cannot open "<<fname<<"\n"; exit(1); }
        vector<pair<int,int>> edges;
        string line;
        int u,v, maxid = -1;
        while (getline(in,line)) {
            if (line.empty()||line[0]=='#') continue;
            istringstream iss(line);
            if (!(iss>>u>>v)) continue;
            edges.emplace_back(u,v);
            maxid = max(maxid, max(u,v));
        }
        Graph G(maxid+1);
        for (auto &e:edges) G.addEdge(e.first,e.second);
        return G;
    }
};

// Enumerate all k-cliques in graph G using recursive backtracking
void enumCliques(const Graph &G, int k, vector<vector<int>> &out) {
    int n = G.n;
    vector<int> clique;
    function<void(vector<int>&)> dfs = [&](vector<int> &cand) {
        if ((int)clique.size() == k) {
            out.push_back(clique);
            return;
        }
        while (!cand.empty()) {
            int v = cand.back(); cand.pop_back();
            vector<int> next;
            for (int u: cand) {
                // check adjacency
                if (find(G.adj[v].begin(), G.adj[v].end(), u) != G.adj[v].end())
                    next.push_back(u);
            }
            clique.push_back(v);
            dfs(next);
            clique.pop_back();
        }
    };
    for (int v=0; v<n; ++v) {
        vector<int> cand;
        for (int u: G.adj[v]) if (u>v) cand.push_back(u);
        clique = {v};
        dfs(cand);
    }
}

// Algorithm 3: compute (k,Ψ)-core numbers via h-clique enumeration
vector<int> computeCliqueCore(const Graph &G, int h) {
    int n = G.n;
    vector<vector<int>> cliques;
    enumCliques(G, h, cliques);
    int m = cliques.size();
    vector<int> deg(n,0);
    vector<vector<int>> inc(n);
    for (int i=0;i<m;++i) {
        for (int v: cliques[i]) {
            deg[v]++;
            inc[v].push_back(i);
        }
    }
    int maxDeg = *max_element(deg.begin(), deg.end());
    vector<vector<int>> bucket(maxDeg+1);
    for (int v=0;v<n;++v) bucket[deg[v]].push_back(v);
    vector<int> core(n,0);
    vector<bool> removed(n,false);
    int rem=0;
    for (int d=0;d<=maxDeg && rem<n;++d) {
        for (int idx=0; idx<bucket[d].size(); ++idx) {
            int v = bucket[d][idx];
            if (removed[v]||deg[v]!=d) continue;
            removed[v]=true;
            core[v]=d;
            rem++;
            for (int ci: inc[v]) for (int u: cliques[ci]) {
                if (!removed[u] && deg[u]>d) {
                    int du = deg[u];
                    deg[u]--;
                    bucket[du-1].push_back(u);
                }
            }
        }
    }
    return core;
}

// Dinic max-flow
struct Edge{int to; long long cap; int rev;};
struct Dinic{
    int N,s,t;
    vector<vector<Edge>> g;
    vector<int> level, ptr;
    Dinic(int _N,int _s,int _t):N(_N),s(_s),t(_t),g(N),level(N),ptr(N){}
    void addEdge(int u,int v,long long c) {
        g[u].push_back({v,c,(int)g[v].size()});
        g[v].push_back({u,0,(int)g[u].size()-1});
    }
    bool bfs() {
        fill(level.begin(),level.end(),-1);
        queue<int>q; q.push(s); level[s]=0;
        while(!q.empty()){
            int u=q.front();q.pop();
            for(auto &e:g[u]) if(e.cap>0 && level[e.to]<0) {
                level[e.to]=level[u]+1; q.push(e.to);
            }
        }
        return level[t]>=0;
    }
    long long dfs(int u,long long f) {
        if (!f||u==t) return f;
        for (int &cid=ptr[u]; cid<g[u].size(); ++cid) {
            auto &e = g[u][cid];
            if (e.cap>0 && level[e.to]==level[u]+1) {
                long long pushed = dfs(e.to, min(f,e.cap));
                if (pushed) { e.cap-=pushed; g[e.to][e.rev].cap+=pushed; return pushed; }
            }
        }
        return 0;
    }
    long long maxflow() {
        long long flow=0;
        while(bfs()){ fill(ptr.begin(),ptr.end(),0);
            while(long long pushed=dfs(s,LLONG_MAX)) flow+=pushed;
        }
        return flow;
    }
};

int main(int argc,char*argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    if (argc!=3) { cerr<<"Usage: "<<argv[0]<<" <edge-file> <h>\n"; return 1; }
    Graph G = Graph::readEdgeList(argv[1]);
    int h = stoi(argv[2]); if (h<2) { cerr<<"h>=2 required\n"; return 1; }

    // 1. Compute core numbers
    auto coreNum = computeCliqueCore(G, h);
    int kmax = *max_element(coreNum.begin(), coreNum.end());

    // 2. Extract (kmax,Ψ)-core vertices
    vector<int> inv;
    vector<int> id(G.n, -1);
    for (int v=0; v<G.n; ++v) if (coreNum[v]>=kmax) { id[v]=inv.size(); inv.push_back(v); }
    int cN = inv.size();
    Graph CG(cN);
    for (int i=0;i<cN;++i) {
        for (int u: G.adj[inv[i]]) if (id[u]>=0 && id[u]>i) CG.addEdge(i,id[u]);
    }

    // 3. Enumerate (h-1)-cliques in CG
    vector<vector<int>> cliques_hm1;
    enumCliques(CG, h-1, cliques_hm1);
    int C = cliques_hm1.size();

    // 4. Binary search
    double l = (double)kmax / h, u = kmax;
    vector<int> best = inv;
    while (u - l > 1.0/(cN*(cN-1))) {
        double alpha = (l+u)/2;
        int Nn = 1 + cN + C + 1;
        int s=0, t=Nn-1;
        Dinic din(Nn,s,t);
        vector<long long> degCD(cN,0);
        long long total=0;
        // compute clique-degree in CG
        for (int i=0;i<C;++i) for (int v: cliques_hm1[i]) degCD[v]++;
        for (int v=0;v<cN;++v) { din.addEdge(s,1+v,degCD[v]); total+=degCD[v]; }
        // v->t
        for (int v=0;v<cN;++v) din.addEdge(1+v,1+cN+C,(long long)(alpha*h));
        // clique->verts
        for (int i=0;i<C;++i) {
            int cid = 1+cN+i;
            for (int v: cliques_hm1[i]) din.addEdge(cid,1+v,LLONG_MAX);
        }
        // verts->clique
        for (int i=0;i<C;++i) {
            int cid = 1+cN+i;
            for (int v=0;v<cN;++v) {
                bool inC=false;
                for (int u2: cliques_hm1[i]) if (u2==v) { inC=true; break; }
                if (inC) continue;
                bool ok=true;
                for (int u2: cliques_hm1[i]) {
                    if (find(CG.adj[v].begin(),CG.adj[v].end(),u2)==CG.adj[v].end()) { ok=false; break; }
                }
                if (ok) din.addEdge(1+v,cid,1);
            }
        }
        long long flow = din.maxflow();
        if (flow < total) {
            // extract reachable
            vector<char> vis(Nn,0);
            queue<int>q; vis[s]=1; q.push(s);
            while(!q.empty()){
                int x=q.front();q.pop();
                for (auto &e: din.g[x]) if (e.cap>0 && !vis[e.to]) { vis[e.to]=1; q.push(e.to); }
            }
            vector<int> cand;
            for (int v=0;v<cN;++v) if (vis[1+v]) cand.push_back(inv[v]);
            best=cand;
            l=alpha;
        } else {
            u=alpha;
        }
    }

    // 5. Output
    cout << "Densest subgraph size = " << best.size() <<endl<< " Vertices:";
    for (int v : best) cout << " " << v;
    cout <<endl;

    // Compute and print densities for verification
    // Build induced subgraph S on 'best'
    int sz = best.size();
    unordered_set<int> inBest(best.begin(), best.end());
    Graph S(sz);
    // map old to new indices
    unordered_map<int,int> rmap;
    for (int i = 0; i < sz; ++i) rmap[best[i]] = i;
    long long edgeCount = 0;
    for (int old_v : best) {
        int v = rmap[old_v];
        for (int nbr : G.adj[old_v]) {
            if (inBest.count(nbr)) {
                int u = rmap[nbr];
                if (u > v) {
                    S.addEdge(v, u);
                    edgeCount++;
                }
            }
        }
    }
    double edgeDensity = (double)edgeCount / sz;
    cout << "Edge-density = " << edgeDensity <<endl;

    // Count h-cliques in S
    vector<vector<int>> clqs;
    enumCliques(S, h, clqs);
    long long mu = clqs.size();
    double cliqueDensity = (double)mu / sz;
    cout << h << "-clique-density = " << cliqueDensity << " (" << mu << " cliques over " << sz << " vertices)"<<endl;

    return 0;
}
