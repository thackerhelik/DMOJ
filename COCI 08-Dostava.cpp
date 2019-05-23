/*
Make a graph of the first and last columns and run dijkstra for multiple source (note floyd-warshall will take n3 whereas dijkstra will take n(nlogm + mlogn))
Then for each source to destination count the limited possibilities in O(1) using prefix sums for inter-row and for one row to other we have calculated by using dijkstra
*/

/*
When you walk through a storm
Hold your head up high
And don't be afraid of the dark
At the end of the storm
There's a golden sky
And the sweet silver song
of the lark
Walk on, through the wind
Walk on, through the rain
Though your dreams be tossed
and blown
Walk on, walk on
With hope in your heart
And you'll never walk alone
YNWA
*/

//hell_hacker
//IT TAKES EVERYTHING and IT IS NOT LUCK
//PICK ME LAST -- COME OUT OF N^WHERE
//WHY NOT?

/*
And you came my way on a winner's day
Shouted loudly come out and play
Can't you tell I got news for you
Sun is shining and so are you
*/

#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <bits/stdc++.h>

using namespace __gnu_pbds;
using namespace std;

#define endl '\n'
#define fo(i,n) for(i=0;i<n;++i)
#define forr(i,n) for(i=n-1;i>=0;--i)
#define rep(i,a,b) for(i=a;i<=b;++i)
#define per(i,a,b) for(i=b;i>=a;--i)

#define ffs(a) __builtin_ffs(a) // find first set
#define clz(a) __builtin_clz(a) // count leading zeroes
#define ctz(a) __builtin_ctz(a) // count trailing zeroes
#define popc(a) __builtin_popcount(a) // count set bits
#define popcll(a) __builtin_popcountll(a) //count set bits for long long int

#define all(x) (x).begin(),(x).end()
#define sz(x) (int)(x).size()
#define fi first
#define se second
#define mp make_pair
#define pb push_back
#define yolo "debug_statement"

typedef long long int ll;
typedef long double ld;
typedef pair<int,int> ii;
typedef vector<ii> vii;
typedef vector<int> vi;

const ll inf = 1e9;
const ld eps = 1e-9;
const ld pi=acos(-1);
const ll mod=1000000007;

ll powmod(ll a,ll b,ll mo=mod){ll res=1;a%=mo; assert(b>=0); for(;b;b>>=1){if(b&1)res=res*a%mo;a=a*a%mo;}return res;}

inline ll inv_(ll a) {return powmod(a, mod-2, mod);}
inline ll add(ll a, ll b){a%=mod; b%=mod; a+=b;if(a>=mod)a-=mod;return a;}
inline ll mul(ll a, ll b){a%=mod; b%=mod; return a*1ll*b%mod;}
inline ll sub(ll a, ll b){a%=mod; b%=mod; a-=b;if(a<0)a+=mod;return a;}

typedef tree<
int,
null_type,
less<int>,
rb_tree_tag,
tree_order_statistics_node_update>
ordered_set;

/*USAGE
ordered_set X;
X.insert(1);
X.insert(3);
cout << *X.find_by_order(1) << endl; //3
cout << X.order_of_key(-5) << endl; //0
cout << X.order_of_key(1) << endl; //0
cout << (end(X) == X.order_of_key(100)) << endl; //true since utne elements nahi hai
*/

#if DEBUG && !ONLINE_JUDGE

    #define debug(args...)     (Debugger()) , args

    class Debugger
    {
        public:
        Debugger(const std::string& _separator = ", ") :
        first(true), separator(_separator){}

        template<typename ObjectType>
        Debugger& operator , (const ObjectType& v)
        {
            if(!first)
                std::cerr << separator;
            std::cerr << v;
            first = false;
            return *this;
        }
        ~Debugger() {  std::cerr << endl;}

        private:
        bool first;
        std::string separator;
    };

    template <typename T1, typename T2>
    inline std::ostream& operator << (std::ostream& os, const std::pair<T1, T2>& p)
    {
        return os << "(" << p.first << ", " << p.second << ")";
    }

    template <typename T1, unsigned long int N>
    inline std::ostream& operator << (std::ostream& os, const std::array<T1, N>& a){
        bool first = true;
        os << "[";
        for(unsigned int i = 0; i < N; i++)
        {
            if(!first)
                os << ", ";
            os << a[i];
            first = false;
        }
        return os << "]";   
    }

    template<typename T>
    inline std::ostream &operator << (std::ostream & os,const std::vector<T>& v)
    {
        bool first = true;
        os << "[";
        for(unsigned int i = 0; i < v.size(); i++)
        {
            if(!first)
                os << ", ";
            os << v[i];
            first = false;
        }
        return os << "]";
    }

    template<typename T>
    inline std::ostream &operator << (std::ostream & os,const std::set<T>& v)
    {
        bool first = true;
        os << "[";
        for (typename std::set<T>::const_iterator iii = v.begin(); iii != v.end(); ++iii)
        {
            if(!first)
                os << ", ";
            os << *iii;
            first = false;
        }
        return os << "]";
    }

    template<typename T1, typename T2>
    inline std::ostream &operator << (std::ostream & os,const std::map<T1, T2>& v)
    {
        bool first = true;
        os << "[";
        for (typename std::map<T1, T2>::const_iterator iii = v.begin(); iii != v.end(); ++iii)
        {
            if(!first)
                os << ", ";
            os << *iii ;
            first = false;
        }
        return os << "]";
    }

#else
    #define debug(args...)                  // Just strip off all debug tokens
#endif

const int naxR = 2005;
const int naxC = 205;

int R, C, D;
int matrix[naxR][naxC];
int prefmatrix[naxR][naxC];
int dist[2*naxR][2*naxR];

vector<pair<int,int>>adj[2*naxR];
vector<int>d;

void dijkstra(int x, int y, int s) {
    d.assign(2*naxR, inf);

    d[s] = matrix[x][y];
    using pii = pair<int, int>;
    priority_queue<pii, vector<pii>, greater<pii>> q;
    q.push({d[s], s});
    while (!q.empty()) {
        int v = q.top().second;
        int d_v = q.top().first;
        q.pop();
        if (d_v != d[v])
            continue;

        for (auto edge : adj[v]) {
            int to = edge.first;
            int len;
            if(to % 2 == 1)
            	len = edge.second + matrix[to/2 + 1][1];
            else
            	len = edge.second + matrix[to/2][C];
            if (d[v] + len < d[to]){
                d[to] = d[v] + len;
                q.push({d[to], to});
            }
        }
    }
    for(int j = 0; j < sz(d); ++j){
    	dist[s][j] = d[j];
    }
}

int min4(int a, int b, int c, int d){
	return min(a, min(b, min(c, d)));
}

int main(){
    #if ONLINE_JUDGE
        ios_base::sync_with_stdio(0);cin.tie(NULL);cout.tie(NULL);
	#endif

    int a, b;

    cin >> R >> C;        

    for(int i = 1; i <= R; ++i){
    	int sum = 0;
    	for(int j = 1; j <= C; ++j){
    		cin >> matrix[i][j];
    		sum = sum + matrix[i][j];
    		if(j > 1)
    			prefmatrix[i][j] = prefmatrix[i][j - 1] + matrix[i][j];
    		else
    			prefmatrix[i][j] = matrix[i][j];
    	}
    	adj[i*2 - 1].push_back({i*2, sum - matrix[i][1] - matrix[i][C]});
    	adj[i*2].push_back({i*2 - 1, sum - matrix[i][1] - matrix[i][C]});
    }

    for(int i = 2; i <= R; ++i){
    	adj[i*2 - 1].push_back({(i - 1)*2 - 1, 0});
    	adj[(i - 1)*2 - 1].push_back({i*2 - 1, 0});
    	adj[i*2].push_back({(i - 1)*2, 0});
    	adj[(i - 1)*2].push_back({i*2, 0});
    }

    for(int i = 1; i <= R; ++i){
    	dijkstra(i, 1, i*2 - 1);
    	dijkstra(i, C, i*2);
    }

    cin >> D;

    int preva = 1, prevb = 1;
    ll sum = 0;

    for(int i = 0; i < D; ++i){
    	cin >> a >> b;
    	if(i < (D - 1))
    		sum = sum - matrix[a][b];
    	if(preva == a){
    		//same row
    		int fx, fy, sx, sy, ans1, ans2;
    		fx = a*2 - 1; fy = 1; sx = a*2; sy = 1;		
    		if(b >= prevb)
    			ans1 = dist[2*a - 1][2*a] + prefmatrix[preva][prevb] - matrix[preva][1] + prefmatrix[a][C - 1] - prefmatrix[a][b - 1];
    		else
    			ans1 = dist[2*a - 1][2*a] + prefmatrix[a][b] - matrix[a][1] + prefmatrix[preva][C - 1] - prefmatrix[preva][prevb - 1];
    		if(b >= prevb)
    			ans2 = prefmatrix[a][b] - prefmatrix[a][prevb - 1];
    		else
    			ans2 = prefmatrix[a][prevb] - prefmatrix[a][b - 1];
    		// cout << min(ans1, ans2) << endl;
    		sum = sum + min(ans1, ans2);
    	}
    	else{
    		//different row
    		int fx, fy, sx, sy, ans1, ans2, ans3, ans4;
    		fx = preva*2 - 1; fy = 1; sx = a*2 - 1; sy = 1;
    		ans1 = dist[2*preva - 1][2*a - 1] + prefmatrix[preva][prevb] - matrix[preva][1] + prefmatrix[a][b] - matrix[a][1];
    		fx = preva*2 - 1; fy = 1; sx = a*2; sy = C;
    		ans2 = dist[2*preva - 1][2*a] + prefmatrix[preva][prevb] - matrix[preva][1] + prefmatrix[a][C - 1] - prefmatrix[a][b - 1];
    		fx = preva*2; fy = C; sx = a*2 - 1; sy = 1;
    		ans3 = dist[2*preva][2*a - 1] + prefmatrix[preva][C - 1] - prefmatrix[preva][prevb - 1] + prefmatrix[a][b] - matrix[a][1];
    		fx = preva*2; fy = C; sx = a*2; sy = C;
    		ans4 = dist[2*preva][2*a] + prefmatrix[preva][C - 1] - prefmatrix[preva][prevb - 1] + prefmatrix[a][C - 1] - prefmatrix[a][b - 1];
			// cout << min4(ans1, ans2, ans3, ans4) << endl;    		
    		sum = sum + min4(ans1, ans2, ans3, ans4);
    	}
    	preva = a; prevb = b;
    }

    cout << sum << endl;

    //#if !ONLINE_JUDGE
    //    cout << fixed << setprecision(12) << "-------------------------------------------------\n";
    //    cout << double(clock())/CLOCKS_PER_SEC << " seconds" << endl;
    //#endif
    return 0;
}

