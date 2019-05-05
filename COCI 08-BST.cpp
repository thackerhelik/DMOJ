/*
While inserting a value in bst will go to the left child of the next bigger element in tree (if present)
or it will go to the right child of the just smaller element in tree(if present)
we can see the max height among the both (even if it has a fixed position) because if right child of just smaller is full
then that element will be the next greater element and its left child will be the place
Answer to be taken in long long since roughly upper bounded by n^2
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

const int nax = 300005;
int height[nax];
set<int>s;
ll cumulative_height;

int main(){
    #if ONLINE_JUDGE
        ios_base::sync_with_stdio(0);cin.tie(NULL);cout.tie(NULL);
	#endif

    int n, x;
    cin >> n;
    cin >> x;
    s.insert(x);
    cout << 0 << endl;
    for(int i = 1; i < n; ++i){
    	cin >> x;
    	auto it = s.lower_bound(x);
    	if(it == s.end()){
    		//right child of biggest element
    		int val = *s.rbegin();
    		height[x] = height[val] + 1;
    	}	
    	else if(it != s.begin()){
    		auto it2 = it;
    		cout << x << " " << *(it2) << " " << *(--it2) << endl;
    		height[x] = height[*it] + 1;
    		height[x] = max(height[x], height[*(--it)] + 1);
    	}
    	else
    		height[x] = height[*it] + 1;
    	s.insert(x);
    	cumulative_height += height[x];
    	cout << cumulative_height << endl;
    }        

    //#if !ONLINE_JUDGE
    //    cout << fixed << setprecision(12) << "-------------------------------------------------\n";
    //    cout << double(clock())/CLOCKS_PER_SEC << " seconds" << endl;
    //#endif
    return 0;
}

