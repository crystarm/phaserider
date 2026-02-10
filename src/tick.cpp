#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <random>
#include <chrono>
#include <iomanip>
#include <cmath>

using namespace std;

typedef long long ll;

static const int OK = 0;
static const int ERR_USAGE = 1;
static const int ERR_FAIL = 2;

static unordered_map<string, string> parse_args(int argc, char** argv)
{
    unordered_map<string, string> m;
    for (int i = 1; i < argc; i++)
    {
        string a = argv[i];
        if (a.rfind("--", 0) == 0)
        {
            string key = a.substr(2);
            string val = "1";
            if (i + 1 < argc)
            {
                string nxt = argv[i + 1];
                if (nxt.rfind("--", 0) != 0)
                {
                    val = nxt;
                    i++;
                }
            }
            m[key] = val;
        }
    }
    return m;
}

static int parse_int(const string& v, int& out)
{
    try
    {
        size_t p = 0;
        int x = stoi(v, &p, 10);
        if (p != v.size()) return -1;
        out = x;
        return 0;
    }
    catch (...)
    {
        return -1;
    }
}

static int parse_double(const string& v, double& out)
{
    try
    {
        size_t p = 0;
        double x = stod(v, &p);
        if (p != v.size()) return -1;
        out = x;
        return 0;
    }
    catch (...)
    {
        return -1;
    }
}

static string bin8(int x)
{
    string s;
    for (int i = 7; i >= 0; i--)
    {
        s.push_back(((x >> i) & 1) ? '1' : '0');
    }
    return s;
}

static void print_usage()
{
    cout << "usage:\n";
    cout << "  ./parametron_tick --a A --b B [--cin 0|1] [options]\n";
    cout << "  ./parametron_tick --test [options]\n\n";

    cout << "options:\n";
    cout << "  --cin 0|1           (default 0)\n";
    cout << "  --ticks N           (default auto)\n";
    cout << "  --maj_delay N       (default 1)\n";
    cout << "  --not_delay N       (default 1)\n";
    cout << "  --p_flip X          (default 0.0)\n";
    cout << "  --seed N            (default 1 for test, time-based for run)\n";
    cout << "  --trace 0|1         (default 1 for run, 0 for test)\n";

    cout << "examples:\n";
    cout << "  ./parametron_tick --a 13 --b 250 --cin 1\n";
    cout << "  ./parametron_tick --a 13 --b 250 --cin 1 --maj_delay 2 --not_delay 1 --ticks 40\n";
    cout << "  ./parametron_tick --test --maj_delay 1 --not_delay 1\n";
}

struct node
{
    int kind;
    int delay;
    vector<int> in;
    vector<int> pipe;
    int fixed;
};

static int out_bit(const vector<node>& ns, int id)
{
    return ns[id].pipe[0] & 1;
}

struct circuit
{
    vector<node> ns;
    int c0;
    int c1;
    vector<int> a;
    vector<int> b;
    int cin;
    vector<int> sum;
    int coutv;
    vector<int> depth;
};

static int add_node(circuit& c, int kind, int delay, const vector<int>& in, int fixed)
{
    node n;
    n.kind = kind;
    n.in = in;
    n.fixed = fixed;

    if (kind == 2 || kind == 3)
    {
        if (delay < 1) delay = 1;
        n.delay = delay;
        n.pipe.assign(delay, 0);
    }
    else
    {
        n.delay = 0;
        n.pipe.assign(1, fixed & 1);
    }

    c.ns.push_back(n);
    return (int)c.ns.size() - 1;
}

static int make_const(circuit& c, int v)
{
    return add_node(c, 0, 0, {}, v);
}

static int make_input(circuit& c, int v)
{
    return add_node(c, 1, 0, {}, v);
}

static int make_not(circuit& c, int x, int not_delay)
{
    return add_node(c, 3, not_delay, {x}, 0);
}

static int make_maj(circuit& c, int a, int b, int d, int maj_delay)
{
    return add_node(c, 2, maj_delay, {a, b, d}, 0);
}

static int op_and(circuit& c, int a, int b, int maj_delay)
{
    return make_maj(c, a, b, c.c0, maj_delay);
}

static int op_or(circuit& c, int a, int b, int maj_delay)
{
    return make_maj(c, a, b, c.c1, maj_delay);
}

static int op_xor(circuit& c, int a, int b, int maj_delay, int not_delay)
{
    int nb = make_not(c, b, not_delay);
    int na = make_not(c, a, not_delay);
    int t1 = op_and(c, a, nb, maj_delay);
    int t2 = op_and(c, na, b, maj_delay);
    return op_or(c, t1, t2, maj_delay);
}

static circuit build_adder8(int aval, int bval, int cin, int maj_delay, int not_delay)
{
    circuit c;
    c.ns.clear();

    c.c0 = make_const(c, 0);
    c.c1 = make_const(c, 1);

    c.a.resize(8);
    c.b.resize(8);
    for (int i = 0; i < 8; i++)
    {
        c.a[i] = make_input(c, (aval >> i) & 1);
        c.b[i] = make_input(c, (bval >> i) & 1);
    }
    c.cin = make_input(c, cin & 1);

    c.sum.clear();
    c.sum.reserve(8);

    int carry = c.cin;
    for (int i = 0; i < 8; i++)
    {
        int t = op_xor(c, c.a[i], c.b[i], maj_delay, not_delay);
        int si = op_xor(c, t, carry, maj_delay, not_delay);
        int co = make_maj(c, c.a[i], c.b[i], carry, maj_delay);
        c.sum.push_back(si);
        carry = co;
    }
    c.coutv = carry;

    c.depth.assign((int)c.ns.size(), 0);
    for (int i = 0; i < (int)c.ns.size(); i++)
    {
        int d = 0;
        for (int j : c.ns[i].in)
        {
            d = max(d, c.depth[j]);
        }
        c.depth[i] = d + c.ns[i].delay;
    }

    return c;
}

struct sim_params
{
    int ticks;
    double p_flip;
    int seed;
    int trace;
};

static void step_tick(circuit& c, mt19937_64& rng, double p_flip)
{
    uniform_real_distribution<double> ud(0.0, 1.0);
    int n = (int)c.ns.size();
    vector<int> nxt(n, 0);

    for (int i = 0; i < n; i++)
    {
        auto& nd = c.ns[i];
        int res = 0;

        if (nd.kind == 0 || nd.kind == 1)
        {
            res = nd.fixed & 1;
        }
        else if (nd.kind == 3)
        {
            int a = out_bit(c.ns, nd.in[0]);
            res = a ^ 1;

            if (p_flip > 0.0)
            {
                double pe = p_flip;
                if (ud(rng) < pe) res ^= 1;
            }
        }
        else
        {
            int a = out_bit(c.ns, nd.in[0]);
            int b = out_bit(c.ns, nd.in[1]);
            int d = out_bit(c.ns, nd.in[2]);
            int s = a + b + d;
            res = (s >= 2) ? 1 : 0;

            if (p_flip > 0.0)
            {
                int margin = (s == 0 || s == 3) ? 3 : 1;
                double pe = p_flip / (double)margin;
                if (pe > 1.0) pe = 1.0;
                if (ud(rng) < pe) res ^= 1;
            }
        }

        nxt[i] = res;
    }

    for (int i = 0; i < n; i++)
    {
        auto& nd = c.ns[i];
        if (nd.kind == 0 || nd.kind == 1)
        {
            nd.pipe[0] = nd.fixed & 1;
            continue;
        }

        int d = nd.delay;
        if (d <= 1)
        {
            nd.pipe[0] = nxt[i] & 1;
        }
        else
        {
            for (int j = 0; j < d - 1; j++)
            {
                nd.pipe[j] = nd.pipe[j + 1];
            }
            nd.pipe[d - 1] = nxt[i] & 1;
        }
    }
}

static int read_sum8(const circuit& c)
{
    int s = 0;
    for (int i = 0; i < 8; i++)
    {
        int b = out_bit(c.ns, c.sum[i]);
        s |= (b << i);
    }
    return s & 255;
}

static int critical_ticks(const circuit& c)
{
    int d = 0;
    for (int id : c.sum) d = max(d, c.depth[id]);
    d = max(d, c.depth[c.coutv]);
    return d + 2;
}

static int run_once(int aval, int bval, int cin, int maj_delay, int not_delay, const sim_params& sp)
{
    circuit c = build_adder8(aval & 255, bval & 255, cin & 1, maj_delay, not_delay);

    int ticks = sp.ticks;
    if (ticks <= 0) ticks = critical_ticks(c);

    unsigned long long seed = 0;
    if (sp.seed >= 0) seed = (unsigned long long)sp.seed;
    else seed = (unsigned long long)chrono::high_resolution_clock::now().time_since_epoch().count();

    mt19937_64 rng(seed);

    if (sp.trace != 0)
    {
        cout << "seed=" << seed << " ticks=" << ticks
             << " maj_delay=" << maj_delay << " not_delay=" << not_delay
             << " p_flip=" << fixed << setprecision(6) << sp.p_flip
             << " crit=" << critical_ticks(c)
             << "\n";
        cout << "t sum cout\n";
    }

    for (int t = 0; t <= ticks; t++)
    {
        int sum = read_sum8(c);
        int coutv = out_bit(c.ns, c.coutv);

        if (sp.trace != 0)
        {
            cout << setw(3) << t << " " << bin8(sum) << " " << coutv << "\n";
        }

        if (t == ticks) break;
        step_tick(c, rng, sp.p_flip);
    }

    int sum = read_sum8(c);
    int coutv = out_bit(c.ns, c.coutv);

    int ref = (aval & 255) + (bval & 255) + (cin & 1);
    int ref_sum = ref & 255;
    int ref_cout = (ref >> 8) & 1;

    if (sp.trace != 0)
    {
        cout << "a=" << (aval & 255) << " (" << bin8(aval & 255) << ")\n";
        cout << "b=" << (bval & 255) << " (" << bin8(bval & 255) << ")\n";
        cout << "cin=" << (cin & 1) << "\n";
        cout << "sum=" << (sum & 255) << " (" << bin8(sum & 255) << ")\n";
        cout << "cout=" << (coutv & 1) << "\n";
        cout << "ref_sum=" << ref_sum << " ref_cout=" << ref_cout << "\n";
    }

    if ((sum & 255) != ref_sum || (coutv & 1) != ref_cout) return ERR_FAIL;
    return OK;
}

int main(int argc, char** argv)
{
    auto args = parse_args(argc, argv);

    int maj_delay = 1;
    int not_delay = 1;
    int ticks = 0;
    double p_flip = 0.0;
    int seed = -1;
    int trace = 1;

    if (args.count("maj_delay") && parse_int(args["maj_delay"], maj_delay) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("not_delay") && parse_int(args["not_delay"], not_delay) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("ticks") && parse_int(args["ticks"], ticks) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("p_flip") && parse_double(args["p_flip"], p_flip) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("seed") && parse_int(args["seed"], seed) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("trace") && parse_int(args["trace"], trace) != 0) { print_usage(); return ERR_USAGE; }

    if (maj_delay < 0 || not_delay < 0 || ticks < 0 || p_flip < 0.0)
    {
        print_usage();
        return ERR_USAGE;
    }

    sim_params sp;
    sp.ticks = ticks;
    sp.p_flip = p_flip;
    sp.seed = seed;
    sp.trace = trace;

    if (args.count("test"))
    {
        if (!args.count("trace")) sp.trace = 0;
        if (!args.count("seed")) sp.seed = 1;

        ll bad = 0;
        auto t0 = chrono::high_resolution_clock::now();

        for (int a = 0; a < 256; a++)
        {
            for (int b = 0; b < 256; b++)
            {
                for (int cin = 0; cin <= 1; cin++)
                {
                    sim_params sp2 = sp;
                    if (sp2.seed >= 0)
                    {
                        ll idx = ((ll)a << 9) ^ ((ll)b << 1) ^ (ll)cin;
                        sp2.seed = sp.seed + (int)(idx & 0x7fffffff);
                    }

                    int rc = run_once(a, b, cin, maj_delay, not_delay, sp2);
                    if (rc != OK)
                    {
                        bad++;
                        if (bad <= 10)
                        {
                            cout << "mismatch a=" << a << " b=" << b << " cin=" << cin << "\n";
                        }
                    }
                }
            }
        }

        auto t1 = chrono::high_resolution_clock::now();
        double ms = chrono::duration<double, milli>(t1 - t0).count();

        if (bad == 0)
        {
            cout << "ok: all cases passed (131072). time_ms=" << fixed << setprecision(2) << ms << "\n";
            return OK;
        }
        else
        {
            cout << "fail: bad=" << bad << " time_ms=" << fixed << setprecision(2) << ms << "\n";
            return ERR_FAIL;
        }
    }

    if (!args.count("a") || !args.count("b"))
    {
        print_usage();
        return ERR_USAGE;
    }

    int aval = 0, bval = 0, cin = 0;
    if (parse_int(args["a"], aval) != 0 || parse_int(args["b"], bval) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("cin"))
    {
        if (parse_int(args["cin"], cin) != 0 || (cin != 0 && cin != 1)) { print_usage(); return ERR_USAGE; }
    }

    if (!args.count("trace")) sp.trace = 1;

    int rc = run_once(aval, bval, cin, maj_delay, not_delay, sp);
    if (rc == OK) return OK;

    cout << "result mismatch (increase ticks, or set p_flip=0 for deterministic)\n";
    return rc;
}
