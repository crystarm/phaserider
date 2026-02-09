#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include <unordered_map>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace std;

typedef long long ll;

static const int OK = 0;
static const int ERR_USAGE = 1;

struct edge_in
{
    int src;
    double w;
};

struct net
{
    int n = 0;
    vector<vector<edge_in>> in;
    vector<int> s;
    vector<double> bias;
};

static int sign_bit(double x)
{
    return (x >= 0.0) ? +1 : -1;
}

static string bits01(const vector<int>& s)
{
    string out;
    out.reserve(s.size());
    for (int v : s)
    {
        out.push_back(v > 0 ? '1' : '0');
    }
    return out;
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

static void add_edge(net& g, int src, int dst, double w)
{
    g.in[dst].push_back(edge_in{src, w});
}

static int step_net(net& g,
                    const vector<double>& ext,
                    const vector<int>& clamp,
                    mt19937_64& rng,
                    double sigma,
                    double mem)
{
    normal_distribution<double> nd(0.0, sigma);
    vector<int> ns(g.n, -1);

    for (int i = 0; i < g.n; i++)
    {
        if (clamp[i] != 0)
        {
            ns[i] = clamp[i];
            continue;
        }

        double drive = g.bias[i] + ext[i];
        for (const auto& e : g.in[i])
        {
            drive += e.w * (double)g.s[e.src];
        }

        double x = mem * (double)g.s[i] + drive + nd(rng);
        ns[i] = sign_bit(x);
    }

    g.s.swap(ns);
    return 0;
}

static net make_demo_majority(vector<int>& clamp_out)
{
    net g;
    g.n = 4;
    g.in.assign(g.n, {});
    g.s.assign(g.n, -1);
    g.bias.assign(g.n, 0.0);

    add_edge(g, 0, 3, +1.0);
    add_edge(g, 1, 3, +1.0);
    add_edge(g, 2, 3, +1.0);

    clamp_out.assign(g.n, 0);
    return g;
}

static net make_demo_inverter(vector<int>& clamp_out)
{
    net g;
    g.n = 2;
    g.in.assign(g.n, {});
    g.s.assign(g.n, -1);
    g.bias.assign(g.n, 0.0);

    add_edge(g, 0, 1, -1.0);

    clamp_out.assign(g.n, 0);
    return g;
}

static void print_usage()
{
    cout << "usage:\n";
    cout << "  ./parametron_lab --demo majority|inverter [options]\n\n";
    cout << "options:\n";
    cout << "  --steps N         (default 30)\n";
    cout << "  --sigma X         (default 0.10) noise stdev\n";
    cout << "  --mem X           (default 0.80) inertia term\n";
    cout << "  --seed N          (default: time-based)\n";
    cout << "  --in a,b,c,...    input bits for demo (0/1), clamped\n";
    cout << "examples:\n";
    cout << "  ./parametron_lab --demo majority --in 1,0,1 --steps 25\n";
    cout << "  ./parametron_lab --demo inverter --in 1 --steps 20\n";
}

static int parse_in_bits(const string& s, vector<int>& bits)
{
    bits.clear();
    string cur;
    for (char c : s)
    {
        if (c == ',')
        {
            if (cur.empty()) return -1;
            bits.push_back((cur == "1") ? +1 : (cur == "0" ? -1 : 0));
            if (bits.back() == 0) return -1;
            cur.clear();
        }
        else
        {
            cur.push_back(c);
        }
    }
    if (!cur.empty())
    {
        bits.push_back((cur == "1") ? +1 : (cur == "0" ? -1 : 0));
        if (bits.back() == 0) return -1;
    }
    return 0;
}

int main(int argc, char** argv)
{
    auto args = parse_args(argc, argv);

    string demo = "majority";
    if (args.count("demo")) demo = args["demo"];

    int steps = 30;
    if (args.count("steps"))
    {
        if (parse_int(args["steps"], steps) != 0 || steps <= 0)
        {
            print_usage();
            return ERR_USAGE;
        }
    }

    double sigma = 0.10;
    if (args.count("sigma"))
    {
        if (parse_double(args["sigma"], sigma) != 0 || sigma < 0.0)
        {
            print_usage();
            return ERR_USAGE;
        }
    }

    double mem = 0.80;
    if (args.count("mem"))
    {
        if (parse_double(args["mem"], mem) != 0 || mem < 0.0)
        {
            print_usage();
            return ERR_USAGE;
        }
    }

    unsigned long long seed = 0;
    bool has_seed = false;
    if (args.count("seed"))
    {
        int tmp = 0;
        if (parse_int(args["seed"], tmp) != 0 || tmp < 0)
        {
            print_usage();
            return ERR_USAGE;
        }
        seed = (unsigned long long)tmp;
        has_seed = true;
    }
    if (!has_seed)
    {
        seed = (unsigned long long)chrono::high_resolution_clock::now().time_since_epoch().count();
    }

    vector<int> clamp;
    net g;
    if (demo == "majority")
    {
        g = make_demo_majority(clamp);
    }
    else if (demo == "inverter")
    {
        g = make_demo_inverter(clamp);
    }
    else
    {
        print_usage();
        return ERR_USAGE;
    }

    vector<int> in_bits;
    if (args.count("in"))
    {
        if (parse_in_bits(args["in"], in_bits) != 0)
        {
            print_usage();
            return ERR_USAGE;
        }
    }

    if (demo == "majority")
    {
        if (in_bits.empty()) in_bits = {+1, -1, +1};
        if ((int)in_bits.size() != 3)
        {
            print_usage();
            return ERR_USAGE;
        }
        clamp[0] = in_bits[0];
        clamp[1] = in_bits[1];
        clamp[2] = in_bits[2];
    }
    else if (demo == "inverter")
    {
        if (in_bits.empty()) in_bits = {+1};
        if ((int)in_bits.size() != 1)
        {
            print_usage();
            return ERR_USAGE;
        }
        clamp[0] = in_bits[0];
    }

    mt19937_64 rng(seed);
    vector<double> ext(g.n, 0.0);

    cout << "demo=" << demo
         << " steps=" << steps
         << " sigma=" << fixed << setprecision(3) << sigma
         << " mem=" << fixed << setprecision(3) << mem
         << " seed=" << seed
         << "\n";

    cout << "t  state\n";
    for (int t = 0; t <= steps; t++)
    {
        cout << setw(2) << t << " " << bits01(g.s) << "\n";
        if (t == steps) break;
        step_net(g, ext, clamp, rng, sigma, mem);
    }

    if (demo == "majority")
    {
        cout << "out(node3)=" << (g.s[3] > 0 ? 1 : 0) << "\n";
    }
    else
    {
        cout << "out(node1)=" << (g.s[1] > 0 ? 1 : 0) << "\n";
    }

    return OK;
}
