#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <random>
#include <chrono>
#include <iomanip>
#include <cmath>

#include "net_tick.h"

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

struct sim_params
{
    int ticks;
    double p_flip;
    long long seed;
    int trace;
};

struct trace_ctx
{
    int width;
};

static void on_trace(void* ctx, int t, ull sum, int coutv)
{
    trace_ctx* tc = (trace_ctx*)ctx;
    (void)tc;
    cout << setw(3) << t << " " << bin8((int)(sum & 255ull)) << " " << (coutv & 1) << "\n";
}

static int run_once(int aval, int bval, int cin, int maj_delay, int not_delay, const sim_params& sp)
{
    net_tick_params p;
    p.width = 8;
    p.cin = cin & 1;
    p.ticks = sp.ticks;
    p.maj_delay = maj_delay;
    p.not_delay = not_delay;
    p.p_flip = sp.p_flip;

    ull seed = 0;
    if (sp.seed >= 0) seed = (ull)sp.seed;
    else seed = (ull)chrono::high_resolution_clock::now().time_since_epoch().count();
    p.seed = (long long)seed;

    int crit = net_tick_crit_ticks(8, maj_delay, not_delay);
    int ticks = p.ticks;
    if (ticks <= 0) ticks = crit;

    trace_ctx tc;
    tc.width = 8;
    net_tick_trace_fn cb = nullptr;

    if (sp.trace != 0)
    {
        cout << "seed=" << seed << " ticks=" << ticks
             << " maj_delay=" << maj_delay << " not_delay=" << not_delay
             << " p_flip=" << fixed << setprecision(6) << sp.p_flip
             << " crit=" << crit
             << "\n";
        cout << "t sum cout\n";
        cb = on_trace;
    }

    net_tick_result rr = net_tick_add((ull)(aval & 255), (ull)(bval & 255), p, cb, &tc);
    int sum = (int)(rr.sum & 255ull);
    int coutv = rr.carry & 1;

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
    long long seed = -1;
    int trace = 1;

    if (args.count("maj_delay") && parse_int(args["maj_delay"], maj_delay) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("not_delay") && parse_int(args["not_delay"], not_delay) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("ticks") && parse_int(args["ticks"], ticks) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("p_flip") && parse_double(args["p_flip"], p_flip) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("seed"))
    {
        int tmp = 0;
        if (parse_int(args["seed"], tmp) != 0) { print_usage(); return ERR_USAGE; }
        seed = (long long)tmp;
    }
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
