#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <random>
#include <chrono>
#include <iomanip>
#include <cstdint>
#include <limits>
#include <algorithm>

#include "circuit_tick.h"

using namespace std;

typedef long long ll;
typedef unsigned __int128 u128;
typedef __int128 s128;

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

static int parse_s64(const string& s, ll& out)
{
    if (s.empty()) return -1;

    int i = 0;
    int sign = 1;
    if (s[i] == '+') i++;
    else if (s[i] == '-')
    {
        sign = -1;
        i++;
    }
    if (i >= (int)s.size()) return -1;

    int base = 10;
    if (s[i] == '0' && i + 1 < (int)s.size() && (s[i + 1] == 'x' || s[i + 1] == 'X'))
    {
        base = 16;
        i += 2;
    }
    else if (s[i] == '0' && i + 1 < (int)s.size() && (s[i + 1] == 'b' || s[i + 1] == 'B'))
    {
        base = 2;
        i += 2;
    }

    s128 acc = 0;
    bool any = false;
    for (; i < (int)s.size(); i++)
    {
        char c = s[i];
        if (c == '_') continue;

        int d = -1;
        if ('0' <= c && c <= '9') d = c - '0';
        else if ('a' <= c && c <= 'f') d = 10 + (c - 'a');
        else if ('A' <= c && c <= 'F') d = 10 + (c - 'A');
        else return -1;

        if (d >= base) return -1;

        acc = acc * base + d;
        any = true;
    }

    if (!any) return -1;

    acc *= (s128)sign;

    s128 mn = (s128)numeric_limits<ll>::min();
    s128 mx = (s128)numeric_limits<ll>::max();
    if (acc < mn || acc > mx) return -1;

    out = (ll)acc;
    return 0;
}

static ull mask_for_width(int width)
{
    if (width >= 64) return ~0ull;
    return (1ull << width) - 1ull;
}

static string u128_to_dec(u128 x)
{
    if (x == 0) return "0";
    string s;
    while (x > 0)
    {
        int d = (int)(x % 10);
        s.push_back(char('0' + d));
        x /= 10;
    }
    reverse(s.begin(), s.end());
    return s;
}

static string bin_n(ull x, int width)
{
    string s;
    s.reserve(width);
    for (int i = width - 1; i >= 0; i--)
    {
        s.push_back(((x >> i) & 1ull) ? '1' : '0');
    }
    return s;
}

static string hex_n(ull x, int width)
{
    static const char* dig = "0123456789abcdef";
    int n = (width + 3) / 4;
    if (n < 1) n = 1;

    string s;
    s.reserve(n);
    for (int i = n - 1; i >= 0; i--)
    {
        int shift = i * 4;
        int v = (int)((x >> shift) & 0xfull);
        s.push_back(dig[v]);
    }
    return s;
}

static ll as_signed(ull x, int width)
{
    if (width >= 64) return (ll)(int64_t)x;
    ull sign = 1ull << (width - 1);
    if (x & sign)
    {
        ull base = 1ull << width;
        ull y = x - base;
        return (ll)(int64_t)y;
    }
    return (ll)(int64_t)x;
}

static int maj3(int a, int b, int c)
{
    int s = a + b + c;
    return (s >= 2) ? 1 : 0;
}

static int inv1(int a)
{
    return a ^ 1;
}

static int and2(int a, int b)
{
    return maj3(a, b, 0);
}

static int or2(int a, int b)
{
    return maj3(a, b, 1);
}

static int xor2(int a, int b)
{
    int nb = inv1(b);
    int na = inv1(a);
    int t1 = and2(a, nb);
    int t2 = and2(na, b);
    return or2(t1, t2);
}

struct add_out
{
    ull sum;
    int carry;
};

static add_out add_logic(ull a, ull b, int cin, int width)
{
    ull res = 0;
    int carry = cin & 1;
    for (int i = 0; i < width; i++)
    {
        int ai = (int)((a >> i) & 1ull);
        int bi = (int)((b >> i) & 1ull);

        int t = xor2(ai, bi);
        int si = xor2(t, carry);
        int co = maj3(ai, bi, carry);

        res |= (ull)si << i;
        carry = co;
    }

    add_out o;
    o.sum = res & mask_for_width(width);
    o.carry = carry & 1;
    return o;
}

struct flags_out
{
    int c;
    int v;
    int z;
    int n;
};

static flags_out calc_flags(ull a, ull b, ull sum, int carry, int width)
{
    ull mask = mask_for_width(width);
    a &= mask;
    b &= mask;
    sum &= mask;

    ull sign = (width >= 64) ? (1ull << 63) : (1ull << (width - 1));

    int z = (sum == 0) ? 1 : 0;
    int n = ((sum & sign) != 0) ? 1 : 0;
    int v = ((~(a ^ b) & (a ^ sum) & sign) != 0) ? 1 : 0;

    flags_out f;
    f.c = carry & 1;
    f.v = v;
    f.z = z;
    f.n = n;
    return f;
}

struct tick_params
{
    int ticks;
    int maj_delay;
    int not_delay;
    double p_flip;
    int seed;
    int trace;
};

struct tick_trace_ctx
{
    int width;
};

static void on_tick_trace(void* ctx, int t, ull sum, ull carry_pack, int coutv)
{
    tick_trace_ctx* tc = (tick_trace_ctx*)ctx;
    (void)carry_pack;
    cout << setw(3) << t << " " << bin_n(sum, tc->width) << " " << (coutv & 1) << "\n";
}

static add_out add_tick(ull a, ull b, int cin, int width, const tick_params& tp)
{
    circuit_tick_params p;
    p.width = width;
    p.cin = cin & 1;
    p.ticks = tp.ticks;
    p.maj_delay = tp.maj_delay;
    p.not_delay = tp.not_delay;
    p.p_flip = tp.p_flip;

    ull seed = 0;
    if (tp.seed >= 0) seed = (ull)tp.seed;
    else seed = (ull)chrono::high_resolution_clock::now().time_since_epoch().count();
    p.seed = (long long)seed;

    int crit = circuit_tick_crit_ticks(width, tp.maj_delay, tp.not_delay);
    int ticks = p.ticks;
    if (ticks <= 0) ticks = crit;

    tick_trace_ctx tc;
    tc.width = width;

    circuit_tick_trace_fn cb = nullptr;
    if (tp.trace != 0)
    {
        cout << "seed=" << seed << " ticks=" << ticks
             << " maj_delay=" << tp.maj_delay << " not_delay=" << tp.not_delay
             << " p_flip=" << fixed << setprecision(6) << tp.p_flip
             << " crit=" << crit << "\n";
        cout << "t sum cout\n";
        cb = on_tick_trace;
    }

    circuit_tick_result rr = circuit_tick_add(a, b, p, cb, &tc);

    add_out o;
    o.sum = rr.sum & mask_for_width(width);
    o.carry = rr.carry & 1;
    return o;
}

static void print_usage()
{
    cout << "usage:\n";
    cout << "  ./phazeride_adc --a VAL --b VAL [--cin 0|1] [options]\n\n";

    cout << "options:\n";
    cout << "  --width N             (default 8, 1..64)\n";
    cout << "  --model logic|tick    (default tick)\n";
    cout << "  --signed 0|1          (default 0; affects decimal printing)\n";
    cout << "  --fmt dec|hex|bin|all (default all)\n";
    cout << "  --strict 0|1          (default 0; if 1, mismatch vs ref -> exit !=0)\n\n";

    cout << "tick options (model=tick):\n";
    cout << "  --maj_delay N         (default 1)\n";
    cout << "  --not_delay N         (default 1)\n";
    cout << "  --ticks N             (default auto)\n";
    cout << "  --p_flip X            (default 0.0)\n";
    cout << "  --seed N              (default time-based; irrelevant if p_flip=0)\n";
    cout << "  --trace 0|1           (default 0)\n\n";

    cout << "examples:\n";
    cout << "  ./phazeride_adc --a 250 --b 13 --cin 1 --width 8 --model tick --trace 1\n";
    cout << "  ./phazeride_adc --a 0xfa --b 0x0d --cin 1 --width 8 --model logic\n";
    cout << "  ./phazeride_adc --a -1 --b 1 --width 8 --signed 1\n";
}

static void print_val(const string& name, ull v, int width, int signed_mode, const string& fmt)
{
    ull mask = mask_for_width(width);
    v &= mask;

    if (fmt == "dec")
    {
        if (signed_mode != 0) cout << name << "=" << as_signed(v, width) << "\n";
        else cout << name << "=" << u128_to_dec((u128)v) << "\n";
        return;
    }

    if (fmt == "hex")
    {
        cout << name << "=0x" << hex_n(v, width) << "\n";
        return;
    }

    if (fmt == "bin")
    {
        cout << name << "=0b" << bin_n(v, width) << "\n";
        return;
    }

    if (signed_mode != 0)
    {
        cout << name << "=" << as_signed(v, width) << " (0x" << hex_n(v, width) << ") (" << bin_n(v, width) << ")\n";
    }
    else
    {
        cout << name << "=" << u128_to_dec((u128)v) << " (0x" << hex_n(v, width) << ") (" << bin_n(v, width) << ")\n";
    }
}

int main(int argc, char** argv)
{
    auto args = parse_args(argc, argv);

    if (args.count("help") || args.count("h"))
    {
        print_usage();
        return OK;
    }

    if (!args.count("a") || !args.count("b"))
    {
        print_usage();
        return ERR_USAGE;
    }

    int width = 8;
    int cin = 0;
    string model = "tick";
    string fmt = "all";
    int signed_mode = 0;
    int strict = 0;

    if (args.count("width") && parse_int(args["width"], width) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("cin") && (parse_int(args["cin"], cin) != 0 || (cin != 0 && cin != 1))) { print_usage(); return ERR_USAGE; }
    if (args.count("model")) model = args["model"];
    if (args.count("fmt")) fmt = args["fmt"];
    if (args.count("signed") && (parse_int(args["signed"], signed_mode) != 0 || (signed_mode != 0 && signed_mode != 1))) { print_usage(); return ERR_USAGE; }
    if (args.count("strict") && (parse_int(args["strict"], strict) != 0 || (strict != 0 && strict != 1))) { print_usage(); return ERR_USAGE; }

    if (width < 1 || width > 64)
    {
        print_usage();
        return ERR_USAGE;
    }

    if (!(fmt == "dec" || fmt == "hex" || fmt == "bin" || fmt == "all"))
    {
        print_usage();
        return ERR_USAGE;
    }

    ll a_in = 0, b_in = 0;
    if (parse_s64(args["a"], a_in) != 0 || parse_s64(args["b"], b_in) != 0)
    {
        print_usage();
        return ERR_USAGE;
    }

    ull mask = mask_for_width(width);
    ull a = ((ull)a_in) & mask;
    ull b = ((ull)b_in) & mask;

    tick_params tp;
    tp.maj_delay = 1;
    tp.not_delay = 1;
    tp.ticks = 0;
    tp.p_flip = 0.0;
    tp.seed = -1;
    tp.trace = 0;

    if (args.count("maj_delay") && parse_int(args["maj_delay"], tp.maj_delay) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("not_delay") && parse_int(args["not_delay"], tp.not_delay) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("ticks") && parse_int(args["ticks"], tp.ticks) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("p_flip") && parse_double(args["p_flip"], tp.p_flip) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("seed") && parse_int(args["seed"], tp.seed) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("trace") && parse_int(args["trace"], tp.trace) != 0) { print_usage(); return ERR_USAGE; }

    if (tp.maj_delay < 0 || tp.not_delay < 0 || tp.ticks < 0 || tp.p_flip < 0.0)
    {
        print_usage();
        return ERR_USAGE;
    }

    add_out r;
    if (model == "logic")
    {
        r = add_logic(a, b, cin, width);
    }
    else if (model == "tick")
    {
        r = add_tick(a, b, cin, width, tp);
    }
    else
    {
        print_usage();
        return ERR_USAGE;
    }

    ull sum = r.sum & mask;
    int carry = r.carry & 1;

    u128 full = (u128)a + (u128)b + (u128)(cin & 1);
    ull ref_sum = (ull)(full & (u128)mask);
    int ref_carry = (int)((full >> width) & 1);

    flags_out fl = calc_flags(a, b, sum, carry, width);

    cout << "adc width=" << width << " model=" << model;
    if (model == "tick")
    {
        cout << " (maj_delay=" << tp.maj_delay << " not_delay=" << tp.not_delay;
        if (tp.ticks > 0) cout << " ticks=" << tp.ticks;
        else cout << " ticks=auto";
        cout << " p_flip=" << fixed << setprecision(6) << tp.p_flip << ")";
    }
    cout << "\n";

    print_val("a", a, width, signed_mode, fmt);
    print_val("b", b, width, signed_mode, fmt);
    cout << "cin=" << (cin & 1) << "\n";

    print_val("sum", sum, width, signed_mode, fmt);
    cout << "C=" << fl.c << " V=" << fl.v << " Z=" << fl.z << " N=" << fl.n << "\n";
    cout << "ref_full=" << u128_to_dec(full) << " ref_sum=0x" << hex_n(ref_sum, width) << " ref_cout=" << ref_carry << "\n";

    if (strict != 0)
    {
        if (sum != ref_sum || carry != ref_carry)
        {
            cerr << "mismatch vs reference (strict=1)\n";
            return ERR_FAIL;
        }
    }

    return OK;
}
