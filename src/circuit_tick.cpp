#include "circuit_tick.h"

#include <vector>
#include <random>
#include <chrono>
#include <algorithm>

using namespace std;

static ull mask_for_width(int width)
{
    if (width >= 64) return ~0ull;
    return (1ull << width) - 1ull;
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
    vector<int> carry_in;
    int coutv;
    vector<int> depth;
    int width;
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

static circuit build_adder(int width, ull aval, ull bval, int cin, int maj_delay, int not_delay)
{
    circuit c;
    c.width = width;
    c.ns.clear();

    c.c0 = make_const(c, 0);
    c.c1 = make_const(c, 1);

    c.a.resize(width);
    c.b.resize(width);
    for (int i = 0; i < width; i++)
    {
        c.a[i] = make_input(c, (int)((aval >> i) & 1ull));
        c.b[i] = make_input(c, (int)((bval >> i) & 1ull));
    }
    c.cin = make_input(c, cin & 1);

    c.sum.clear();
    c.sum.reserve(width);
    c.carry_in.clear();
    c.carry_in.reserve(width);

    int carry = c.cin;
    for (int i = 0; i < width; i++)
    {
        c.carry_in.push_back(carry);
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
                if (ud(rng) < p_flip) res ^= 1;
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

static ull read_sum(const circuit& c)
{
    ull s = 0;
    for (int i = 0; i < c.width; i++)
    {
        int b = out_bit(c.ns, c.sum[i]);
        s |= (ull)b << i;
    }
    return s & mask_for_width(c.width);
}

static ull read_carry_pack(const circuit& c)
{
    ull x = 0;
    for (int i = 0; i < c.width; i++)
    {
        int b = out_bit(c.ns, c.carry_in[i]);
        x |= (ull)b << i;
    }
    return x & mask_for_width(c.width);
}

static int critical_ticks(const circuit& c)
{
    int d = 0;
    for (int id : c.sum) d = max(d, c.depth[id]);
    d = max(d, c.depth[c.coutv]);
    return d + 2;
}

int circuit_tick_crit_ticks(int width, int maj_delay, int not_delay)
{
    if (width < 1) width = 1;
    if (width > 64) width = 64;
    if (maj_delay < 0) maj_delay = 0;
    if (not_delay < 0) not_delay = 0;

    circuit c = build_adder(width, 0ull, 0ull, 0, maj_delay, not_delay);
    return critical_ticks(c);
}

circuit_tick_result circuit_tick_add(ull a, ull b, const circuit_tick_params& p, circuit_tick_trace_fn trace, void* trace_ctx)
{
    circuit_tick_result r;
    r.sum = 0;
    r.carry = 0;
    r.used_seed = 0;
    r.used_ticks = 0;
    r.crit_ticks = 0;

    int width = p.width;
    if (width < 1) width = 1;
    if (width > 64) width = 64;

    ull mask = mask_for_width(width);
    a &= mask;
    b &= mask;

    int maj_delay = p.maj_delay;
    int not_delay = p.not_delay;
    if (maj_delay < 0) maj_delay = 0;
    if (not_delay < 0) not_delay = 0;

    circuit c = build_adder(width, a, b, p.cin & 1, maj_delay, not_delay);

    int crit = critical_ticks(c);
    int ticks = p.ticks;
    if (ticks <= 0) ticks = crit;

    ull seed = 0;
    if (p.seed >= 0) seed = (ull)p.seed;
    else seed = (ull)chrono::high_resolution_clock::now().time_since_epoch().count();

    mt19937_64 rng(seed);

    for (int t = 0; t <= ticks; t++)
    {
        ull sum = read_sum(c);
        ull carry_pack = read_carry_pack(c);
        int coutv = out_bit(c.ns, c.coutv);
        if (trace) trace(trace_ctx, t, sum, carry_pack, coutv);
        if (t == ticks) break;
        step_tick(c, rng, p.p_flip);
    }

    r.sum = read_sum(c);
    r.carry = out_bit(c.ns, c.coutv) & 1;
    r.used_seed = seed;
    r.used_ticks = ticks;
    r.crit_ticks = crit;
    return r;
}
