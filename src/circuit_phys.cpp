#include <cstdint>
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>

#include "parametron_phys.h"

using std::size_t;

typedef unsigned long long ull;

namespace parametron_phys_circuit {

// Physical adder configuration.
// This builds a full adder ripple chain using MAJ/NOT/XOR decomposition.
struct PhysAdderParams {
    int width = 8;

    // Base physical parameters (shared for all nodes).
    parametron_phys::PhysParams phys;

    // Coupling strengths.
    parametron_phys::CouplingConfig coupling;

    // Input drive for source nodes (A/B/CIN).
    double input_kick_amp = 0.20;
    double input_kick_dur = 0.005;

    // Constant nodes drive (for OR/AND construction).
    double const_kick_amp = 0.20;
    double const_kick_dur = 0.2; // default to full simulation horizon

    // If true, non-input logic nodes get a small kick to break symmetry.
    bool logic_kick = false;
    double logic_kick_amp = 0.02;
    double logic_kick_dur = 0.002;
};

// Result of the physical adder.
struct PhysAdderResult {
    ull sum = 0;
    int carry = 0;

    // Raw physics output for inspection.
    parametron_phys::CoupledResult phys;
};

enum NodeKind {
    NODE_CONST = 0,
    NODE_INPUT = 1,
    NODE_NOT = 2,
    NODE_MAJ = 3,
};

struct Node {
    NodeKind kind;
    std::vector<int> in;
    int init_bit = 0; // used for CONST/INPUT
};

struct AdderGraph {
    std::vector<Node> nodes;
    int c0 = -1;
    int c1 = -1;
    std::vector<int> a;
    std::vector<int> b;
    int cin = -1;
    std::vector<int> sum;
    int coutv = -1;
};

static int add_node(AdderGraph& g, NodeKind kind, const std::vector<int>& in, int init_bit)
{
    Node n;
    n.kind = kind;
    n.in = in;
    n.init_bit = init_bit;
    g.nodes.push_back(n);
    return (int)g.nodes.size() - 1;
}

static int make_const(AdderGraph& g, int v)
{
    return add_node(g, NODE_CONST, {}, v & 1);
}

static int make_input(AdderGraph& g, int v)
{
    return add_node(g, NODE_INPUT, {}, v & 1);
}

static int make_not(AdderGraph& g, int x)
{
    return add_node(g, NODE_NOT, {x}, 0);
}

static int make_maj(AdderGraph& g, int a, int b, int c)
{
    return add_node(g, NODE_MAJ, {a, b, c}, 0);
}

static int op_and(AdderGraph& g, int a, int b)
{
    return make_maj(g, a, b, g.c0);
}

static int op_or(AdderGraph& g, int a, int b)
{
    return make_maj(g, a, b, g.c1);
}

static int op_xor(AdderGraph& g, int a, int b)
{
    int nb = make_not(g, b);
    int na = make_not(g, a);
    int t1 = op_and(g, a, nb);
    int t2 = op_and(g, na, b);
    return op_or(g, t1, t2);
}

static AdderGraph build_adder_graph(int width, ull aval, ull bval, int cinv)
{
    AdderGraph g;
    g.nodes.clear();

    g.c0 = make_const(g, 0);
    g.c1 = make_const(g, 1);

    g.a.resize(width);
    g.b.resize(width);

    for (int i = 0; i < width; i++)
    {
        int abit = (int)((aval >> i) & 1ull);
        int bbit = (int)((bval >> i) & 1ull);
        g.a[i] = make_input(g, abit);
        g.b[i] = make_input(g, bbit);
    }

    g.cin = make_input(g, cinv & 1);

    g.sum.clear();
    g.sum.reserve(width);

    int carry = g.cin;
    for (int i = 0; i < width; i++)
    {
        int t = op_xor(g, g.a[i], g.b[i]);
        int si = op_xor(g, t, carry);
        int co = make_maj(g, g.a[i], g.b[i], carry);
        g.sum.push_back(si);
        carry = co;
    }
    g.coutv = carry;

    return g;
}

static void build_coupling_matrix(
    const AdderGraph& g,
    const parametron_phys::CouplingConfig& cfg,
    std::vector<std::vector<double>>& k)
{
    const size_t n = g.nodes.size();
    k.assign(n, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; i++)
    {
        const Node& nd = g.nodes[i];
        if (nd.kind == NODE_MAJ)
        {
            const double gain = cfg.k_majority * parametron_phys::coupling_gain(nd.in.size(), cfg);
            for (int src : nd.in)
            {
                k[i][(size_t)src] += gain;
            }
        }
        else if (nd.kind == NODE_NOT)
        {
            const double gain = cfg.k_invert * parametron_phys::coupling_gain(nd.in.size(), cfg);
            if (!nd.in.empty())
            {
                k[i][(size_t)nd.in[0]] += -gain;
            }
        }
    }
}

static std::vector<parametron_phys::PhysParams> build_node_params(
    const AdderGraph& g,
    const PhysAdderParams& p)
{
    std::vector<parametron_phys::PhysParams> ps;
    ps.reserve(g.nodes.size());

    for (const Node& nd : g.nodes)
    {
        parametron_phys::PhysParams cp = p.phys;

        if (nd.kind == NODE_INPUT)
        {
            cp.kick_bit = nd.init_bit;
            cp.kick_amp = p.input_kick_amp;
            cp.kick_dur = p.input_kick_dur;
        }
        else if (nd.kind == NODE_CONST)
        {
            cp.kick_bit = nd.init_bit;
            cp.kick_amp = p.const_kick_amp;
            cp.kick_dur = p.const_kick_dur;
        }
        else
        {
            cp.kick_amp = 0.0;
            cp.kick_dur = 0.0;
            if (p.logic_kick)
            {
                cp.kick_bit = 1;
                cp.kick_amp = p.logic_kick_amp;
                cp.kick_dur = p.logic_kick_dur;
            }
        }

        ps.push_back(cp);
    }

    return ps;
}

PhysAdderResult phys_adder_add(ull a, ull b, int cinv, const PhysAdderParams& p)
{
    if (p.width < 1 || p.width > 64)
    {
        throw std::runtime_error("phys_adder_add: width out of range");
    }

    AdderGraph g = build_adder_graph(p.width, a, b, cinv);

    std::vector<std::vector<double>> k;
    build_coupling_matrix(g, p.coupling, k);

    std::vector<parametron_phys::PhysParams> ps = build_node_params(g, p);

    parametron_phys::CoupledResult cr = parametron_phys::run_coupled(ps, k, false);

    PhysAdderResult out;
    out.phys = cr;

    if (cr.final_states.size() != g.nodes.size())
    {
        return out;
    }

    ull sum = 0;
    for (int i = 0; i < p.width; i++)
    {
        int node_id = g.sum[i];
        int bit = parametron_phys::bit_from_state(cr.final_states[(size_t)node_id]) & 1;
        sum |= (ull)bit << i;
    }

    int cout_bit = parametron_phys::bit_from_state(cr.final_states[(size_t)g.coutv]) & 1;

    out.sum = sum;
    out.carry = cout_bit;
    return out;
}

} // namespace parametron_phys_circuit