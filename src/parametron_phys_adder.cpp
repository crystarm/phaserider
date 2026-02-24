#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <iomanip>
#include <cmath>
#include <cstdint>
#include <limits>

#include "parametron_phys.h"

using namespace std;

typedef long long ll;
typedef unsigned long long ull;

static const int OK = 0;
static const int ERR_USAGE = 1;

struct Node {
    enum Kind {
        CONST,
        INPUT,
        MAJ,
        NOT
    } kind;
    vector<int> in;
    int fixed_bit;
};

static unordered_map<string, string> parse_args(int argc, char** argv) {
    unordered_map<string, string> m;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a.rfind("--", 0) == 0) {
            string key = a.substr(2);
            string val = "1";
            if (i + 1 < argc) {
                string nxt = argv[i + 1];
                if (nxt.rfind("--", 0) != 0) {
                    val = nxt;
                    i++;
                }
            }
            m[key] = val;
        }
    }
    return m;
}

static int parse_int(const string& v, int& out) {
    try {
        size_t p = 0;
        int x = stoi(v, &p, 10);
        if (p != v.size()) return -1;
        out = x;
        return 0;
    } catch (...) {
        return -1;
    }
}

static int parse_double(const string& v, double& out) {
    try {
        size_t p = 0;
        double x = stod(v, &p);
        if (p != v.size()) return -1;
        out = x;
        return 0;
    } catch (...) {
        return -1;
    }
}

static int parse_s64(const string& s, ll& out) {
    if (s.empty()) return -1;

    int i = 0;
    int sign = 1;
    if (s[i] == '+') i++;
    else if (s[i] == '-') {
        sign = -1;
        i++;
    }
    if (i >= (int)s.size()) return -1;

    int base = 10;
    if (s[i] == '0' && i + 1 < (int)s.size() && (s[i + 1] == 'x' || s[i + 1] == 'X')) {
        base = 16;
        i += 2;
    } else if (s[i] == '0' && i + 1 < (int)s.size() && (s[i + 1] == 'b' || s[i + 1] == 'B')) {
        base = 2;
        i += 2;
    }

    __int128 acc = 0;
    bool any = false;
    for (; i < (int)s.size(); i++) {
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

    acc *= (__int128)sign;

    __int128 mn = (__int128)numeric_limits<ll>::min();
    __int128 mx = (__int128)numeric_limits<ll>::max();
    if (acc < mn || acc > mx) return -1;

    out = (ll)acc;
    return 0;
}

static ull mask_for_width(int width) {
    if (width >= 64) return ~0ull;
    return (1ull << width) - 1ull;
}

static int add_node(vector<Node>& nodes, Node::Kind kind, const vector<int>& in = {}, int fixed = 0) {
    Node n;
    n.kind = kind;
    n.in = in;
    n.fixed_bit = fixed & 1;
    nodes.push_back(n);
    return (int)nodes.size() - 1;
}

static int make_const(vector<Node>& nodes, int v) {
    return add_node(nodes, Node::CONST, {}, v & 1);
}

static int make_input(vector<Node>& nodes, int v) {
    return add_node(nodes, Node::INPUT, {}, v & 1);
}

static int make_not(vector<Node>& nodes, int x) {
    return add_node(nodes, Node::NOT, {x}, 0);
}

static int make_maj(vector<Node>& nodes, int a, int b, int c) {
    return add_node(nodes, Node::MAJ, {a, b, c}, 0);
}

static int op_and(vector<Node>& nodes, int a, int b, int c0) {
    return make_maj(nodes, a, b, c0);
}

static int op_or(vector<Node>& nodes, int a, int b, int c1) {
    return make_maj(nodes, a, b, c1);
}

static int op_xor(vector<Node>& nodes, int a, int b, int c0, int c1) {
    int nb = make_not(nodes, b);
    int na = make_not(nodes, a);
    int t1 = op_and(nodes, a, nb, c0);
    int t2 = op_and(nodes, na, b, c0);
    return op_or(nodes, t1, t2, c1);
}

static void print_usage() {
    cout << "usage:\n";
    cout << "  ./parametron_phys_adder --a A --b B [--cin 0|1] [options]\n\n";
    cout << "options:\n";
    cout << "  --width N            (default 8)\n";
    cout << "  --k_majority X       (default 0.05)\n";
    cout << "  --k_invert X         (default 0.05)\n";
    cout << "  --normalize 0|1      (default 0)\n";
    cout << "  --hard 0|1           (default 1)\n";
    cout << "  --emit 0|1           (default 0)\n\n";
    cout << "physics:\n";
    cout << "  --f0 HZ              (default 50000)\n";
    cout << "  --dt SECONDS         (default 1e-6)\n";
    cout << "  --t SECONDS          (default 0.2)\n";
    cout << "  --zeta X             (default 0.010)\n";
    cout << "  --h X                (default 0.60)\n";
    cout << "  --alpha_ratio X      (default 1.0)\n";
    cout << "  --det_tau SECONDS    (default 0.010)\n";
    cout << "  --noise X            (default 0.0)\n";
    cout << "  --kick_amp X         (default 0.20)\n";
    cout << "  --kick_phase RAD     (default 0.0)\n";
    cout << "  --kick_dur SEC       (default 0.005)\n";
    cout << "  --out_every N        (default 10)\n";
    cout << "  --seed N             (default -1 time-based)\n";
    cout << "  --min_periods N      (default 3)\n";
    cout << "  --amp_eps X          (default 1e-4)\n";
    cout << "  --max_periods N      (default 200)\n";
    cout << "  --amp_min X          (default 1e-3)\n\n";
    cout << "examples:\n";
    cout << "  ./parametron_phys_adder --a 13 --b 250 --cin 1 --width 8\n";
}

static double effective_k(double base, const parametron_phys::CouplingConfig& cfg) {
    double scale = cfg.hard_logic ? 2.0 : 1.0;
    return base * scale;
}

int main(int argc, char** argv) {
    auto args = parse_args(argc, argv);

    if (!args.count("a") || !args.count("b") || args.count("help")) {
        print_usage();
        return ERR_USAGE;
    }

    int width = 8;
    int cin = 0;
    int emit = 0;

    parametron_phys::CouplingConfig cfg;
    int normalize = 0;
    int hard = 1;

    parametron_phys::PhysParams base;

    if (args.count("width") && parse_int(args["width"], width) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("cin") && (parse_int(args["cin"], cin) != 0 || (cin != 0 && cin != 1))) { print_usage(); return ERR_USAGE; }
    if (args.count("emit") && (parse_int(args["emit"], emit) != 0 || (emit != 0 && emit != 1))) { print_usage(); return ERR_USAGE; }

    if (args.count("k_majority") && parse_double(args["k_majority"], cfg.k_majority) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("k_invert") && parse_double(args["k_invert"], cfg.k_invert) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("normalize") && (parse_int(args["normalize"], normalize) != 0 || (normalize != 0 && normalize != 1))) { print_usage(); return ERR_USAGE; }
    if (args.count("hard") && (parse_int(args["hard"], hard) != 0 || (hard != 0 && hard != 1))) { print_usage(); return ERR_USAGE; }

    cfg.normalize_inputs = (normalize != 0);
    cfg.hard_logic = (hard != 0);

    ll a_in = 0, b_in = 0;
    if (parse_s64(args["a"], a_in) != 0 || parse_s64(args["b"], b_in) != 0) { print_usage(); return ERR_USAGE; }

    if (width < 1 || width > 64) { print_usage(); return ERR_USAGE; }

    if (args.count("f0") && parse_double(args["f0"], base.f0) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("dt") && parse_double(args["dt"], base.dt) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("t") && parse_double(args["t"], base.t_end) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("zeta") && parse_double(args["zeta"], base.zeta) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("h") && parse_double(args["h"], base.h) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("alpha_ratio") && parse_double(args["alpha_ratio"], base.alpha_ratio) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("det_tau") && parse_double(args["det_tau"], base.det_tau) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("noise") && parse_double(args["noise"], base.noise) != 0) { print_usage(); return ERR_USAGE; }

    if (args.count("kick_amp") && parse_double(args["kick_amp"], base.kick_amp) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("kick_phase") && parse_double(args["kick_phase"], base.kick_phase) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("kick_dur") && parse_double(args["kick_dur"], base.kick_dur) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("out_every") && parse_int(args["out_every"], base.out_every) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("seed") && parse_int(args["seed"], base.seed) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("min_periods") && parse_int(args["min_periods"], base.min_periods) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("amp_eps") && parse_double(args["amp_eps"], base.amp_eps) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("max_periods") && parse_int(args["max_periods"], base.max_periods) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("amp_min") && parse_double(args["amp_min"], base.amp_min) != 0) { print_usage(); return ERR_USAGE; }

    ull mask = mask_for_width(width);
    ull a = (ull)a_in & mask;
    ull b = (ull)b_in & mask;

    vector<Node> nodes;
    int c0 = make_const(nodes, 0);
    int c1 = make_const(nodes, 1);

    vector<int> a_in_nodes(width);
    vector<int> b_in_nodes(width);

    for (int i = 0; i < width; i++) {
        int abit = (int)((a >> i) & 1ull);
        int bbit = (int)((b >> i) & 1ull);
        a_in_nodes[i] = make_input(nodes, abit);
        b_in_nodes[i] = make_input(nodes, bbit);
    }

    int cin_node = make_input(nodes, cin & 1);

    vector<int> sum_nodes;
    sum_nodes.reserve(width);

    int carry = cin_node;
    for (int i = 0; i < width; i++) {
        int t = op_xor(nodes, a_in_nodes[i], b_in_nodes[i], c0, c1);
        int si = op_xor(nodes, t, carry, c0, c1);
        int co = make_maj(nodes, a_in_nodes[i], b_in_nodes[i], carry);
        sum_nodes.push_back(si);
        carry = co;
    }

    int node_count = (int)nodes.size();
    vector<vector<double>> k(node_count, vector<double>(node_count, 0.0));

    for (int out = 0; out < node_count; out++) {
        const Node& n = nodes[out];
        if (n.kind == Node::MAJ) {
            double base_k = effective_k(cfg.k_majority, cfg);
            double gain = parametron_phys::coupling_gain(n.in.size(), cfg);
            for (int in : n.in) {
                k[out][in] += base_k * gain;
            }
        } else if (n.kind == Node::NOT) {
            double base_k = effective_k(cfg.k_invert, cfg);
            double gain = parametron_phys::coupling_gain(n.in.size(), cfg);
            for (int in : n.in) {
                k[out][in] -= base_k * gain;
            }
        }
    }

    vector<parametron_phys::PhysParams> ps(node_count, base);

    for (int idx = 0; idx < node_count; idx++) {
        if (nodes[idx].kind == Node::CONST || nodes[idx].kind == Node::INPUT) {
            ps[idx].kick_bit = nodes[idx].fixed_bit;
        } else {
            ps[idx].kick_bit = 1;
        }
    }

    auto result = parametron_phys::run_coupled(ps, k, emit != 0);

    ull sum = 0;
    for (int i = 0; i < width; i++) {
        int node_id = sum_nodes[i];
        int bit = parametron_phys::bit_from_state(result.final_states[node_id]) & 1;
        sum |= (ull)bit << i;
    }

    int coutv = parametron_phys::bit_from_state(result.final_states[carry]) & 1;

    ull sign = (width >= 64) ? (1ull << 63) : (1ull << (width - 1));
    int C = coutv & 1;
    int V = ((~(a ^ b) & (a ^ sum) & sign) != 0) ? 1 : 0;
    int Z = (sum == 0) ? 1 : 0;
    int N = ((sum & sign) != 0) ? 1 : 0;

    cout << "phys_adder width=" << width
         << " k_maj=" << cfg.k_majority
         << " k_inv=" << cfg.k_invert
         << " normalize=" << (cfg.normalize_inputs ? 1 : 0)
         << " hard=" << (cfg.hard_logic ? 1 : 0)
         << "\n";

    cout << "a=" << (unsigned long long)a << " b=" << (unsigned long long)b << " cin=" << (cin & 1) << "\n";
    cout << "sum=" << (unsigned long long)(sum & mask) << " cout=" << coutv << "\n";
    cout << "C=" << C << " V=" << V << " Z=" << Z << " N=" << N << "\n";

    return OK;
}