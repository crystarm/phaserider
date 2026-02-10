#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <chrono>
#include <iomanip>

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

static void add8(int a, int b, int cin, int& sum, int& coutv)
{
    int carry = cin;
    int res = 0;

    for (int i = 0; i < 8; i++)
    {
        int ai = (a >> i) & 1;
        int bi = (b >> i) & 1;

        int t = xor2(ai, bi);
        int si = xor2(t, carry);
        int co = maj3(ai, bi, carry);

        res |= (si << i);
        carry = co;
    }

    sum = res & 255;
    coutv = carry & 1;
}

static void print_usage()
{
    cout << "usage:\n";
    cout << "  ./parametron_logic --a A --b B [--cin 0|1]\n";
    cout << "  ./parametron_logic --test\n\n";
    cout << "examples:\n";
    cout << "  ./parametron_logic --a 13 --b 250 --cin 1\n";
    cout << "  ./parametron_logic --test\n";
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

int main(int argc, char** argv)
{
    auto args = parse_args(argc, argv);

    if (args.count("test"))
    {
        ll bad = 0;
        auto t0 = chrono::high_resolution_clock::now();

        for (int a = 0; a < 256; a++)
        {
            for (int b = 0; b < 256; b++)
            {
                for (int cin = 0; cin <= 1; cin++)
                {
                    int sum = 0, coutv = 0;
                    add8(a, b, cin, sum, coutv);

                    int ref = a + b + cin;
                    int ref_sum = ref & 255;
                    int ref_cout = (ref >> 8) & 1;

                    if (sum != ref_sum || coutv != ref_cout)
                    {
                        bad++;
                        if (bad <= 10)
                        {
                            cout << "mismatch a=" << a << " b=" << b << " cin=" << cin
                                 << " got(sum=" << sum << " cout=" << coutv << ")"
                                 << " ref(sum=" << ref_sum << " cout=" << ref_cout << ")\n";
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

    int a = -1, b = -1, cin = 0;
    if (!args.count("a") || !args.count("b"))
    {
        print_usage();
        return ERR_USAGE;
    }
    if (parse_int(args["a"], a) != 0 || parse_int(args["b"], b) != 0)
    {
        print_usage();
        return ERR_USAGE;
    }
    if (args.count("cin"))
    {
        if (parse_int(args["cin"], cin) != 0 || (cin != 0 && cin != 1))
        {
            print_usage();
            return ERR_USAGE;
        }
    }

    a &= 255;
    b &= 255;

    int sum = 0, coutv = 0;
    add8(a, b, cin, sum, coutv);

    int full = a + b + cin;

    cout << "a=" << a << " (" << bin8(a) << ")\n";
    cout << "b=" << b << " (" << bin8(b) << ")\n";
    cout << "cin=" << cin << "\n";
    cout << "sum=" << sum << " (" << bin8(sum) << ")\n";
    cout << "cout=" << coutv << "\n";
    cout << "ref=" << full << "\n";

    return OK;
}
