//  Rigorous enclosures of minimal detectable differences for ANOVA models
//
//  Copyright (C) 2015   Ali Baharev
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Library General Public
//  License as published by the Free Software Foundation; either
//  version 2 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Library General Public License for more details.
//
//  You should have received a copy of the GNU Library General Public
//  License along with this library; if not, write to the Free
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//-----------------------------------------------------------------------------
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdexcept>
#include "ddf_ari.hpp"  // automatic differentiation in C-XSC
#include "nlfzero.hpp"  // univariate interval Newton method in C-XSC

namespace {

using cxsc::interval;

// I_x(a, b)
interval beta(const interval x,
              const interval a,
              const int      b)
{
    const interval one_minus_x(1-x);
    interval sum(0);
    interval prod(1);
    for (int n=1; n<b; ++n)  {
        prod *= (a+(n-1))/n;
        sum  += prod*power(one_minus_x, n);
    }
    return pow(x, a)*(1+sum);
}

//  I_x(a, b; lambda), currently unused function
interval nc_beta(const interval x,
                 const interval a,
                 const int      b,
                 const interval lambda)
{
    const interval lam_one_minus_x((lambda/2)*(1-x));
    interval sum = beta(x, a, b);
    interval fakt(1);
    for (int i=1; i<b; ++i) {
        fakt *= lam_one_minus_x/i;
        sum  += fakt*beta(x, a+i, b-i);
    }
    return exp(-lam_one_minus_x)*sum;
}

//-----------------------------------------------------------------------------
// Setting up the two univariate equations to be solved. The C-XSC's interval
// Newton solver expects a function pointer. Therefore, we first wrap up the
// state in two classes, corresponding to the two equations to be solved.

// For computing the p quantile
class quantile {
    const interval a;
    const int      b;
    const interval p;

public:

    quantile(interval a, int b, interval p) : a(a), b(b), p(p) { }

    DerivType residual(const DerivType& x) {
        const DerivType one_minus_x(1-x);
        DerivType sum  = DerivConst(0);
        DerivType prod = DerivConst(1);
        for (int n=1; n<b; ++n)  {
            prod = prod*((a+(n-1))/n);
            sum  = sum+prod*power(one_minus_x, n);
        }
        return power(x, a)*(1+sum)-p;
    }
};

// For computing the noncentrality parameter
class noncent_param {
    const interval x;
    const interval a;
    const int      b;
    const interval prob2;

public:

    noncent_param(interval x, interval a, int b, interval prob2)
        : x(x), a(a), b(b), prob2(prob2) { }

    DerivType residual(const DerivType& lambda) {
        const DerivType lam_one_minus_x((lambda/2)*(1-x));
        DerivType sum  = DerivConst(beta(x, a, b));
        DerivType fakt = DerivConst(1);
        for (int i=1; i<b; ++i) {
            fakt = fakt*(lam_one_minus_x/i);
            sum  = sum + fakt*beta(x, a+i, b-i);
        }
        return exp(-lam_one_minus_x)*sum-prob2;
    }
};

// We have to pass a function pointer to C-XSC, not a pointer to a member
// function or a functor. We have to circumvent it with global variables.

quantile* q_ptr;
noncent_param* np_ptr;

// These two functions will be passed to the solver as a function pointer.

DerivType beta_quantile(const DerivType& x) {
    return q_ptr->residual(x);
}

DerivType ncbeta_noncent(const DerivType& lambda) {
    return np_ptr->residual(lambda);
}

//-----------------------------------------------------------------------------
// Invoking the solver `AllZeros` in C-XSC, and checking the existence and
// uniqueness of the returned solution.

struct newton_result {
    enum outcome { HAS_SOLUTION, NO_SOLUTION, FAILED };
    outcome outcome;
    interval value;

    bool has_solution() const { return outcome == HAS_SOLUTION; }
    bool no_solution()  const { return outcome == NO_SOLUTION; }
};

// The f argument is either the beta_quantile or the ncbeta_noncent from above
newton_result solve_for_unique_solution(ddf_FctPtr f,
                                        interval   search_interval,
                                        cxsc::real accuracy)
{
    cxsc::ivector Zero;
    cxsc::intvector Unique;
    int NumberOfZeros, Error;
    AllZeros(f, search_interval, accuracy, Zero, Unique, NumberOfZeros, Error);
    newton_result s;
    s.outcome = newton_result::FAILED;
    if (!Error && NumberOfZeros==1 && Unique[1]) {
        s.outcome = newton_result::HAS_SOLUTION;
        s.value = Zero[1];
    }
    if (!Error && NumberOfZeros==0) {
        s.outcome = newton_result::NO_SOLUTION;
    }
    return s;
}

//-----------------------------------------------------------------------------
// Parsing one line of input. The expected format is:
//   a  b  x_0  lambda_0  alpha  beta  eps_x  eps_lambda
// where the elements are separated by arbitrary whitespace.
// The parameter ranges are also checked.

struct input {
    bool        error;
    std::string error_msg;
    interval    a;
    std::string a_str; // The string representation is used for the output
    int         b;
    interval    x0;
    interval    lambda0;
    interval    alpha;
    std::string alpha_str;
    interval    beta;
    std::string beta_str;
};

const int X_DIGITS = 12;
const cxsc::real TOL_X = std::pow(10.0, -X_DIGITS);

const int LAM_DIGITS = 10;
const cxsc::real TOL_LAMBDA = std::pow(10.0, -LAM_DIGITS);

std::string accuracy_x() {
    using namespace std;
    ostringstream oss;
    oss << scientific << setprecision(1) << cxsc::_double(TOL_X);
    return oss.str();
}

std::string accuracy_lambda() {
    using namespace std;
    ostringstream oss;
    oss << scientific << setprecision(1) << cxsc::_double(TOL_LAMBDA);
    return oss.str();
}

// trim whitespace from the end of s
void rtrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

interval to_interval(const std::string& s, const std::string& what) {
    using namespace std;
    istringstream iss(s);
    double dummy;
    iss >> dummy;
    char c;
    if (iss.fail() || iss.get(c)) {
        throw invalid_argument("Failed to parse \'"+what+"\' as double, got: \""+s+"\"");
    }
    string str = '[' + s + ',' + s + ']';
    interval x;
    str >> x;
    return x;
}

void check_ranges(input& i) {

    if (i.error)
        return;

    i.error = true;

    using cxsc::Inf;
    using cxsc::Sup;

    if (Inf(i.a) <= 0.0) {
        i.error_msg = "Failed, the shape parameter \'a\' must be strictly positive";
        return;
    }
    if (i.b < 1) {
        i.error_msg = "Failed, the integer shape parameter \'b\' must be at least one";
        return;
    }
    if (Inf(i.x0) <= 0.0 || Sup(i.x0) >= 1.0) {
        i.error_msg = "Failed, the inflated search interval \'x_0\' must be strictly between 0 and 1";
        return;
    }
    if (Inf(i.lambda0) <= 0.0) {
        i.error_msg = "Failed, the inflated search interval \'lambda_0\' must be strictly positive";
        return;
    }
    if (Inf(i.alpha) <= 0.0 || Sup(i.alpha) >= 1.0) {
        i.error_msg = "Failed, \'alpha\' must be strictly between 0 and 1";
        return;
    }
    if (Inf(i.beta) <= 0.0 || Sup(i.beta) >= 1.0) {
        i.error_msg = "Failed, \'beta\' must be strictly between 0 and 1";
        return;
    }

    i.error = false;
}

input parse_line(std::string line) {
    using namespace std;

    rtrim(line);

    // split the line to pieces
    istringstream iss(line);
    string a, alpha, beta;
    int b;
    double x0, lambda0, eps_x, eps_lambda;
    input i;
    i.error = true;
    iss >> a;
    if (iss.fail()) {
        i.error_msg = "Failed to parse the shape parameter \'a\' on line \""+line+"\"";
        return i;
    }
    iss >> b;
    char ws;
    iss.get(ws);
    if (iss.fail() || !isspace(ws)) {
        i.error_msg = "Failed to parse the shape parameter \'b\' as integer on line"
                      " \"" + line +"\"";
        return i;
    }
    iss >> x0 >> lambda0 >> alpha >> beta >> eps_x >> eps_lambda;
    char unused;
    if (iss.fail() || iss.get(unused)) {
        i.error_msg = "Failed to parse \""+line+"\", wrong format";
        return i;
    }

    // The ranges of eps_x and eps_lambda checked separately as we use them
    // in the construction of the inflated intervals
    if (eps_x < TOL_X*0.9999 || eps_x >= 1.0) {
        i.error_msg = "Failed, eps_x is expected to be between " + accuracy_x() +
                      " and 1";
        return i;
    }
    if (eps_lambda < TOL_LAMBDA*0.9999 || eps_lambda >= 1.0) {
        i.error_msg = "Failed, eps_lambda is expected to be between " + accuracy_lambda() +
                      " and 1";
        return i;
    }

    // convert the pieces to the corresponding type
    try {
        i.a      = to_interval(a, "a");
        i.b      = b;
        i.x0     = interval(x0*(1-eps_x), x0*(1+eps_x));
        i.lambda0= interval(lambda0*(1-eps_lambda), lambda0*(1+eps_lambda));
        i.alpha  = to_interval(alpha, "alpha");
        i.beta   = to_interval(beta, "beta");
        i.a_str     = a;
        i.alpha_str = alpha;
        i.beta_str  = beta;
        i.error = false;
    }
    catch(invalid_argument& e) {
        i.error_msg = e.what();
    }
    catch(cxsc::ERROR_INTERVAL_EMPTY_INTERVAL& ) {
        i.error_msg = "Failed to create the inflated interval for line \""+line+"\"";
    }
    catch(...) {
        i.error_msg = "Failed to parse \""+line+"\", unknown cause";
    }

    check_ranges(i);

    return i;
}

//-----------------------------------------------------------------------------

struct result {
    bool error;
    std::string error_msg;
    std::string x;
    std::string lambda;
};

// The result values are properly formatted strings
std::string format(double x, int digits) {
    using namespace std;
    ostringstream oss;
    oss << setprecision(digits) << scientific << x;
    return oss.str();
}

double midpoint(interval i) {
    return cxsc::_double(cxsc::mid(i));
}

std::string format(const std::string& name, interval a) {
    std::ostringstream oss;
    oss << "The search interval " << cxsc::SetPrecision(23,15)
        << cxsc::Scientific << a << " for \'" << name << "\' "
        << "is verified NOT to contain a zero";
    return oss.str();
}

result verify(interval a,
              int      b,
              interval x0,
              interval lambda0,
              interval alpha,
              interval beta)
{
    result res;
    res.error = true;

    quantile q(a, b, 1-alpha);
    q_ptr = &q; // save the state for the callback function beta_quantile

    newton_result s1 = solve_for_unique_solution(beta_quantile, x0, TOL_X);

    interval p_quantile;
    if (s1.has_solution())
        p_quantile = s1.value;
    else if (s1.no_solution()) {
        res.error_msg = format("x_0", x0);
        return res;
    }
    else {
        res.error_msg = "Failed to verify the quantile, try increasing eps_x";
        return res;
    }

    noncent_param np(p_quantile, a, b, beta);
    np_ptr = &np; // save the state for the callback function ncbeta_noncent

    newton_result s2 = solve_for_unique_solution(ncbeta_noncent,lambda0,TOL_LAMBDA);

    interval lambda;
    if (s2.has_solution())
        lambda = s2.value;
    else if (s2.no_solution()) {
        res.error_msg = format("lambda_0", lambda0);
        return res;
    }
    else {
        res.error_msg = "Failed to verify lambda, try increasing eps_lambda";
        return res;
    }
    res.x = format(midpoint(p_quantile), X_DIGITS);
    res.lambda = format(midpoint(lambda), LAM_DIGITS);
    res.error = false;
    return res;
}

//-----------------------------------------------------------------------------

std::string verify(const std::string& line) {

    static const std::string accu_x(accuracy_x()), accu_lambda(accuracy_lambda());

    input i = parse_line(line);

    if (i.error) {
       assert(!i.error_msg.empty());
       return i.error_msg + ", please check the documentation";
    }

    assert(i.error_msg.empty());

    result r = verify(i.a, i.b, i.x0, i.lambda0, i.alpha, i.beta);

    if (r.error) {
        assert(!r.error_msg.empty());
        return r.error_msg;
    }

    assert(r.error_msg.empty());

    std::ostringstream oss;
    oss << i.a_str << '\t' <<i.b << '\t' << r.x << '\t' << r.lambda << '\t' <<
            i.alpha_str << '\t' << i.beta_str << '\t' <<
            accu_x << '\t' << accu_lambda;

    return oss.str();
}

} // end of anonymous namespace

int main(int argc, char* argv[]) {

    using namespace std;

    if (argc > 1) {
        cout << "This program is distributed in the hope that it will be useful,\n"
                "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
                "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n"

                "This executable is an open source software. The source code is\n"
                "available at https://github.com/baharev/mindiffver under the \n"
                "GNU LGPL license.\n"
                "This program depends on C-XSC (GNU LGPL), see http://www.xsc.de/\n\n"

                "The program reads from the standard input, and writes to the\n"
                "standard output. If you have a table of the data to verify in\n"
                "a text file called input.txt, then you can run this program\n"
                "like this:\n\n";

        cout << "    " << argv[0] << " <input.txt >output.txt\n\n"

                "The verified results will be written to the output.txt file.\n";
        return 1;
    }

    // A hack: C-XSC logs to cerr when an exception is thrown, we "disable" it
    cerr.setstate(ios::failbit);

    const char* EOL =
#ifdef __CYGWIN__
                        "\r\n"; // Cygwin only write \n on Windows
#else
                        "\n";
#endif

    // Each line of input:
    //   a  b  x_0  lambda_0  alpha  beta  eps_x  eps_lambda
    // is mapped to the verified output:
    //   a  b  x   lambda   alpha  beta  eps_x  eps_lambda
    // where x and lambda has the relative accuracy TOL_X and TOL_LAMBDA, or a
    // single line error message is given if the parsing or the verification
    // fails.

    string line;
    while (getline(cin, line)) {
        try {
            cout << verify(line) << EOL << flush;
        }
        catch(...) {
            cout << "Failed with an unhandled exception, please report this bug" << EOL << flush;
        }
    }
}
