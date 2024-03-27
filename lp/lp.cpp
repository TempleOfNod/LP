/*
    Revised Dual Simplex Method
    Tommy Zhong
    V00928055
*/

#include <iostream>
#include <string>
#include <functional>
#include <vector>
#include <algorithm>

#define BOOST_UBLAS_TYPE_CHECK 0                    // float type checks aren't appliable here
// #include <boost/rational.hpp>                    // fast, but prone to overflow during LU factorization
#include <boost/multiprecision/cpp_int.hpp>         // safe choice but slow as hell
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>

//typedef boost::rational<long long> rational;
typedef boost::multiprecision::cpp_rational rational;
typedef boost::numeric::ublas::matrix<rational> matrix_r;
typedef boost::numeric::ublas::vector<rational> vector_r;
typedef boost::numeric::ublas::matrix_column<matrix_r> column_r;

// finds an entering variable by Bland's rule
// returns j where xN[j] is chosen
// if no variable can be chosen i.e., LP is optimal, return -1
// for dual leaving variable, replace (zN, N) with (xB, B)
inline int enter_blands(const vector_r &zN, const std::vector<int> &N)
{
    int min = std::numeric_limits<int>::max();
    int enter = -1;

    for (unsigned j = 0; j < zN.size(); j++)
    {
        if (zN(j) < 0)
        {
            if (N[j] < min)
            {
                min = N[j];
                enter = j;
            }
        }
    }
    return enter;
}

// finds an entering variable by largest coefficient rule
// (really should be called "smallest coeff rule" since sign is flipped)
// returns j where xN[j] is chosen
// if no variable can be chosen i.e., LP is optimal, return -1
// for dual leaving variable, replace zN with xB
inline int enter_largest_coeff(const vector_r &zN)
{
    rational coeff = 0;
    int enter = -1;

    for (unsigned j = 0; j < zN.size(); j++)
    {
        if (zN(j) < coeff)
        {
            coeff = zN(j);
            enter = j;
        }
    }
    return enter;
}

// finds a primal leaving variable given an entering variable
// returns i where xB[i] is chosen
// if no variable can be chosen i.e., LP is unbounded, return -1
int primal_leave(const matrix_r &AB, const matrix_r &AN, vector_r &xB,
        vector_r &xN, const std::vector<int> &B, const int enter)
{
    vector_r dxB(xB.size());
    rational t = -1;
    rational temp;
    int min = std::numeric_limits<int>::max();
    int leave = -1;
    
    // calculate delta xB = AB^-1 Aj
    for (unsigned i = 0; i < xB.size(); i++)
    {
        dxB(i) = AN(i, enter);
    }
    matrix_r AB_copy = AB;
    boost::numeric::ublas::permutation_matrix<size_t> pm(AB_copy.size1());
    boost::numeric::ublas::lu_factorize(AB_copy, pm);
    boost::numeric::ublas::lu_substitute(AB_copy, pm, dxB);

    // find maximum increase t
    for (unsigned i = 0; i < xB.size(); i++)
    {
        if (dxB(i) > 0)
        {
            temp = xB(i) / dxB(i);
            if (t < 0 || temp < t)
            {
                t = temp;
                leave = 0;
            }
        }
    }

    // choose leaving variable by Bland's rule
    if (!leave)
    {
        for (unsigned i = 0; i < xB.size(); i++)
        {
            if (dxB(i) > 0 && xB(i) / dxB(i) == t && B[i] < min)
            {
                min = B[i];
                leave = i;
            }
        }

        // update xB and xN
        xB -= dxB * t;
        xN(enter) = t;
    }
    return leave;
}

// finds a dual entering variable given a leaving variable
// returns j where zN[j] is chosen
// if no variable can be chosen i.e., LP is unbounded, return -1
int dual_enter(const matrix_r &AB, const matrix_r &AN, vector_r &zB,
        vector_r &zN, const std::vector<int> &N, const int leave)
{
    vector_r dzN(zN.size());
    vector_r v(zB.size());
    rational s = -1;
    rational temp;
    int min = std::numeric_limits<int>::max();
    int enter = -1;
    
    // calculate delta zN = -AN^T v where v = (AB^-1)^T u where u[leave] = 1
    v(leave) = 1;   // leave is an index of B, already what u needs
    matrix_r ABT = boost::numeric::ublas::trans(AB);
    boost::numeric::ublas::permutation_matrix<size_t> pmT(ABT.size1());
    boost::numeric::ublas::lu_factorize(ABT, pmT);
    boost::numeric::ublas::lu_substitute(ABT, pmT, v);

    dzN = boost::numeric::ublas::prod(-boost::numeric::ublas::trans(AN), v);

    // find maximum increase s
    for (unsigned j = 0; j < zN.size(); j++)
    {
        if (dzN(j) > 0)
        {
            temp = zN(j) / dzN(j);
            if (s < 0 || temp < s)
            {
                s = temp;
                enter = 0;
            }
        }
    }

    // choose entering variable by Bland's rule
    if (!enter)
    {
        for (unsigned j = 0; j < zN.size(); j++)
        {
            if (dzN(j) > 0 && zN(j) / dzN(j) == s && N[j] < min)
            {
                min = N[j];
                enter = j;
            }
        }

        // update zN and zB
        zN -= dzN * s;
        zB(leave) = s;
    }
    return enter;
}

// simplex method on primal LP
// returns optimum when optimal, -1 when unbounded, -2 when infeasible
rational primal_simplex(matrix_r &AB, matrix_r &AN, const vector_r &b, 
        vector_r &cB, vector_r &cN, std::vector<int> &B, std::vector<int> &N, bool print)
{
    vector_r xB;
    vector_r xN(N.size());
    vector_r zN(N.size());
    int enter;
    int leave;
    bool degen = false;

    // compute xB = (AB)^-1 b using LU factorization
    // use a copy to avoid changing AB
    // thanks to https://stackoverflow.com/questions/1225411/ for providing
    // an example that Boost documentation should have provided
    {
        xB = b;
        matrix_r AB_copy = AB;
        boost::numeric::ublas::permutation_matrix<size_t> pm(AB.size1());
        boost::numeric::ublas::lu_factorize(AB_copy, pm);
        boost::numeric::ublas::lu_substitute(AB_copy, pm, xB);
    }
    
    // pivot loop
    for (;;)
    {   
        // compute zN = AN^T v - cN where v = (AB^-1)^T cB    
        vector_r v = cB;
        matrix_r ABT = boost::numeric::ublas::trans(AB);
        boost::numeric::ublas::permutation_matrix<size_t> pmT(ABT.size1());
        boost::numeric::ublas::lu_factorize(ABT, pmT);
        boost::numeric::ublas::lu_substitute(ABT, pmT, v);
        zN = boost::numeric::ublas::prod(boost::numeric::ublas::trans(AN), v) - cN;
        
        // check degeneracy to decide pivoting rule
        degen = false;
        for (const auto &r : xB)
        {
            if (r < 0)
            {
                // error: pivot reaches infeasible point
                // variable overflow can cause this, but might not be the only factor
                std::cerr << "primal simplex reached infeasible point" << std::endl;
                for (auto &xb : xB) std::cerr << xb << std::endl;
                return -2;   // infeasible
            }
            else if (r == 0)
            {
                degen = true;   // degenerate
            }
        }

        // find leaving variable and check if zN is optimal
        // use Bland's rule if degeneracy is found
        // use largest coefficient rule otherwise
        enter = degen? enter_blands(zN, N) : enter_largest_coeff(zN);

        if (enter < 0)
        {
            rational optimum = boost::numeric::ublas::inner_prod(cB, xB);
            if (print)
            {
                double dbuf = optimum.convert_to<double>();
                // double dbuf = boost::rational_cast<double>(optimum);

                // print optimum and x
                // need to search from B which is unsorted
                printf("optimal\n%.7g\n", dbuf);
                for (unsigned i = 0; i < N.size(); i++)
                {
                    std::vector<int>::iterator pos = std::find(B.begin(), B.end(), i);
                    if (i > 0)
                    {
                        printf(" ");
                    }
                    if (pos != B.end())
                    {
                        dbuf = xB(pos - B.begin()).convert_to<double>();
                        // dbuf = boost::rational_cast<double>(xB(pos - B.begin()));
                        printf("%.7g", dbuf);
                    }
                    else
                    {
                        printf("0");
                    }
                }
                std::cout << std::endl;
            }
            return optimum;     // optimal
        }
        
        // find leaving variable
        leave = primal_leave(AB, AN, xB, xN, B, enter);

        if (leave < 0)
        {
            if (print)
            {
                std::cout << "unbounded" << std::endl;
            }
            return -1;   // unbounded
        }

        // pivot by swapping values
        int temp = B[leave];
        B[leave] = N[enter];
        N[enter] = temp;

        rational temp_r = xB(leave);
        xB(leave) = xN(enter);
        xN(enter) = temp_r;

        temp_r = cB(leave);
        cB(leave) = cN(enter);
        cN(enter) = temp_r;

        column_r colB(AB, leave);
        column_r colN(AN, enter);
        colB.swap(colN);
    }
}

// simplex method on dual LP
// returns optimum when optimal, -1 when unbounded, -2 when infeasible
rational dual_simplex(matrix_r &AB, matrix_r &AN, const vector_r &b, 
        vector_r &cB, vector_r &cN, std::vector<int> &B, std::vector<int> &N, bool print)
{
    vector_r xB(B.size());
    vector_r zB(B.size());
    vector_r zN(N.size());
    int enter;
    int leave;
    bool degen = false;

    // compute zN = AN^T v - cN where v = (AB^-1)^T cB 
    {   
        vector_r v = cB;
        matrix_r ABT = boost::numeric::ublas::trans(AB);
        boost::numeric::ublas::permutation_matrix<size_t> pmT(ABT.size1());
        boost::numeric::ublas::lu_factorize(ABT, pmT);
        boost::numeric::ublas::lu_substitute(ABT, pmT, v);
        zN = boost::numeric::ublas::prod(boost::numeric::ublas::trans(AN), v) - cN;
    }
    
    // pivot loop
    for (;;)
    {   
        // compute xB = (AB)^-1 b
        xB = b;
        matrix_r AB_copy = AB;
        boost::numeric::ublas::permutation_matrix<size_t> pm(AB.size1());
        boost::numeric::ublas::lu_factorize(AB_copy, pm);
        boost::numeric::ublas::lu_substitute(AB_copy, pm, xB);
        
        // check degeneracy to decide pivoting rule
        degen = false;
        for (const auto &r : zN)
        {
            if (r < 0)
            {
                // error: pivot reaches infeasible point
                std::cerr << "dual simplex reached infeasible point" << std::endl;
                for (auto &zn : zN) std::cerr << zn << std::endl;
                return -2;   // infeasible
            }
            else if (r == 0)
            {
                degen = true;   // degenerate
            }
        }

        // find leaving variable and check if xB is optimal
        // use Bland's rule if degeneracy is found
        // use largest coefficient rule otherwise
        leave = degen? enter_blands(xB, B) : enter_largest_coeff(xB);

        if (leave < 0)
        {
            rational optimum = boost::numeric::ublas::inner_prod(cB, xB);
            if (print)
            {
                double dbuf = optimum.convert_to<double>();
                // double dbuf = boost::rational_cast<double>(optimum);

                // print optimum and x
                // need to search from B which is unsorted
                printf("optimal\n%.7g\n", dbuf);
                for (unsigned i = 0; i < N.size(); i++)
                {
                    std::vector<int>::iterator pos = std::find(B.begin(), B.end(), i);
                    if (i > 0)
                    {
                        printf(" ");
                    }
                    if (pos != B.end())
                    {
                        dbuf = xB(pos - B.begin()).convert_to<double>();
                        // dbuf = boost::rational_cast<double>(xB(pos - B.begin()));
                        printf("%.7g", dbuf);
                    }
                    else
                    {
                        printf("0");
                    }
                }
                std::cout << std::endl;
            }
            return optimum;     // optimal
        }
        
        // find leaving variable
        enter = dual_enter(AB, AN, zB, zN, N, leave);

        if (enter < 0)
        {
            if (print)
            {
                std::cout << "infeasible" << std::endl;
            }
            return -1;   // unbounded (primal infeasible)
        }

        // pivot by swapping values
        int temp = B[leave];
        B[leave] = N[enter];
        N[enter] = temp;

        rational temp_r = zB(leave);
        zB(leave) = zN(enter);
        zN(enter) = temp_r;

        temp_r = cB(leave);
        cB(leave) = cN(enter);
        cN(enter) = temp_r;

        column_r colB(AB, leave);
        column_r colN(AN, enter);
        colB.swap(colN);
    }
}

int main()
{
    std::vector<double> inputs;
    double dbuf;
    int bsize = -1;
    int nsize = 0;
    int count = 0;
    bool feasible_p = false;
    bool feasible_d = false;

    // find number of non-basic variables
    while (std::cin.peek() != '\n')
    {
        std::cin >> dbuf;
        inputs.push_back(dbuf);
        nsize++;
    }

    // find number of basic variables
    {
        std::string strbuf;
        while (std::getline(std::cin, strbuf))  // skip a line
        {
            bsize++;    // start at -1 because the loop reads the first line
            if (std::cin.peek() != EOF)
            {
                for (int i = 0; i < nsize + 1; i++)
                {
                    std::cin >> dbuf;
                    inputs.push_back(dbuf);
                }
            }

        }
    }

    // initialize matrices
    matrix_r AB = boost::numeric::ublas::identity_matrix<rational>(bsize);
    matrix_r AN(bsize, nsize);
    vector_r b(bsize);
    vector_r cB(bsize);
    vector_r cN(nsize);

    // initialize column indices
    std::vector<int> B(bsize);
    std::vector<int> N(nsize); 
    
    // read objective function
    for (int j = 0; j < nsize; j++)
    {
        cN(j) = inputs[count++];
        // cN(j) = rational((long long)(inputs[count++] * 10000000LL), 10000000LL);
        N[j] = j;
    }

    // read constraints
    for (int i = 0; i < bsize; i++)
    {
        for (int j = 0; j < nsize; j++)
        {
            AN(i, j) = inputs[count++];
            // AN(i, j) = rational((long long)(inputs[count++] * 10000000LL), 10000000LL);
        }
        b(i) = inputs[count++];
        // b(i) = rational((long long)(inputs[count++] * 10000000LL), 10000000LL);
        B[i] = i + nsize;
    }

    // check primal feasibility
    feasible_p = true;
    for (auto &i : b)
    {
        if (i < 0)
        {
            feasible_p = false;
            // break;
        }
    }
    
    if (feasible_p)
    {
        // run primal simplex method
        primal_simplex(AB, AN, b, cB, cN, B, N, true);
    }
    else
    {
        // check dual feasibility
        feasible_d = true;
        for (auto &j : cN)
        {
            if (j > 0)
            {
                feasible_d = false;
                break;
            }
        }

        if (feasible_d)
        {
            // run dual simplex method
            dual_simplex(AB, AN, b, cB, cN, B, N, true);
        }
        else
        {
            // call dual simplex with c = 0 to find a feasible point
            vector_r zeroB(bsize);
            vector_r zeroN(nsize);
            auto ret = dual_simplex(AB, AN, b, zeroB, zeroN, B, N, false);

            if (!ret)
            {
                // feasible point found
                // update original cB, cN to correspond new B and N
                zeroB = cB;     // reused as temp storage
                zeroN = cN;
                for (int i = 0; i < bsize; i++)
                {
                    if (B[i] >= nsize)
                    {
                        cB(i) = zeroB(B[i] - nsize);
                    }
                    else
                    {
                        cB(i) = zeroN(B[i]);
                    }
                }

                for (int j = 0; j < nsize; j++)
                {
                    if (N[j] >= nsize)
                    {
                        cN(j) = zeroB(N[j] - nsize);
                    }
                    else
                    {
                        cN(j) = zeroN(N[j]);
                    }
                }

                // run simplex with new B and N
                primal_simplex(AB, AN, b, cB, cN, B, N, true);
            }
            else if (ret == -1)
            {
                // no feasible point found
                std::cout << "infeasible" << std::endl;
            }
            else
            {
                // unexpected failure
                std::cerr << "failed dual simplex with c = 0" << std::endl;
            }
        }
    }

    return 0;
}