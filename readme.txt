===============================================================
C++ LP Solver Using Revised Dual Simplex Method
===============================================================
Tommy Zhong
tzrunm@gmail.com

---------------------------------------------------------------
How to Run
---------------------------------------------------------------
1.  Build source code in lp folder if needed.
2.  Run lp.exe with stdin input which is a standard form LP
    repesented by floating-point values separated by spaces.

    sample input format:
		c_1	c_2	c_3
		a_11	a_12	a_13	b_1
		a_21	a_22	a_23	b_2

    represents the following LP:
	max.	c_1x_1 + c_2x_2 + c_3x_3
	s.t.	a_11x_1 + a_12x_2 + a_13x_3 <= b_1
		a_21x_1 + a_22x_2 + a_23x_3 <= b_2
		x_1, x_2 >= 0

3.  Result would be printed in stdout. The first line is one
    of the following: optimal, infeasible, or unbounded. If
    optimal, the second line contains optimal objective value,
    and the third line contains optimal assignments for the
    variables separated by spaces.

---------------------------------------------------------------
Architecture
---------------------------------------------------------------
lp.cpp uses dense matrices of arbitrary-precise rational values
to implement revised dual simplex method. It is very slow, but
at least it outputs correct answers, but again it is slow.

This implementation is based on the pesudocode provided by UVic
linear programming course with some modifications.

The main function reads input and then call primal_simplex(), 
dual_simplex(), or both. A simplex function calls enter_blands()
or enter_largest_coeff() to find one pivoting variable, and then
either primal_leave() or dual_enter() to find another.

Pivoting rules include largest coefficient and Bland's.
The latter is only used whenever degeneracy appears to
avoid cycles while maintaining performance.

A subset of Boost library including Multiprecision and uBLAS is
used. This subset is extracted using the BCP tool and included
in the source code.

Linear systems are solved using LU factorization functions
provided by uBLAS library. These functions solve Ax = b
without explicitly creating A^-1 matrix.

---------------------------------------------------------------
Known problems
---------------------------------------------------------------
LU factorizaion shows poor performance when used with rational
values. I was just trying out the boost library and speed was
not a part of my consideration. This program is slow, and
testing with large inputs is not advised.