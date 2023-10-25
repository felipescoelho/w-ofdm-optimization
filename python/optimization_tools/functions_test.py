"""test_function.py

Unit test for functions. Based on examples from Antoniou, A. and Lu W.
"Practical Optimization: Algorithms and Engineering Applications", 2nd
Ed., Springer, 2021.

luizfelipe.coelho@smt.ufrj.br
Jul 31, 2023
"""


import numpy as np
import unittest
from optimization_core import (
    eliminate_equality, pd_pf_cqp, convert_inequality, nfi_pd_pf_cqp,
    p_sqp
)


class TestOptimizationAlgorithms(unittest.TestCase):

    def test_eliminate_equality(self):
        # QP problem:
        H = np.array(((1, 0, 0), (0, 1, 0), (0, 0, 0)), ndmin=2)
        p = np.array((2, 1, -1), ndmin=2).T
        A = np.array((0, 1, 1), ndmin=2)
        b = np.array((1), ndmin=2)
        # Constraint elimination:
        hat_H, hat_p = eliminate_equality(H, p, A, b)
        # Solve QP using hat_H and hat_p
        phi_opt = np.linalg.solve(hat_H, -hat_p)
        n_rows, n_cols = A.shape
        Q, hat_R = np.linalg.qr(A.T, mode='complete')
        R = hat_R[:n_rows, 0:n_rows]
        Q_1 = Q[0:n_cols, 0:n_rows]
        Q_2 = Q[0:n_cols, n_rows:n_cols]
        opt_sol = np.matmul(Q_2, phi_opt) + np.matmul(Q_1, np.matmul(
            np.linalg.inv(R).T, b
        ))
    
        self.assertTrue(np.allclose(
            opt_sol, np.array((-2., -2., 3.), ndmin=2).T
        ))

    def test_pd_pf_cqp(self):
        H = np.array(((4., 0., 0.), (0., 1., -1.), (0., -1., 1.)),
                     dtype=np.float64, ndmin=2)
        p = np.array((-8., -6., -6.), dtype=np.float64, ndmin=2).T
        A = np.array((1., 1., 1.), ndmin=2, dtype=np.float64)
        b = np.array((3.), ndmin=2, dtype=np.float64)
        x0 = np.array((1., 1., 1.), ndmin=2, dtype=np.float64).T
        lamb0 = np.array((7.), ndmin=2, dtype=np.float64)
        mu0 = np.array((3., 1., 1.), ndmin=2, dtype=np.float64).T
        rho = 3.
        x_sol, y_sol, n_iter = pd_pf_cqp(H, p, A, b, (x0, lamb0, mu0), rho)
        opt_xsol = np.array((.5, 1.25, 1.25), dtype=np.float64, ndmin=2).T
        opt_ysol = np.array((-18.5), dtype=np.float64, ndmin=2)
        
        self.assertTrue(np.allclose(y_sol, opt_ysol))
        self.assertTrue(np.allclose(x_sol, opt_xsol))

    # def test_nfi_pd_pf_cqp(self):
    #     H = np.array(((1., 0., -1., 0.), (0., 1., 0., -1.), (-1., 0., 1., 0.),
    #                   (0., -1., 0., 1.)), dtype=np.float64, ndmin=2)
    #     p = np.array((), ndmin=2)
    #     A = np.array(((-1., 0., 0., 0.), (0., -1., 0., 0.), (1., 2., 0., 0.),
    #                   (0., 0., 0., -1.), (0., 0., -1., -1.), (0., 0., 1., 2.)),
    #                  dtype=np.float64, ndmin=2)
    #     b = np.array((0., 0., 2., -2., -3., 6.), dtype=np.float64, ndmin=2).T
    #     hat_H, hat_p, hat_A = convert_inequality(H, p, A)
    #     x_sol, y_sol, n_iter = nfi_pd_pf_cqp(hat_H, hat_p, hat_A, b)
    #     x_sol_recovered = recover_variable(x_sol, 4)
    #     expected_x = np.array((.400002, .799999, 1.000001, 2.000003),
    #                           dtype=np.float64, ndmin=2).T
    #     self.assertTrue(np.allclose(x_sol_recovered, expected_x))
    #     self.assertAlmostEqual(1.341644, round(np.sqrt(2*y_sol)[0,0], 5), 5)

    def test_p_sqp(self):
        def f(x):
            f_val = np.log(1+x[0]**2) + x[1]
            a_val = (1 + x[0]**2)**2 + x[1]**2 - 4
            c_val = np.array([x[0] + x[1]**2 + .3], ndmin=2,
                             dtype=np.float64)

            return f_val, a_val, c_val
        
        def g(x):
            x = x.flatten()
            g_k = np.array([[2*x[0]/(1+x[0]**2)],
                            [1]], ndmin=2, dtype=np.float64)
            A_ek = np.array([4*x[0]*(1+x[0]**2), 2*x[1]], ndmin=2,
                            dtype=np.float64)
            A_ik = np.array([1, 2*x[1]], ndmin=2, dtype=np.float64)

            return g_k, A_ek, A_ik
        
        x0 = np.array([[-1.5], [1]], ndmin=2, dtype=np.float64)
        x_sol, _ = p_sqp(f, g, x0)

        print(x_sol)


        
if __name__ == '__main__':
    # TestOptimizationAlgorithms()
    unittest.main()


# EoF
