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
    var_trans_equality, recover_var_eq_trans, pd_pf_cqp, nfi_pd_pf_cqp,
    slack_var_inequality, recover_var_inequality, pd_pf_sdp, pc_sdp
)


class TestOptimizationAlgorithms(unittest.TestCase):

    def test_example_13_1(self):
        """
        Example 13.1 from Antoniou, A. and Lu W. "Practical
        Optimization: Algorithms and Engineering Applications", 2nd Ed.,
        Springer, 2021.
        """

        # QP problem:
        H = np.array(((1, 0, 0), (0, 1, 0), (0, 0, 0)), ndmin=2)
        p = np.array((2, 1, -1), ndmin=2).T
        A = np.array((0, 1, 1), ndmin=2)
        b = np.array((1), ndmin=2)
        # Computing:
        hat_H, hat_p = var_trans_equality(H, p, A, b)
        # Solve QP using hat_H and hat_p
        phi_opt = np.linalg.solve(hat_H, -hat_p)
        x_sol = recover_var_eq_trans(phi_opt, A, b)

        self.assertTrue(np.allclose(
            x_sol, np.array([[-2.], [-2.], [3.]], ndmin=2, dtype=np.float64)
        ))
        print('Example 13.1: OK.')

    def test_example_13_3(self):
        """
        Example 13.3 using Algorithm 13.2 and Algorithm 13.3 from
        Antoniou, A. and Lu W. "Practical Optimization: Algorithms and
        Engineering Applications", 2nd Ed., Springer, 2021.
        """

        # Problem definition:
        H = np.array([[4., 0., 0.], [0., 1., -1.], [0., -1., 1.]],
                     dtype=np.float64, ndmin=2)
        p = np.array([[-8.], [-6.], [-6.]], dtype=np.float64, ndmin=2)
        A = np.array([[1., 1., 1.]], ndmin=2, dtype=np.float64)
        b = np.array([[3.]], ndmin=2, dtype=np.float64)
        # Solve for Alg. 13.2
        x0_0 = np.array([[1.], [1.], [1.]], ndmin=2, dtype=np.float64)
        lamb0_0 = np.array([[7.]], ndmin=2, dtype=np.float64)
        mu0_0 = np.array([[3.], [1.], [1.]], ndmin=2, dtype=np.float64)
        x_sol_0, _ = pd_pf_cqp(H, p, A, b, (x0_0, lamb0_0, mu0_0))
        # Solve for Alg. 13.3
        x0_1 = np.array([[1], [2], [2]], ndmin=2, dtype=np.float64)
        lamb0_1 = np.array([[1]], ndmin=2, dtype=np.float64)
        mu0_1 = np.array([[.2], [.2], [.2]], ndmin=2, dtype=np.float64)
        rho_1 = 3+2*np.sqrt(3)
        x_sol_1, _ = nfi_pd_pf_cqp(H, p, A, b, (x0_1, lamb0_1, mu0_1), rho_1,
                                   epsilon=1e-5)
        x_exp = np.array([[.5], [1.25], [1.25]], dtype=np.float64, ndmin=2)

        self.assertTrue(np.allclose(x_sol_1, x_exp))
        self.assertTrue(np.allclose(x_sol_0, x_exp))
        print('Example 13.3: OK.')

    def test_example_13_4(self):
        """
        Example 13.4 from Antoniou, A. and Lu W. "Practical
        Optimization: Algorithms and Engineering Applications", 2nd Ed.,
        Springer, 2021.
        """

        # Definitions:
        n = 4
        H = np.array([[1., 0., -1., 0.],
                      [0., 1., 0., -1.],
                      [-1., 0., 1., 0.],
                      [0., -1., 0., 1.]], ndmin=2, dtype=np.float64)
        p = np.zeros((n, 1), dtype=np.float64)
        A = np.array([[-1., 0., 0., 0.],
                      [0., -1., 0., 0.],
                      [1., 2., 0., 0.],
                      [0., 0., 0., -1.],
                      [0., 0., -1., -1.],
                      [0., 0., 1., 2.]], ndmin=2, dtype=np.float64)
        b = np.array([[0.], [0.], [2.], [-2.], [-3.], [6.]], ndmin=2,
                     dtype=np.float64)
        hat_H, hat_p, hat_A = slack_var_inequality(H, p, A)
        x0 = np.ones((14, 1), dtype=np.float64)
        lamb0 = np.ones((6, 1), dtype=np.float64)
        mu0 = np.ones((14, 1), dtype=np.float64)
        rho = 14 + 20*np.sqrt(14)
        epsilon = 1e-5
        x_sol_hat, _ = nfi_pd_pf_cqp(hat_H, hat_p, hat_A, b, (x0, lamb0, mu0),
                                     rho, epsilon)
        x_sol = recover_var_inequality(x_sol_hat, n)
        x_exp = np.array([[.400002], [.799999], [1.000001], [2.000003]],
                         ndmin=2, dtype=np.float64)
        y_sol = .5* x_sol.T @ H @ x_sol + x_sol.T@ p
        min_dist = np.sqrt(2*y_sol.flatten())[0]
        dist_exp = 1.341644

        self.assertTrue(np.allclose(x_sol, x_exp, rtol=1e-5, atol=1e-5))
        self.assertAlmostEqual(dist_exp, min_dist, places=5)
        print('Example 13.4: OK.')

    def test_example_13_5(self):
        """
        Example 13.5 from Antoniou, A. and Lu W. "Practical
        Optimization: Algorithms and Engineering Applications", 2nd Ed.,
        Springer, 2021.
        """

        # Definitions:
        A0 = np.array([[2, -.5, -.6],
                       [-.5, 2, .4],
                       [-.6, .4, 3]], ndmin=2, dtype=np.float64)
        A1 = np.array([[0, 1, 0],
                       [1, 0, 0],
                       [0, 0, 0]], ndmin=2, dtype=np.float64)
        A2 = np.array([[0, 0, 1],
                       [0, 0, 0],
                       [1, 0, 0]], ndmin=2, dtype=np.float64)
        A3 = np.array([[0, 0, 0],
                       [0, 0, 1],
                       [0, 1, 0]], ndmin=2, dtype=np.float64)
        A4 = np.eye(3)
        b = np.array([[0], [0], [0], [1]], ndmin=2, dtype=np.float64)
        C = -A0
        X0 = np.eye(3)/3
        y0 = np.array([[-.2], [-.2], [-.2], [4]], ndmin=2, dtype=np.float64)
        S0 = np.array([[2, .3, .4],
                       [.3, 2, -.6],
                       [.4, -.6, 1]], ndmin=2, dtype=np.float64)
        epsilon = 1e-3
        X, y, S, n_iter = pd_pf_sdp([A1, A2, A3, A4], b, C, (X0, y0, S0),
                                    .9, None, epsilon)
        y_exp = np.array([[-.392921], [-.599995], [.399992], [3.000469]],
                         dtype=np.float64)
        
        self.assertTrue(np.allclose(y, y_exp, rtol=1e-5, atol=1e-5))
        print('Example 13.5: OK.')

    def test_example_13_6(self):
        """
        Example 13.6 from Antoniou, A. and Lu W. "Practical
        Optimization: Algorithms and Engineering Applications", 2nd Ed.,
        Springer, 2021.
        """

        def diag(x0, x1, x2):
            n_row0, n_col0 = x0.shape
            n_row1, n_col1 = x1.shape
            n_row2, n_col2 = x2.shape

            y = np.vstack((
                np.hstack((x0, np.zeros((n_row0, n_col1+n_col2)))),
                np.hstack((np.zeros((n_row1, n_col0)), x1, np.zeros((n_row1, n_col2)))),
                np.hstack((np.zeros((n_row2, n_col0+n_col1)), x2))
            ))

            return y
        
        # Define problem as SDP:
        q1 = np.array([[-1/2], [0]], ndmin=2, dtype=np.float64)
        q2 = np.array([[-44], [-52]], ndmin=2, dtype=np.float64)
        H_hat = np.array([[1, 0, -1, 0], [0, -1, 0, 1], [0, 0, 0, 0],
                          [0, 0, 0, 0]], ndmin=2, dtype=np.float64)
        Q1_hat = np.array([[1/2, 0], [0, 1]], ndmin=2, dtype=np.float64)
        Q2_hat = np.array([[2, 2], [-1, 1]], ndmin=2, dtype=np.float64)
        b = np.array([[0], [0], [0], [0], [1]], ndmin=2, dtype=np.float64)
        A1 = -diag(
            np.vstack((
                np.hstack((np.zeros((4, 4)), H_hat[:, 0:1])),
                np.hstack((H_hat[:, 0:1].T, np.array([[0]], ndmin=2,
                                                      dtype=np.float64)))
            )),
            np.vstack((
                np.hstack((np.zeros((2, 2)), Q1_hat[:, 0:1])),
                np.hstack((Q1_hat[:, 0:1].T, -q1[0:1, 0:1]))
            )),
            np.zeros((3, 3))
        )
        A2 = -diag(
            np.vstack((
                np.hstack((np.zeros((4, 4)), H_hat[:, 1:2])),
                np.hstack((H_hat[:, 1:2].T, np.array([[0]], ndmin=2,
                                                     dtype=np.float64)))
            )),
            np.vstack((
                np.hstack((np.zeros((2, 2)), Q1_hat[:, 1:2])),
                np.hstack((Q1_hat[:, 1:2].T, -q1[1:2, 0:1]))
            )),
            np.zeros((3, 3))
        )
        A3 = -diag(
            np.vstack((
                np.hstack((np.zeros((4, 4)), H_hat[:, 2:3])),
                np.hstack((H_hat[:, 2:3].T, np.array([[0]], ndmin=2,
                                                     dtype=np.float64)))
            )),
            np.zeros((3, 3)),
            np.vstack((
                np.hstack((np.zeros((2, 2)), Q2_hat[:, 0:1])),
                np.hstack((Q2_hat[:, 0:1].T, -q2[0:1, 0:1]))
            ))
        )
        A4 = -diag(
            np.vstack((
                np.hstack((np.zeros((4, 4)), H_hat[:, 3:4])),
                np.hstack((H_hat[:, 3:4].T, np.array([[0]], ndmin=2,
                                                     dtype=np.float64)))
            )),
            np.zeros((3, 3)),
            np.vstack((
                np.hstack((np.zeros((2, 2)), Q2_hat[:, 1:2])),
                np.hstack((Q2_hat[:, 1:2].T, -q2[1:2, 0:1]))
            ))
        )
        A5 = -diag(
            np.vstack((
                np.hstack((np.zeros((4, 4)), np.zeros((4, 1)))),
                np.hstack((np.zeros((1, 4)), np.array([[1]], ndmin=2,
                                                      dtype=np.float64)))
            )), np.zeros((3, 3)), np.zeros((3, 3))
        )
        C = -diag(
            np.vstack((
                np.hstack((np.eye(4), np.zeros((4, 1)))),
                np.zeros((1, 5))
            )),
            np.vstack((
                np.hstack((np.eye(2), np.zeros((2, 1)))),
                np.hstack((np.zeros((1, 2)), np.array([[3/4]], ndmin=2,
                                                      dtype=np.float64)))
            )),
            np.vstack((
                np.hstack((np.eye(2), np.zeros((2, 1)))),
                np.hstack((np.zeros((1, 2)), np.array([[-140]], ndmin=2,
                                                      dtype=np.float64)))
            ))
        )
        X02 = np.array([[1, 0, -.5], [0, 1, 0], [-.5, 0, 1]], ndmin=2,
                       dtype=np.float64)
        X03 = np.array([[180, 0, -12], [0, 60, -2], [-12, -2, 1]], ndmin=2,
                       dtype=np.float64)
        X0 = diag(np.eye(5), X02, X03)
        y0 = np.array([[1], [0], [2], [4], [20]], ndmin=2, dtype=np.float64)
        A = [A1, A2, A3, A4, A5]
        S0 = C - sum([y0[i]*A[i] for i in range(5)])
        
        X, y, S, n_iter = pc_sdp(A, b, C, (X0, y0, S0), gamma=.9, epsilon=1e-3)
        print(y)

    def _example_15_1(self):
        """
        Example 15.1 from Antoniou, A. and Lu W. "Practical
        Optimization: Algorithms and Engineering Applications", 2nd Ed.,
        Springer, 2021.
        """

        # Definitions
        def f(x):
            x = x.flatten()
            f_val = np.array([[np.log(1+x[0]**2) + x[1]]], ndmin=2,
                             dtype=np.float64)
            a_val = np.array([[(1 + x[0]**2)**2 + x[1]**2 - 4]], ndmin=2,
                             dtype=np.float64)
            c_val = np.array([[x[0] + x[1]**2 + .3]], ndmin=2,
                             dtype=np.float64)

            return f_val, a_val, c_val
        
        def g(x):
            x = x.flatten()
            f_val = np.array([[2*x[0]/(1+x[1]**2)], [1]], ndmin=2,
                             dtype=np.float64)
            a_val = np.array([[4*x[0]*(1+x[0]**2), 2*x[1]]], ndmin=2,
                             dtype=np.float64)
            c_val = np.array([[1, 2*x[1]]], ndmin=2, dtype=np.float64)
            
            return f_val, a_val, c_val
        
        def h(x):
            x = x.flatten()
            f_val = np.array([[np.max((0, 2*(1-x[0]**2)/((1+x[0]**2)**2))), 0],
                              [0, 0]], ndmin=2, dtype=np.float64)
            a_val = np.array([[4+12*x[0]**2, 0], [0, 2]], ndmin=2,
                             dtype=np.float64)
            c_val = np.array([[0, 0], [0, 2]], ndmin=2, dtype=np.float64)

            return f_val, a_val, c_val
        
        x0 = np.array([[-1.5], [1]], ndmin=2, dtype=np.float64)
        rho = np.array([[.4], [.4]], ndmin=2, dtype=np.float64)
        epsilon = 1e-6
        tau = 1
        rho = 1
        x_sol, _ = sqp_basic(f, g, h, x0, 0, 0, tau, epsilon)
        print(x_sol)
        # f_sol, a_sol, c_sol = f(x_sol)
        # f_exp = -1.755174*1e-1
        # a_exp = 1e-15
        # c_exp = -2.060264*1e-9

        # self.assertAlmostEqual(f_sol, f_exp)
        # self.assertLess(a_sol, a_exp)
        print('Example 15.1: OK.')
        

    # def test_p_sqp(self):
    #     def f(x):
    #         f_val = np.log(1+x[0]**2) + x[1]
    #         a_val = (1 + x[0]**2)**2 + x[1]**2 - 4
    #         c_val = x[0] + x[1]**2 + .3

    #         return f_val, a_val, c_val
        
    #     def g(x):
    #         x = x.flatten()
    #         g_k = np.array([[2*x[0]/(1+x[0]**2)],
    #                         [1]], ndmin=2, dtype=np.float64)
    #         A_ek = np.array([4*x[0]*(1+x[0]**2), 2*x[1]], ndmin=2,
    #                         dtype=np.float64)
    #         A_ik = np.array([1, 2*x[1]], ndmin=2, dtype=np.float64)

    #         return g_k, A_ek, A_ik
        
    #     x0 = np.array([[-.9162420], [-.7850108]], ndmin=2, dtype=np.float64)

    #     x_sol, _ = sqp_1e(f, g, x0)
    #     x_exp = np.array([[-.9162420], [-.7850108]], ndmin=2, dtype=np.float64)

    #     print(x_sol)
    #     print(f(x_sol))


        
if __name__ == '__main__':
    # TestOptimizationAlgorithms()
    unittest.main()


# EoF
