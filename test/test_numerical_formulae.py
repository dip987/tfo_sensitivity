import unittest
from tfo_sensitivity.jacobian.numerical_formulae import MuANumericalJC, OperatingPoint
from tfo_sensitivity.jacobian import PartialBloodJacobianMuAEqn, FullBloodJacobianMuAEqn


class MuANumericalJCTest(unittest.TestCase):
    def setUp(self) -> None:
        op = OperatingPoint(11.0, 1.0, 11.0, 0.5, 2)
        dx = "FS"
        mu_a_eqn = PartialBloodJacobianMuAEqn(0.2, 0.1, 0.75, 0.2, 0.1, 0.75)
        self.jc = MuANumericalJC(op, dx, mu_a_eqn)

    def test_runs_without_error(self):  # type: ignore
        try:
            print(self.jc.calculate_jacobian())
        except Exception as e:  # type: ignore
            self.fail(e)


if __name__ == "__main__":
    unittest.main()
