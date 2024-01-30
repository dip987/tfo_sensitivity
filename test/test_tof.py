from pathlib import Path
import unittest
import pandas as pd
import numpy as np
from tfo_sensitivity.tof.base import ToF
from tfo_sensitivity.tof.optimization import two_pointer_discrete_optimize, two_pointer_brute_force_optimize

MU_MAP = {1: 0.1, 2: 0.2, 3: 0.3, 4: 0.4}
TIME_RES = 2e-13


class ToFTest(unittest.TestCase):
    def setUp(self):
        test_data_path = Path(__file__).parent.resolve() / "raw_data_test.pkl"
        self.data = pd.read_pickle(test_data_path)
        self.tof = ToF(self.data, TIME_RES, MU_MAP, None)

    def test_ToF_data_non_zero_length(self):
        self.assertTrue(len(self.tof.data) > 0)

    def test_ToF_interpolation(self):
        index_to_drop = self.tof.data.index.values[2]
        # remove a point
        self.tof.data.drop(index_to_drop, inplace=True)
        self.tof._fill_data_gaps()  # type: ignore
        self.assertTrue(index_to_drop in self.tof.data.index.values)
        self.assertTrue(self.tof.interpolation_needed)

    def test_ToF_create_class_from_data(self):
        index_array = pd.RangeIndex(2, 5, name="ToF_quantized")
        data = pd.Series([1e-13, 2e-15, 1e-20], index_array, name="Intensity")
        tof = ToF.from_data(data, 5e-9)
        self.assertTrue(tof.data.equals(data))

    def test_ToF_truediv(self):
        index_array = pd.RangeIndex(2, 5, name="ToF_quantized")
        data1 = pd.Series([1e-13, 2e-15, 1e-20], index_array, name="Intensity")
        data2 = pd.Series([1e-13, 2e-15, 1e-20], index_array, name="Intensity")
        tof2 = ToF.from_data(data2, 3e-9)
        tof1 = ToF.from_data(data1, 3e-9)
        tof3 = tof1 / tof2
        self.assertTrue(tof3.data.equals(pd.Series([1.0, 1.0, 1.0], index_array, name="Intensity")))

    def test_ToF_mul(self):
        index_array = pd.RangeIndex(2, 5, name="ToF_quantized")
        data1 = pd.Series([1e-13, 2e-15, 1e-20], index_array, name="Intensity")
        data2 = pd.Series([1e-13, 2e-15, 1e-20], index_array, name="Intensity")
        tof2 = ToF.from_data(data2, 3e-9)
        tof1 = ToF.from_data(data1, 3e-9)
        tof3 = tof1 * tof2
        self.assertTrue(tof3.data.equals(pd.Series([1e-26, 4e-30, 1e-40], index_array, name="Intensity")))

    def test_two_ToF_operation_compatibility(self):
        index_array = pd.RangeIndex(2, 5, name="ToF_quantized")
        index_array2 = pd.RangeIndex(6, 9, name="ToF_quantized")
        data1 = pd.Series([1e-13, 2e-15, 1e-20], index_array, name="Intensity")
        data2 = pd.Series([1e-13, 2e-15, 1e-20], index_array2, name="Intensity")
        tof1 = ToF.from_data(data1, 3e-9)
        tof2 = ToF.from_data(data2, 3e-9)
        self.assertTrue(tof1.check_operation_compatibility(tof1))
        self.assertFalse(tof1.check_operation_compatibility(tof2))


class OptimizationTest(unittest.TestCase):
    def setUp(self):
        self.left = 2
        self.right = 7
        index_array = pd.RangeIndex(self.left, self.right, name="ToF_quantized")
        self.data1 = pd.Series([1e-13, 4e-14, 2e-15, 4e-18, 1e-20], index_array, name="Intensity")
        self.data2 = pd.Series([-1e-13, 4e-14, 2e-15, 4e-18, -1e-13], index_array, name="Intensity")

    def test_two_pointer_discrete_optimize_max_is_sum(self):
        def target_func(left, right):
            return self.data1.loc[left : right + 1].sum()

        left, right, optimum = two_pointer_discrete_optimize(target_func, self.left, self.right)
        self.assertEqual(left, self.left)
        self.assertEqual(right, self.right)
        self.assertEqual(optimum, self.data1.sum())

    def test_two_pointer_discrete_optimize_min_is_same_pointer(self):
        def target_func(left, right):
            return self.data1.loc[left : right + 1].sum()

        left, right, optimum = two_pointer_discrete_optimize(target_func, self.left, self.right, "min")
        self.assertEqual(left, right)
        self.assertEqual(optimum, 0.0)

    def test_two_pointer_discrete_optimize_max_slice_around_center(self):
        def target_func(left, right):
            return self.data2.iloc[left : right + 1].sum()

        left, right, optimum = two_pointer_discrete_optimize(target_func, 0, len(self.data2) - 1, "max")
        self.assertEqual(left, 1)
        self.assertEqual(right, 3)
        self.assertEqual(optimum, self.data2.iloc[1:-1].sum())

    def test_two_pointer_brute_force_max(self):
        dummy_data = np.zeros((20, 20))
        max_x, max_y = 5, 7
        max_value = 1.0
        dummy_data[max_x, max_y] = max_value

        def target_func(x, y):
            return dummy_data[x, y]

        left, right, optimum = two_pointer_brute_force_optimize(target_func, 0, 19, "max")
        self.assertEqual(left, max_x)
        self.assertEqual(right, max_y)
        self.assertEqual(optimum, max_value)

    def test_two_pointer_brute_force_min(self):
        dummy_data = np.zeros((20, 20))
        min_x, min_y = 5, 7
        min_value = -1.0
        dummy_data[min_x, min_y] = min_value

        def target_func(x, y):
            return dummy_data[x, y]

        left, right, optimum = two_pointer_brute_force_optimize(target_func, 0, 19, "min")
        self.assertEqual(left, min_x)
        self.assertEqual(right, min_y)
        self.assertEqual(optimum, min_value)


if __name__ == "__main__":
    unittest.main()
