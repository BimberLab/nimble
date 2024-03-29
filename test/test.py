import unittest
import pandas as pd
import numpy as np

import nimble
from nimble.reporting import _min_count, _min_pct


class TestFiltering(unittest.TestCase):
    def test_min_count_100(self):
        data = {
            0: [25, 25, 25, 25],
            1: ["C", "B", "A", "A"],
            2: ["B", "A", "B", np.nan],
            3: ["D", "C", np.nan, np.nan],
            4: ["A", np.nan, np.nan, np.nan],
        }

        data = pd.DataFrame(data)
        data = _min_count(data, 100)

        expected = {
            0: [25, 25, 25, 25],
            1: [np.nan, np.nan, "A", "A"],
            2: [np.nan, "A", np.nan, np.nan],
            3: [np.nan, np.nan, np.nan, np.nan],
            4: ["A", np.nan, np.nan, np.nan],
        }
        expected = pd.DataFrame(expected)
        pd.testing.assert_frame_equal(data, expected)

    def test_min_count_75(self):
        data = {
            0: [25, 25, 25, 25],
            1: ["C", "B", "A", "A"],
            2: ["B", "A", "B", np.nan],
            3: ["D", "C", np.nan, np.nan],
            4: ["A", np.nan, np.nan, np.nan],
        }

        data = pd.DataFrame(data)
        data = _min_count(data, 75)

        expected = {
            0: [25, 25, 25, 25],
            1: [np.nan, "B", "A", "A"],
            2: ["B", "A", "B", np.nan],
            3: [np.nan, np.nan, np.nan, np.nan],
            4: ["A", np.nan, np.nan, np.nan],
        }
        expected = pd.DataFrame(expected)
        pd.testing.assert_frame_equal(data, expected)

    def test_min_count_50(self):
        data = {
            0: [25, 25, 25, 25],
            1: ["C", "B", "A", "A"],
            2: ["B", "A", "B", np.nan],
            3: ["D", "C", np.nan, np.nan],
            4: ["A", np.nan, np.nan, np.nan],
        }

        data = pd.DataFrame(data)
        data = _min_count(data, 50)

        expected = {
            0: [25, 25, 25, 25],
            1: ["C", "B", "A", "A"],
            2: ["B", "A", "B", np.nan],
            3: [np.nan, "C", np.nan, np.nan],
            4: ["A", np.nan, np.nan, np.nan],
        }
        expected = pd.DataFrame(expected)
        pd.testing.assert_frame_equal(data, expected)

    def test_min_pct_100(self):
        data = {
            0: [25, 25, 25, 25],
            1: ["C", "B", "A", "A"],
            2: ["B", "A", "B", np.nan],
            3: ["D", "C", np.nan, np.nan],
            4: ["A", np.nan, np.nan, np.nan],
        }

        data = pd.DataFrame(data)
        data = _min_pct(data, 1)

        expected = {
            0: [25, 25, 25, 25],
            1: [np.nan, np.nan, "A", "A"],
            2: [np.nan, "A", np.nan, np.nan],
            3: [np.nan, np.nan, np.nan, np.nan],
            4: ["A", np.nan, np.nan, np.nan],
        }
        expected = pd.DataFrame(expected)
        pd.testing.assert_frame_equal(data, expected)

    def test_min_pct_75(self):
        data = {
            0: [25, 25, 25, 25],
            1: ["C", "B", "A", "A"],
            2: ["B", "A", "B", np.nan],
            3: ["D", "C", np.nan, np.nan],
            4: ["A", np.nan, np.nan, np.nan],
        }

        data = pd.DataFrame(data)
        data = _min_pct(data, 0.75)

        expected = {
            0: [25, 25, 25, 25],
            1: [np.nan, "B", "A", "A"],
            2: ["B", "A", "B", np.nan],
            3: [np.nan, np.nan, np.nan, np.nan],
            4: ["A", np.nan, np.nan, np.nan],
        }
        expected = pd.DataFrame(expected)
        pd.testing.assert_frame_equal(data, expected)

    def test_min_pct_50(self):
        data = {
            0: [25, 25, 25, 25],
            1: ["C", "B", "A", "A"],
            2: ["B", "A", "B", np.nan],
            3: ["D", "C", np.nan, np.nan],
            4: ["A", np.nan, np.nan, np.nan],
        }

        data = pd.DataFrame(data)
        data = _min_pct(data, 0.5)

        expected = {
            0: [25, 25, 25, 25],
            1: ["C", "B", "A", "A"],
            2: ["B", "A", "B", np.nan],
            3: [np.nan, "C", np.nan, np.nan],
            4: ["A", np.nan, np.nan, np.nan],
        }
        expected = pd.DataFrame(expected)
        pd.testing.assert_frame_equal(data, expected)


if __name__ == "__main__":
    unittest.main()
