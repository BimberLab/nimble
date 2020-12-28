import unittest
import pandas as pd
import numpy as np

from immunogenotyper.reporting import min_count, min_pct


class TestFiltering(unittest.TestCase):
  def test_min_count_4(self):
    data = {'1': ["C", "B", "A", "A"],
         '2': ["B", "A", "B", np.nan],
         '3': ["D", "C", np.nan, np.nan],
         '4': ["A", np.nan, np.nan, np.nan]
        }

    data = pd.DataFrame(data)
    data = min_count(data, 4)

    expected = {'1': [np.nan, np.nan, "A", "A"],
         '2': [np.nan, "A", np.nan, np.nan],
         '3': [np.nan, np.nan, np.nan, np.nan],
         '4': ["A", np.nan, np.nan, np.nan]
        }
    expected = pd.DataFrame(expected)
    pd.testing.assert_frame_equal(data, expected)


  def test_min_count_3(self):
    data = {'1': ["C", "B", "A", "A"],
         '2': ["B", "A", "B", np.nan],
         '3': ["D", "C", np.nan, np.nan],
         '4': ["A", np.nan, np.nan, np.nan]
        }

    data = pd.DataFrame(data)
    data = min_count(data, 3)

    expected = {'1': [np.nan, "B", "A", "A"],
         '2': ["B", "A", "B", np.nan],
         '3': [np.nan, np.nan, np.nan, np.nan],
         '4': ["A", np.nan, np.nan, np.nan]
        }
    expected = pd.DataFrame(expected)
    pd.testing.assert_frame_equal(data, expected)


  def test_min_count_2(self):
    data = {'1': ["C", "B", "A", "A"],
         '2': ["B", "A", "B", np.nan],
         '3': ["D", "C", np.nan, np.nan],
         '4': ["A", np.nan, np.nan, np.nan]
        }

    data = pd.DataFrame(data)
    data = min_count(data, 2)

    expected = {'1': ["C", "B", "A", "A"],
         '2': ["B", "A", "B", np.nan],
         '3': [np.nan, "C", np.nan, np.nan],
         '4': ["A", np.nan, np.nan, np.nan]
        }
    expected = pd.DataFrame(expected)
    pd.testing.assert_frame_equal(data, expected)


  def test_min_pct_100(self):
    data = {'1': ["C", "B", "A", "A"],
         '2': ["B", "A", "B", np.nan],
         '3': ["D", "C", np.nan, np.nan],
         '4': ["A", np.nan, np.nan, np.nan]
        }

    data = pd.DataFrame(data)
    data = min_pct(data, 1)

    expected = {'1': [np.nan, np.nan, "A", "A"],
         '2': [np.nan, "A", np.nan, np.nan],
         '3': [np.nan, np.nan, np.nan, np.nan],
         '4': ["A", np.nan, np.nan, np.nan]
        }
    expected = pd.DataFrame(expected)
    pd.testing.assert_frame_equal(data, expected)


  def test_min_pct_75(self):
    data = {'1': ["C", "B", "A", "A"],
         '2': ["B", "A", "B", np.nan],
         '3': ["D", "C", np.nan, np.nan],
         '4': ["A", np.nan, np.nan, np.nan]
        }

    data = pd.DataFrame(data)
    data = min_pct(data, .75)

    expected = {'1': [np.nan, "B", "A", "A"],
         '2': ["B", "A", "B", np.nan],
         '3': [np.nan, np.nan, np.nan, np.nan],
         '4': ["A", np.nan, np.nan, np.nan]
        }
    expected = pd.DataFrame(expected)
    pd.testing.assert_frame_equal(data, expected)


  def test_min_pct_50(self):
    data = {'1': ["C", "B", "A", "A"],
         '2': ["B", "A", "B", np.nan],
         '3': ["D", "C", np.nan, np.nan],
         '4': ["A", np.nan, np.nan, np.nan]
        }

    data = pd.DataFrame(data)
    data = min_pct(data, .5)

    expected = {'1': ["C", "B", "A", "A"],
         '2': ["B", "A", "B", np.nan],
         '3': [np.nan, "C", np.nan, np.nan],
         '4': ["A", np.nan, np.nan, np.nan]
        }
    expected = pd.DataFrame(expected)
    pd.testing.assert_frame_equal(data, expected)

if __name__ == "__main__":
  unittest.main()