import unittest
import pandas as pd
import numpy as np

import nimble
from nimble.utils import per_umi_thresholding, umi_intersection, intersect_lists

class TestDataProcessing(unittest.TestCase):

    def test_per_umi_thresholding_basic(self):
        """
        Basic test to ensure per_umi_thresholding correctly filters features based on threshold.
        """
        data = pd.DataFrame({
            'cb': ['cell1', 'cell1'],
            'umi': ['UMI1', 'UMI1'],
            'features': ['A,B', 'A,C'],
            'nimble_score': [10, 20]
        })
        threshold = 0.2
        expected_output = pd.DataFrame({
            'cb': ['cell1', 'cell1'],
            'umi': ['UMI1', 'UMI1'],
            'features': ['A,B', 'A,C'],
            'nimble_score': [10, 20],
            'filtered_features': ['A', 'A,C']
        })

        result = per_umi_thresholding(data, threshold)
        pd.testing.assert_frame_equal(result.reset_index(drop=True), expected_output)

    def test_per_umi_thresholding_edge_case_all_below_threshold(self):
        """
        Test where all features are below the threshold and should be filtered out.
        """
        data = pd.DataFrame({
            'cb': ['cell1'],
            'umi': ['UMI1'],
            'features': ['A,B,C'],
            'nimble_score': [3]  # Each feature gets a score of 1
        })
        threshold = 0.4  # Each feature has a ratio of 1/3, which is ~0.33
        result = per_umi_thresholding(data, threshold)
        # Since all features are below the threshold, the result should be empty
        self.assertTrue(result.empty)

    def test_per_umi_thresholding_single_feature(self):
        """
        Test with UMIs that have only one feature; should not be filtered out regardless of threshold.
        """
        data = pd.DataFrame({
            'cb': ['cell1'],
            'umi': ['UMI1'],
            'features': ['A'],
            'nimble_score': [10]
        })
        threshold = 0.9
        result = per_umi_thresholding(data, threshold)
        # The single feature should remain
        self.assertTrue(result['filtered_features'].iloc[0] == 'A')

    def test_umi_intersection_basic(self):
        """
        Basic test to ensure umi_intersection correctly computes the intersection of features.
        """
        data = pd.DataFrame({
            'cb': ['cell1', 'cell1', 'cell1'],
            'umi': ['UMI1', 'UMI1', 'UMI1'],
            'filtered_features': ['A,B', 'A,C', 'A,D']
        })
        expected_output = pd.DataFrame({
            'cb': ['cell1'],
            'umi': ['UMI1'],
            'filtered_features': [['A']]
        })

        result = umi_intersection(data)
        self.assertEqual(result['filtered_features'].iloc[0], expected_output['filtered_features'].iloc[0])

    def test_umi_intersection_no_common_features(self):
        """
        Test where there are no common features among the read-mates; intersection should be empty.
        """
        data = pd.DataFrame({
            'cb': ['cell1', 'cell1'],
            'umi': ['UMI1', 'UMI1'],
            'filtered_features': ['A,B', 'C,D']
        })
        result = umi_intersection(data)
        # The intersection is empty
        self.assertEqual(result['filtered_features'].iloc[0], [])

    def test_intersect_lists_empty_input(self):
        """
        Test intersect_lists with empty input.
        """
        result = intersect_lists([])
        self.assertEqual(result, [])

    def test_intersect_lists_single_list(self):
        """
        Test intersect_lists with a single list.
        """
        result = intersect_lists([['A', 'B', 'C']])
        self.assertEqual(sorted(result), ['A', 'B', 'C'])

    def test_data_integrity_multiple_cells(self):
        """
        Pipeline test with multiple cells and UMIs to ensure data integrity.
        """
        initial_data = pd.DataFrame({
            'r1_CB': ['cell1', 'cell1', 'cell2', 'cell2', 'cell3'],
            'r1_UB': ['UMI1', 'UMI1', 'UMI2', 'UMI2', 'UMI3'],
            'features': ['A,B', 'A,C', 'D,E', 'D,F', 'G'],
            'nimble_score': [10, 20, 30, 40, 50]
        })

        # Simulate the full pipeline
        initial_data.rename(columns={'r1_CB': 'cb', 'r1_UB': 'umi'}, inplace=True)
        df = initial_data[['features', 'umi', 'cb', 'nimble_score']].copy()
        df['features'] = df['features'].apply(lambda x: ','.join(sorted(x.split(','))))
        df = df.groupby(['cb', 'umi', 'features'])['nimble_score'].sum().reset_index()
        df = per_umi_thresholding(df, threshold=0.2)
        df_grouped = umi_intersection(df)
        df_grouped = df_grouped[~df_grouped['filtered_features'].apply(lambda x: len(x) == 0)]
        df_grouped['filtered_features'] = df_grouped['filtered_features'].apply(lambda x: ','.join(x))
        df_grouped.columns = ['cell_barcode', 'umi', 'feature']
        df_counts = df_grouped.groupby(['cell_barcode', 'feature']).size().reset_index(name='count')
        df_counts = df_counts[['feature', 'count', 'cell_barcode']]

        # Expected output after processing
        expected_counts = pd.DataFrame({
            'feature': ['A', 'D', 'G'],
            'count': [1, 1, 1],
            'cell_barcode': ['cell1', 'cell2', 'cell3']
        })

        pd.testing.assert_frame_equal(df_counts.reset_index(drop=True), expected_counts)

    def test_data_integrity_no_features_remaining(self):
        """
        Test the pipeline when no features remain after thresholding and intersection.
        """
        initial_data = pd.DataFrame({
            'r1_CB': ['cell1'],
            'r1_UB': ['UMI1'],
            'features': ['A,B,C'],
            'nimble_score': [3]
        })

        # Each feature will have a nimble_score contributing less than the threshold
        initial_data.rename(columns={'r1_CB': 'cb', 'r1_UB': 'umi'}, inplace=True)
        df = initial_data[['features', 'umi', 'cb', 'nimble_score']].copy()
        df['features'] = df['features'].apply(lambda x: ','.join(sorted(x.split(','))))
        df = df.groupby(['cb', 'umi', 'features'])['nimble_score'].sum().reset_index()
        df = per_umi_thresholding(df, threshold=0.4)
        df_grouped = umi_intersection(df)
        df_grouped = df_grouped[~df_grouped['filtered_features'].apply(lambda x: len(x) == 0)]
        # Since no features remain, df_grouped should be empty
        self.assertTrue(df_grouped.empty)

    def test_data_integrity_duplicate_reads(self):
        """
        Test the pipeline with duplicate reads to ensure they are merged correctly.
        """
        initial_data = pd.DataFrame({
            'r1_CB': ['cell1', 'cell1'],
            'r1_UB': ['UMI1', 'UMI1'],
            'features': ['A,B', 'A,B'],
            'nimble_score': [10, 10]
        })

        # Simulate the full pipeline
        initial_data.rename(columns={'r1_CB': 'cb', 'r1_UB': 'umi'}, inplace=True)
        df = initial_data[['features', 'umi', 'cb', 'nimble_score']].copy()
        df['features'] = df['features'].apply(lambda x: ','.join(sorted(x.split(','))))
        df = df.groupby(['cb', 'umi', 'features'])['nimble_score'].sum().reset_index()
        df = per_umi_thresholding(df, threshold=0.1)
        df_grouped = umi_intersection(df)
        df_grouped = df_grouped[~df_grouped['filtered_features'].apply(lambda x: len(x) == 0)]
        df_grouped['filtered_features'] = df_grouped['filtered_features'].apply(lambda x: ','.join(x))
        df_grouped.columns = ['cell_barcode', 'umi', 'feature']
        df_counts = df_grouped.groupby(['cell_barcode', 'feature']).size().reset_index(name='count')
        df_counts = df_counts[['feature', 'count', 'cell_barcode']]

        # Expected output after processing
        expected_counts = pd.DataFrame({
            'feature': ['A,B'],
            'count': [1],
            'cell_barcode': ['cell1']
        })

        pd.testing.assert_frame_equal(df_counts.reset_index(drop=True), expected_counts)

    def test_per_umi_thresholding_high_threshold(self):
        """
        Test per_umi_thresholding with a high threshold to see if features are correctly dropped.
        """
        data = pd.DataFrame({
            'cb': ['cell1', 'cell1'],
            'umi': ['UMI1', 'UMI1'],
            'features': ['A,B', 'A,B,C,D'],
            'nimble_score': [100, 100]
        })
        threshold = 0.3
        result = per_umi_thresholding(data, threshold)
        self.assertTrue(set(result['filtered_features'].iloc[0].split(',')) == {'A', 'B'})

    def test_per_umi_thresholding_non_uniform_scores(self):
        """
        Test per_umi_thresholding with non-uniform nimble_scores across features.
        """
        data = pd.DataFrame({
            'cb': ['cell1', 'cell1'],
            'umi': ['UMI1', 'UMI1'],
            'features': ['A', 'B,C'],
            'nimble_score': [80, 20]
        })
        threshold = 0.25
        result = per_umi_thresholding(data, threshold)
        expected_features = {'A'}
        self.assertTrue(set(result['filtered_features'].str.cat(sep=',').split(',')) == expected_features)

    def test_umi_intersection_single_read(self):
        """
        Test umi_intersection with only one read for a UMI.
        """
        data = pd.DataFrame({
            'cb': ['cell1'],
            'umi': ['UMI1'],
            'filtered_features': ['A,B,C']
        })
        result = umi_intersection(data)
        # The intersection should be the features of the single read
        self.assertEqual(set(result['filtered_features'].iloc[0]), {'A', 'B', 'C'})



class TestDataProcessing(unittest.TestCase):

    # Existing tests...

    def test_per_umi_thresholding_complex_scores(self):
        """
        Test per_umi_thresholding with complex score distributions and multiple features per UMI.
        """
        data = pd.DataFrame({
            'cb': ['cell1', 'cell1', 'cell1', 'cell1'],
            'umi': ['UMI1', 'UMI1', 'UMI1', 'UMI1'],
            'features': ['A,B', 'A,C', 'B,C,D', 'D,E'],
            'nimble_score': [10, 15, 5, 20]
        })
        threshold = 0.2
        result = per_umi_thresholding(data, threshold)
        # Expected features after thresholding: 'A', 'E', 'D'
        expected_features = {'A', 'E', 'D'}
        filtered_features = set(result['filtered_features'].str.cat(sep=',').split(','))

        self.assertEqual(filtered_features, expected_features)

    def test_umi_intersection_complex(self):
        """
        Test umi_intersection with complex feature sets and multiple read-mates.
        """
        data = pd.DataFrame({
            'cb': ['cell1', 'cell1', 'cell1', 'cell1'],
            'umi': ['UMI1', 'UMI1', 'UMI1', 'UMI1'],
            'filtered_features': ['A,B,C', 'A,C', 'B,C,D', 'C,D,E']
        })
        result = umi_intersection(data)
        # Expected intersection is 'C'
        self.assertEqual(result['filtered_features'].iloc[0], ['C'])

    def test_pipeline_integration_complex(self):
        """
        Integration test simulating the entire pipeline with complex data, multiple cells and UMIs,
        varying scores, and complex patterns of ambiguity.
        """
        initial_data = pd.DataFrame({
            'r1_CB': ['cell1', 'cell1', 'cell1', 'cell2', 'cell2', 'cell3', 'cell3', 'cell3'],
            'r1_UB': ['UMI1', 'UMI1', 'UMI2', 'UMI3', 'UMI3', 'UMI4', 'UMI5', 'UMI5'],
            'features': ['A,B', 'A,C', 'B,D', 'E,F', 'F,G', 'H,I', 'I,J', 'H,J'],
            'nimble_score': [10, 20, 15, 5, 35, 25, 15, 10]
        })

        # Simulate the full pipeline
        initial_data.rename(columns={'r1_CB': 'cb', 'r1_UB': 'umi'}, inplace=True)
        df = initial_data[['features', 'umi', 'cb', 'nimble_score']].copy()

        # Ensure features are sorted
        df['features'] = df['features'].apply(lambda x: ','.join(sorted(x.split(','))))

        # Group by cb, umi, features and sum nimble_score
        df = df.groupby(['cb', 'umi', 'features'])['nimble_score'].sum().reset_index()

        # Apply per-UMI thresholding
        df = per_umi_thresholding(df, threshold=0.2)

        # Apply UMI intersection
        df_grouped = umi_intersection(df)
        df_grouped = df_grouped[df_grouped['filtered_features'].apply(lambda x: len(x) > 0)]

        # Convert filtered_features back to strings
        df_grouped['filtered_features'] = df_grouped['filtered_features'].apply(lambda x: ','.join(x))
        df_grouped.columns = ['cell_barcode', 'umi', 'feature']

        # Group by cell_barcode and feature to get counts
        df_counts = df_grouped.groupby(['cell_barcode', 'feature']).size().reset_index(name='count')
        df_counts = df_counts[['feature', 'count', 'cell_barcode']]

        # Expected output after processing
        expected_counts = pd.DataFrame({
            'feature': ['A', 'B,D', 'F', 'H,I', 'J'],
            'count': [1, 1, 1, 1, 1],
            'cell_barcode': ['cell1', 'cell1', 'cell2', 'cell3', 'cell3']
        })

        # Sort the dataframes for consistent comparison
        df_counts = df_counts.sort_values(by=['cell_barcode', 'feature']).reset_index(drop=True)
        expected_counts = expected_counts.sort_values(by=['cell_barcode', 'feature']).reset_index(drop=True)

        pd.testing.assert_frame_equal(df_counts, expected_counts)

    def test_per_umi_thresholding_tie_scores(self):
        """
        Test per_umi_thresholding when features have equal scores exactly at the threshold boundary.
        """
        data = pd.DataFrame({
            'cb': ['cell1'],
            'umi': ['UMI1'],
            'features': ['A,B'],
            'nimble_score': [10]
        })
        threshold = 0.5  # Both features have a ratio of 0.5
        result = per_umi_thresholding(data, threshold)
        # Both features should remain
        self.assertEqual(set(result['filtered_features'].iloc[0].split(',')), {'A', 'B'})

    def test_per_umi_thresholding_zero_scores(self):
        """
        Test per_umi_thresholding when some features have zero nimble_score.
        """
        data = pd.DataFrame({
            'cb': ['cell1', 'cell1'],
            'umi': ['UMI1', 'UMI1'],
            'features': ['A,B', 'C,D'],
            'nimble_score': [0, 20]
        })
        threshold = 0.1
        result = per_umi_thresholding(data, threshold)
        # Features 'C' and 'D' should remain; 'A' and 'B' should be dropped due to zero score
        expected_features = {'C', 'D'}
        filtered_features = set(result['filtered_features'].str.cat(sep=',').split(','))
        self.assertEqual(filtered_features, expected_features)

    def test_umi_intersection_no_reads(self):
        """
        Test umi_intersection when there are no reads for a UMI (empty DataFrame).
        """
        data = pd.DataFrame(columns=['cb', 'umi', 'filtered_features'])
        result = umi_intersection(data)
        # Result should be an empty DataFrame
        self.assertTrue(result.empty)

    def test_pipeline_integration_no_thresholding(self):
        """
        Integration test simulating the entire pipeline with thresholding disabled (threshold=0).
        """
        initial_data = pd.DataFrame({
            'r1_CB': ['cell1', 'cell2', 'cell2', 'cell3'],
            'r1_UB': ['UMI1', 'UMI2', 'UMI2', 'UMI3'],
            'features': ['A,B', 'C,D', 'D,E', 'F,G'],
            'nimble_score': [10, 20, 30, 40]
        })

        # Simulate the full pipeline
        initial_data.rename(columns={'r1_CB': 'cb', 'r1_UB': 'umi'}, inplace=True)
        df = initial_data[['features', 'umi', 'cb', 'nimble_score']].copy()

        # Ensure features are sorted
        df['features'] = df['features'].apply(lambda x: ','.join(sorted(x.split(','))))

        # Group by cb, umi, features and sum nimble_score
        df = df.groupby(['cb', 'umi', 'features'])['nimble_score'].sum().reset_index()

        # Apply per-UMI thresholding with threshold=0 (no thresholding)
        df = per_umi_thresholding(df, threshold=0)

        # Apply UMI intersection
        df_grouped = umi_intersection(df)
        df_grouped = df_grouped[df_grouped['filtered_features'].apply(lambda x: len(x) > 0)]

        # Convert filtered_features back to strings
        df_grouped['filtered_features'] = df_grouped['filtered_features'].apply(lambda x: ','.join(x))
        df_grouped.columns = ['cell_barcode', 'umi', 'feature']

        # Group by cell_barcode and feature to get counts
        df_counts = df_grouped.groupby(['cell_barcode', 'feature']).size().reset_index(name='count')
        df_counts = df_counts[['feature', 'count', 'cell_barcode']]

        # Expected output after processing
        expected_counts = pd.DataFrame({
            'feature': ['A,B', 'D', 'F,G'],
            'count': [1, 1, 1],
            'cell_barcode': ['cell1', 'cell2', 'cell3']
        })

        # Sort the dataframes for consistent comparison
        df_counts = df_counts.sort_values(by=['cell_barcode', 'feature']).reset_index(drop=True)
        expected_counts = expected_counts.sort_values(by=['cell_barcode', 'feature']).reset_index(drop=True)

        pd.testing.assert_frame_equal(df_counts, expected_counts)

    def test_per_umi_thresholding_duplicate_features(self):
        """
        Test per_umi_thresholding with duplicate features in a read-mate to ensure they are handled correctly.
        """
        data = pd.DataFrame({
            'cb': ['cell1'],
            'umi': ['UMI1'],
            'features': ['A,A,B'],
            'nimble_score': [15]
        })
        threshold = 0.2
        result = per_umi_thresholding(data, threshold)
        # Expected features after thresholding: 'A', 'B'
        expected_features = {'A', 'B'}
        filtered_features = set(result['filtered_features'].iloc[0].split(','))
        self.assertEqual(filtered_features, expected_features)

    def test_pipeline_integration_realistic_data(self):
        """
        Integration test with a realistic dataset simulating potential real-world complexity.
        """
        # Simulate a dataset with varying UMIs, cells, features, and scores
        initial_data = pd.DataFrame({
            'r1_CB': ['cell1'] * 5 + ['cell2'] * 4 + ['cell3'] * 3,
            'r1_UB': ['UMI1', 'UMI1', 'UMI2', 'UMI2', 'UMI2', 'UMI3', 'UMI3', 'UMI4', 'UMI4', 'UMI5', 'UMI5', 'UMI5'],
            'features': ['A', 'B', 'A,B', 'B,C', 'C', 'D', 'E', 'F', 'F,G', 'H,I', 'I,J', 'H,J'],
            'nimble_score': [10, 5, 8, 12, 3, 20, 15, 25, 5, 10, 15, 5]
        })

        # Simulate the full pipeline
        initial_data.rename(columns={'r1_CB': 'cb', 'r1_UB': 'umi'}, inplace=True)
        df = initial_data[['features', 'umi', 'cb', 'nimble_score']].copy()

        # Ensure features are sorted
        df['features'] = df['features'].apply(lambda x: ','.join(sorted(x.split(','))))

        # Group by cb, umi, features and sum nimble_score
        df = df.groupby(['cb', 'umi', 'features'])['nimble_score'].sum().reset_index()

        # Apply per-UMI thresholding
        df = per_umi_thresholding(df, threshold=0.15)

        # Apply UMI intersection
        df_grouped = umi_intersection(df)
        df_grouped = df_grouped[df_grouped['filtered_features'].apply(lambda x: len(x) > 0)]

        # Convert filtered_features back to strings
        df_grouped['filtered_features'] = df_grouped['filtered_features'].apply(lambda x: ','.join(x))
        df_grouped.columns = ['cell_barcode', 'umi', 'feature']

        # Group by cell_barcode and feature to get counts
        df_counts = df_grouped.groupby(['cell_barcode', 'feature']).size().reset_index(name='count')
        df_counts = df_counts[['feature', 'count', 'cell_barcode']]

        # Expected output after processing
        expected_counts = pd.DataFrame({
            'feature': ['F'],
            'count': [1],
            'cell_barcode': ['cell2']
        })

        # Sort the dataframes for consistent comparison
        df_counts = df_counts.sort_values(by=['cell_barcode', 'feature']).reset_index(drop=True)
        expected_counts = expected_counts.sort_values(by=['cell_barcode', 'feature']).reset_index(drop=True)

        pd.testing.assert_frame_equal(df_counts, expected_counts)

    def test_per_umi_thresholding_large_scores(self):
        """
        Test per_umi_thresholding with very large nimble_scores to check for any overflow issues.
        """
        data = pd.DataFrame({
            'cb': ['cell1', 'cell1'],
            'umi': ['UMI1', 'UMI1'],
            'features': ['A,B,C', 'C,D,E'],
            'nimble_score': [1e12, 1e12]
        })
        threshold = 0.2
        result = per_umi_thresholding(data, threshold)
        expected_features = {'C'}
        filtered_features = set(result['filtered_features'].str.cat(sep=',').split(','))
        self.assertEqual(filtered_features, expected_features)

    def test_per_umi_thresholding_decimal_scores(self):
        """
        Test per_umi_thresholding with decimal (float) nimble_scores.
        """
        data = pd.DataFrame({
            'cb': ['cell1', 'cell1'],
            'umi': ['UMI1', 'UMI1'],
            'features': ['A,B', 'A,C'],
            'nimble_score': [0.6, 0.4]
        })
        threshold = 0.5
        result = per_umi_thresholding(data, threshold)
        expected_features = {'A'}
        filtered_features = set(result['filtered_features'].str.cat(sep=',').split(','))
        self.assertEqual(filtered_features, expected_features)


if __name__ == "__main__":
    unittest.main()
