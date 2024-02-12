import unittest
from protein2genome.__main__ import protein2genome
from protein2genome.converter import Converter


class TestSomeFunction(unittest.TestCase):
    def test_map_protein_to_genomic(self):
        converter=Converter('protein2genome/resources/gencode_CDS_DataFrame.parquet')
        
        test_cases=[
            {
                'Transcript ID': 'ENST00000575354',
                'Genome Build': 'GRCh37',
                'Start Position': 199,
                'End Position': 269,
            },
            {
                'Transcript ID': 'ENST00000681038',
                'Genome Build': 'GRCh38',
                'Start Position': 199,
                'End Position': 269
            },
            {
                'Transcript ID': 'ENST00000471181',
                'Genome Build': 'GRCh37',
                'Start Position': 24,
                'End Position': 658
            },
            {
                'Transcript ID': 'ENST00000471181',
                'Genome Build': 'GRCh37',
                'Start Position': 345,
                'End Position': 507
            },
            {
                'Transcript ID': 'ENST00000357654',
                'Genome Build': 'GRCh38',
                'Start Position': 24,
                'End Position': 658
            },
            {
                'Transcript ID': 'ENST00000357654',
                'Genome Build': 'GRCh38',
                'Start Position': 345,
                'End Position': 507
            },
            {
                'Transcript ID': 'ENST00000646891',
                'Genome Build': 'GRCh38',
                'Start Position': 600,
                'End Position': 745
            }
        ]
        
        expected_outputs=[
            ("chr19",42791709,42792003),
            ("chr19",42272378,42272590),
            ("chr17",41245574,41276044),
            ("chr17",41246027,41246515),
            ("chr17",43093557,43124027),
            ("chr17",43094010,43094498),
            ("chr7",140734663,140753337)
        ]
        
        for case_id in range(len(test_cases)):
            test_case=test_cases[case_id]
            test_result=converter.map_protein_to_genomic(test_case['Transcript ID'], test_case['Start Position'], test_case['End Position'], test_case['Genome Build'])
            
            expected_output=expected_outputs[case_id]
            self.assertEqual(test_result,expected_output)
        
if __name__ == '__main__':
    unittest.main()