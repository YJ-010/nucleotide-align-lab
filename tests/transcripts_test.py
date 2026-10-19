import pathlib
import sys
import unittest
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]/'scripts'))
from build_transcripts import assemble, rc


class TranscriptTests(unittest.TestCase):
    def setUp(self):
        self.dna = 'AAAATGCCCCGCTTAATTTT'
        self.record = {'chrom': 'chr1', 'strand': '+', 'cds': 'ATGGCTTAA',
                       'cdsSegments': [[4, 6, 0], [11, 16, 0]]}
        self.exons = [('chr1', 1, 6, '+'), ('chr1', 11, 20, '+')]
        self.tx = {'start': 1, 'end': 20}

    def test_positive_and_negative_introns_are_removed(self):
        for strand in ['+', '-']:
            record = dict(self.record, strand=strand)
            exons, dna = self.exons, self.dna
            if strand == '-':
                dna = rc(dna)
                exons = [('chr1', 21-e, 21-s, '-') for _, s, e, _ in exons]
                record['cdsSegments'] = sorted([[21-e, 21-s, p] for s, e, p in self.record['cdsSegments']])
            result = assemble(record, {'chr1': dna}, exons, [], self.tx)
            self.assertEqual(result['status'], 'available')
            self.assertEqual(result['prefix'], 'AAA')
            self.assertEqual(result['suffix'], 'TTTT')
            self.assertEqual(result['prefix']+record['cds']+result['suffix'], 'AAAATGGCTTAATTTT')
            self.assertEqual((result['cdsStart'], result['cdsEnd']), (4, 12))

    def test_cds_and_utr_features_reconstruct_white_style_transcript(self):
        result = assemble(self.record, {'chr1': self.dna}, [], [('chr1', 1, 3, '+'), ('chr1', 17, 20, '+')], self.tx)
        self.assertEqual(result['method'], 'cds_and_utr')
        self.assertEqual(result['exons'], [[1, 6], [11, 20]])
        self.assertEqual(result['length'], 16)

    def test_unannotated_utr_is_not_invented_and_bad_models_are_refused(self):
        self.assertEqual(assemble(self.record, {}, [], [], self.tx)['reason'], 'no_exon_or_utr')
        bad = dict(self.record, cds='ATGGCCTAA')
        self.assertEqual(assemble(bad, {'chr1': self.dna}, self.exons, [], self.tx)['reason'], 'cds_sequence_mismatch')
        overlap = self.exons+[('chr1', 5, 12, '+')]
        self.assertEqual(assemble(self.record, {'chr1': self.dna}, overlap, [], self.tx)['reason'], 'overlapping_exons')

    def test_phase_bases_remain_in_transcript_without_trimming_exons(self):
        record = dict(self.record, cds='TGGCTTAA', cdsSegments=[[4, 6, 1], [11, 16, 0]])
        result = assemble(record, {'chr1': self.dna}, self.exons, [], self.tx)
        self.assertEqual(result['prefix'], 'AAAA')
        self.assertIn('partial_cds_phase', result['warnings'])
        self.assertEqual(result['length'], 16)


if __name__ == '__main__':
    unittest.main()
