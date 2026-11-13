import io
import unittest
from unittest.mock import patch
from primer_helper import core

class FastaVariants(unittest.TestCase):
    def test_multisequence_alignment_and_grouped_variants(self):
        records=core.parse_fasta('>ref\nAAC--GTT\n>alt\nAATCCGTT\n>third\nAAC--GTT\n',aligned=True)
        found=core.variants(records[0][1],records[1][1])
        self.assertEqual([(v.kind,v.ref,v.alt,v.ref_start) for v in found],[('SNP','C','T',2),('insertion in comparison','','CC',3)])
    def test_reject_nonmatching_alignment_and_ignore_unknown(self):
        with self.assertRaises(core.InputError): core.parse_fasta('>a\nAAG\n>b\nAA\n',aligned=True)
        self.assertEqual(core.variants('AANC','AATC'),[])
    def test_binding_check_includes_insertions(self):
        rec=[('ref','AC--GTAA'),('alt','ACTTGTAA')]
        checked=core.binding_check(rec,1,4,'CGT')
        self.assertTrue(checked[0]['match']);self.assertFalse(checked[1]['match'])
    def test_reverse_complement(self): self.assertEqual(core.revcomp('ACGTTA'),'TAACGT')

class Design(unittest.TestCase):
    def setUp(self):
        self.ref='A'*70;self.alt='A'*30+'G'+'A'*39
        self.records=[('reference',self.ref),('comparison',self.alt)]
        self.v=core.variants(self.ref,self.alt)[0]
    def fake_output(self,left_start=5):
        return ('PRIMER_PAIR_NUM_RETURNED=1\n'
                f'PRIMER_LEFT_0={left_start},20\nPRIMER_RIGHT_0=55,20\n'
                'PRIMER_LEFT_0_SEQUENCE=AAAAAAAAAAAAAAAAAAAA\n'
                'PRIMER_RIGHT_0_SEQUENCE=TTTTTTTTTTTTTTTTTTTT\n'
                'PRIMER_LEFT_0_TM=60.0\nPRIMER_RIGHT_0_TM=60.0\n=\n')
    @patch('primer_helper.core.run_tool')
    def test_design_uses_primer3_right_coordinate_and_detects_binding(self,run):
        run.return_value=self.fake_output()
        p=core.design(self.records,self.v,40,65,3)[0]
        self.assertEqual(p['sizes'],{'reference':51,'comparison':51})
        self.assertEqual(p['ref_left_1based'],6);self.assertEqual(p['ref_right_1based'],56)
        self.assertTrue(all(x['match'] for x in p['reverse_match']))
        self.assertIn('SEQUENCE_TARGET=30,1',run.call_args.kwargs['input_text'])
    @patch('primer_helper.core.run_tool')
    def test_snp_kasp_places_allele_at_three_prime_end(self,run):
        run.return_value=self.fake_output(11)
        p=list(core.kasp_candidates(self.records,self.v,20,40,65,2))[0]
        self.assertEqual(p['allele_1'],'A'*20);self.assertEqual(p['allele_2'],'A'*19+'G')
        request=run.call_args.kwargs['input_text']
        self.assertIn('SEQUENCE_PRIMER='+('A'*20),request)
        self.assertIn('SEQUENCE_FORCE_LEFT_END=30',request)
        self.assertNotIn('SEQUENCE_TARGET=',request)
    def test_kasp_rejects_other_binding_variant(self):
        alt='A'*20+'C'+'A'*9+'G'+'A'*39
        records=[('reference',self.ref),('comparison',alt)]
        v=core.variants(self.ref,alt)[1]
        with self.assertRaises(core.InputError): list(core.kasp_candidates(records,v,20,40,65,2))
    @patch('primer_helper.core.run_tool')
    def test_reject_product_that_cannot_cover_whole_insertion(self,run):
        ref='A'*30+'---'+'A'*40;alt='A'*30+'CCC'+'A'*40
        v=core.variants(ref,alt)[0]
        run.return_value=self.fake_output(31)
        with self.assertRaises(core.InputError): core.design([('ref',ref),('alt',alt)],v,40,65,1)
    def test_bad_product_range(self):
        with self.assertRaises(core.InputError): core.primer3_input(self.ref,self.v,100,50,5)

class Bench(unittest.TestCase):
    def test_reaction_and_overage(self):
        n,parts=core.reaction_mix(25,2,10,10,.2,2,96,10)
        self.assertEqual(n,106)
        self.assertEqual([v for _,v in parts],[12.5,.5,.5,9.5,2])
    def test_refuses_negative_water(self):
        with self.assertRaises(core.InputError):core.reaction_mix(10,2,1,1,5,1,96,10)
    def test_oligo_and_dilution(self): self.assertEqual(core.oligo_resuspension(9,100,100,10,100),(90,10,90))
    def test_plate_is_column_major(self):
        s='\n'.join(f'S{i}' for i in range(1,11));out=core.plate_csv(s)
        self.assertIn('A1,S1',out);self.assertIn('H1,S8',out);self.assertIn('A2,S9',out)
        with self.assertRaises(core.InputError):core.plate_csv('same\nsame')

if __name__=='__main__':unittest.main()
