import uncalled4 as unc4

#placeholder for more comprehensive unit testing
#see uncalled4/example/ for integration tests

def test_pore_models():
    dna_r9 = unc4.PoreModel("dna_r9.4.1_400bps_6mer")
    assert(dna_r9.k == 6)

    dna_r10 = unc4.PoreModel("dna_r10.4.1_400bps_9mer")
    assert(dna_r10.k == 9)

    rna_r9 = unc4.PoreModel("rna_r9.4.1_70bps_5mer")
    assert(rna_r9.k == 5)
