from domain_mapper import Transcript, parse_gff

"""
Setup and teardown create a clean and isolated environment for EACH test function or method, ensuring that the state of one test does not interfere with that of the another.
This is why the argument is just a generic variable 'function' - each of the tests below will be fed in as this 'funciton' object
"""
def setup_function(function):
    print("Running setup: %s" % {function.__name__})
    function.gff_in = 'test/mini-hs38.gff3'
    function.build = 'GRCh38'

def teardown_function(function):
    print("Running Teardown: %s" % {function.__name__})
    del function.gff_in


#-------------------------------------------------

def test_gff_parse():
    df = parse_gff(test_gff_parse.gff_in)
    assert df.shape[0] == 20

def test_map_domain():
    #bypass any pandas functions and test the coordinate conversion directly by parsing the mini gff file
    exp_start = 42272378
    exp_end = 42272590
    exp_len = 213
    
    t = Transcript('CIC','ENST00000681038','GRCh38','\"Highmobility groupboxdomain:199-269\"')
    t.strand = "+"
    t.chrom = "19"
    with open(test_map_domain.gff_in, 'r') as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            ll = line.rstrip().split("\t")
            if ll[2] == 'CDS':
                t.cds_coords.append(( int(ll[3]),int(ll[4]) ))
    t.map_domains_to_genome()

    for dm in t.domains.values():
        assert dm.genomic_start == exp_start
        assert dm.genomic_end == exp_end
        assert dm.genomic_len == exp_len
        
            

    
