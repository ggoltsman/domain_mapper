from domain_mapper import parse_gff

"""
Setup and teardown create a clean and isolated environment for EACH test function or method, ensuring that the state of one test does not interfere with that of the another.
This is why the argument is just a generic variable 'function' - each of the tests below will be fed in as this 'funciton' object
"""
def setup_function(function):
    print("Running setup: %s" % {function.__name__})
    function.gff_in = 'test/mini-gff.gff3'
    function.build = 'GRCh38'

def teardown_function(function):
    print("Running Teardown: %s" % {function.__name__})
    del function.gff_in
    del function.bad_gff_in


#-------------------------------------------------

def test_gff_parse():
    df = parse_gff(test_gff_parse.gff_in)
    assert df.shape[0] == 20

def test_parse_cds():
    bad_gff_in = 'test/bad-gff-2.gff3'
    

    
