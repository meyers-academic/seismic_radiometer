import unittest
def get_suite():
    import seispy.tests
    loader=unittest.TestLoader()
    suite=loader.loadTestsFromModule(seispy.tests)
    return suite
