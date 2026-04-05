"""
Report predispositions subpackage unit test.
"""

import os

from tsoppy.reportPredispositions.predispositions import load_data_from_cancer_susceptibility_genes_table, report_predispositions


def test_load_data_from_cancer_susceptibility_genes_table(tmp_path):

    # create temporary file
    directory = tmp_path / 'tsoppy_test'
    directory.mkdir()
    test_file = directory / 'csg.csv'

    # write content to the temp file
    csg_table_content = "Gene\tActionability\tAge\nTP53\tHA-CSG\tAge<30\n"
    test_file.write_text(csg_table_content)

    # create result
    result = dict()
    result['TP53'] = dict()
    result['TP53']['Actionability'] = 'HA-CSG'
    result['TP53']['Age'] = 'Age<30'

    # pass the file to the function
    assert load_data_from_cancer_susceptibility_genes_table(
        test_file) == result


# TODO


def test_report_predispositions():
    """
    Unit test for the report predisposition function of the subpackage.
    """
    assert True
