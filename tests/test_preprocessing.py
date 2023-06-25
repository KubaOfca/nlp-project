'''
This code aims to test the preprocessing function.
'''
import pandas as pd
from title_generator.preprocesing import pre_processing


def test_preprocessing():
    # define data to check
    test_df = pd.DataFrame(
        {
            "Column1": [
                "<p>Hello, <b>https://github.com/oasdn676/asda.git  World!</b></p>",
                "Some <b>text ( bla bla)</b>",
            ],
            "Column2": [
                "<h1>Title!!!! 1983 aby żeby ?! > ? ! zeby [ { ( ) } ].</h1>",
                "Another !<!i>paragraph 1832r {bla bla (bla!) [bla bla]} </i>",
            ],
        }
    )

    df = pre_processing(test_df)

    # define expected result 
    test_df_after_pre_processing = pd.DataFrame(
        {
            "Column1": [
                ["hello", "world"],
                ["some", "text"],
            ],
            "Column2": [
                ["title"],
                ["another", "paragraph"],
            ],
        }
    )

    # check if the actual result is the same as expected result
    assert df.equals(test_df_after_pre_processing)
