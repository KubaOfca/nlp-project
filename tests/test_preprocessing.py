import pandas as pd
from title_generator.preprocesing import pre_processing


def test_preprocessing():
    test_df = pd.DataFrame(
        {
            "Column1": [
                "<p>Hello, <b>https://github.com/oasdn676/asda.git  World!</b></p>",
                "Some <b>text ( bla bla)</b>",
            ],
            "Column2": [
                "<h1>Title!!!! aby żeby ?! > ? ! zeby [ { ( ) } ].</h1>",
                "Another !<!i>paragraph {bla bla (bla!) [bla bla]} </i>",
            ],
        }
    )

    df = pre_processing(test_df)

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
    assert df.equals(test_df_after_pre_processing)
