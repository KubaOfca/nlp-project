import os
import re
import pandas as pd
import string
import spacy


def pre_processing(df: pd.DataFrame) -> pd.DataFrame:
    pipeline = [
        lower_case_preprocessing,
        html_tags_preprocessing,
        https_pre_processing,
        brackets_preprocessing,
        polish_words_tokenize,
        punctuation_preprocessing,
        polish_stopwords_preprocessing,
    ]

    df_pre_processing = df.copy()
    for step in pipeline:
        try:
            df_pre_processing = step(df_pre_processing)
        except Exception as e:
            print("[PREPROCESSING] ", f"'{step.__name__}'", "error! (omitted)")
            print(e)
        else:
            print("[PREPROCESSING] ", f"'{step.__name__}'", "done!")
    return df_pre_processing


def lower_case_preprocessing(df: pd.DataFrame) -> pd.DataFrame:
    return df.applymap(lambda x: x.lower() if isinstance(x, str) else x)


def html_tags_preprocessing(df: pd.DataFrame) -> pd.DataFrame:
    return df.applymap(
        lambda x: remove_html_tag(x) if isinstance(x, str) else x
    )


def remove_html_tag(text: str) -> str:
    pattern = r"<.*?>"
    matched_html_tags = re.findall(pattern, text)
    for html_tag in matched_html_tags:
        text = text.replace(html_tag, "")
    return text


def https_pre_processing(df: pd.DataFrame) -> pd.DataFrame:
    return df.applymap(lambda x: remove_https(x) if isinstance(x, str) else x)


def remove_https(text: str) -> str:
    pattern = r"https:.*?\s"
    matched_https = re.findall(pattern, text)
    for html_tag in matched_https:
        text = text.replace(html_tag, "")
    return text


def brackets_preprocessing(df: pd.DataFrame) -> pd.DataFrame:
    return df.applymap(
        lambda x: remove_text_inside_brackets(x) if isinstance(x, str) else x
    )


def remove_text_inside_brackets(text: str) -> str:
    for brackets_pattern in [
        r"\([^()]*\)",
        r"\[[^\[\]]*\]",
        r"{[^{}]*\}",
    ]:
        matched_brackets = re.findall(brackets_pattern, text)
        for match in matched_brackets:
            text = text.replace(match, "")
    return text


def polish_words_tokenize(df: pd.DataFrame) -> pd.DataFrame:
    tokenizer = spacy.load("pl_core_news_sm")
    return df.applymap(
        lambda x: list(map(str, tokenizer(x))) if isinstance(x, str) else x
    )


def punctuation_preprocessing(df):
    return df.applymap(
        lambda x: remove_punctuation(x) if isinstance(x, list) else x
    )


def remove_punctuation(data: list):
    translator = str.maketrans(dict.fromkeys(string.punctuation))
    result = []
    for text in data:
        if text.translate(translator).strip():  # if text not empty
            result.append(text.translate(translator).strip())
    return result


def polish_stopwords_preprocessing(df):
    with open(
        os.path.join(os.getcwd(), "data", "polish.stopwords.txt"), "r"
    ) as file:
        stopwords = [line.strip() for line in file]
    return df.applymap(
        lambda x: list(filter(lambda y: y not in stopwords, x))
        if isinstance(x, list)
        else x
    )
