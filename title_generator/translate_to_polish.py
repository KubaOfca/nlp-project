from deep_translator import GoogleTranslator

batch_size = 4999


def translate_to_polish(text):
    translated_text = []
    for i in range(0, len(text), batch_size):
        translated_text.append(
            GoogleTranslator(source="english", target="polish").translate(
                text[i : i + batch_size]
            )
        )
    return "".join(translated_text)
