# nlp-project
<h1> Title generator </h1>

<h2><em>Celem tego projektu było utworzenie modelu LSTM (long-short term memory), który na podstawie streszczeń artykułów będzie generował ich tytuły.</em></h2>
<p>Z bazy PubMed pobrano tytuły i abstrakty publikacji naukowych, następnie poddano je preprocessingowi, a następnie po podzieleniu ich na zbiory testowy i treningowy, trenowano model.</p>
<p>Projekt zawierał wykorzystanie platformy <a href="https://www.kaggle.com">kaggle.com</a> do budowy i trenowania modelu (w trzech różnych wariantach - dla uproszczenia nie zawarto ich wszystkich w niniejszym podsumowaniu) przy wykorzystaniu GPU P100.</p><br>
<h3>Architektura projektu:</h3>
<ul>
    <li>data:
        <ul>
            <li>abstracts_and_titles.csv - zawiera pobrane z bazy PubMed streszczenia (w cudzysłowie) oraz tytuły artykułów</li>
            <li>polish.stopwords.txt - lista stopwords w języku polskim</li>
            <li>pubmed_keywords.py - lista keywords jako zmienna języka Python wykorzystywanych do przeszukiwania bazy danych PubMed</li>
        </ul>
    </li>
    <li>test: zawiera plik test_preprocessing.py, czyli zdefiniowaną funkcję sprawdzającą poprawność funkcji preprocessing</li>
    <li>title_generator:
        <ul>
            <li>_init_.py</li>
            <li>data_fetch.py - funkcja do pobrania streszczeń i tytułów artykułów z bazy danych PubMed</li>
            <li>preprocessing.py - obróbka danych obejmująca: zamianę na małe litery, usunięcie tagów html, tekstów wewnątrz nawiasów i samych nawiasów, cyfr, znaków interpunkcyjnych oraz polskich stopwords, a także tokenizację</li>
            <li>translate_to_polish.py - przetłumaczenie pobranych danych z języka angielskiego na język polski</li>
        </ul>
    </li>
    <li>project-nlp.ipynb - główna część projektu, która przede wszystkim zawiera model, jego trenowanie i komentarz do jego wyników</li>
    <li>req.txt - plik zawierający informację o wersjach używanych bibliotek w języku Python</li>
</ul>
