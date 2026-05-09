# 🧬 Wyszukiwarka wspólnego motywu w sekwencjach: podejście oparte na grafach

### Algorytm wyszukiwania wspólnych motywów DNA w sekwencjach biologicznych z uwzględnieniem jakości odczytu (Phred score).

## 📖 O projekcie
Ten program w C++ służy do identyfikacji powtarzających się podciągów (motywów) w wielu sekwencjach DNA jednocześnie. Projekt wykorzystuje teorię grafów, reprezentując potencjalne motywy jako wierzchołki, a ich podobieństwo jako krawędzie. Celem jest znalezienie **kliki**, która reprezentuje motyw obecny w każdym z badanych łańcuchów.

## 🚀 Główne funkcjonalności
- **Filtr Jakości (Quality Control):** Program wczytuje pliki `.qual` i usuwa nukleotydy o niskiej wiarygodności odczytu (poniżej progu zdefiniowanego przez użytkownika).
- **Pozycjonowanie nukleotydów:** Zachowuje oryginalną numerację pozycji w sekwencji, nawet po usunięciu słabych ogniw.
- **Optymalizacja grafu (Hash Maps):** Budowa grafu wykorzystuje mapy hashujące (`unordered_map`), co przyspiesza łączenie identycznych oligonukleotydów.
- **Detekcja Kliki:** Algorytm szuka zestawu połączonych wierzchołków, z których każdy pochodzi z innej sekwencji wejściowej.
- **Walidacja biologiczna:** Uwzględnia dystans fizyczny między nukleotydami (parametr `10 * dl_podciag`).

## 🛠️ Technologie
- **Język:** C++
- **Struktury danych:** Grafy (listy sąsiedztwa), Wektory, Mapy asocjacyjne.
- **Algorytmy:** Sortowanie oparte na stopniu wierzchołka, znajdowanie kliki w grafie rzadkim.

## 📋 Struktura plików wejściowych
Program wymaga dwóch plików w tym samym katalogu:
1. `sekwencje.fasta` – plik z sekwencjami DNA w formacie FASTA (dostępny przykładowy)
2. `sekwencje.qual` – plik z wartościami jakości (Phred scores) dla każdej bazy (dostępny przykładowy)

## ⚙️ Kompilacja i uruchomienie

### Kompilacja:
Użyj kompilatora `g++` (zalecana flaga optymalizacji `-O3` dla dużych zbiorów danych):
```bash
g++ -O3 main.cpp -o motif_finder