#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <unordered_map>
using namespace std;

struct Wierzcholek {
    int nr_sek;  
    int pozycja;      
    string podciag;   
};

vector<vector<int>> tworzenie_grafu(const vector<vector<pair<char, int>>>& sekwencje, int dl_podciag, vector<Wierzcholek>& wierzcholki) {
    wierzcholki.clear();
    for (int ind = 0; ind < (int)sekwencje.size(); ++ind) {
        const auto& sekwencja = sekwencje[ind];
        for (int i = 0; i + dl_podciag <= (int)sekwencja.size(); ++i) {
            string podciag = "";
            for (int j = 0; j < dl_podciag; ++j) {
                podciag += sekwencja[i + j].first;
            }
            wierzcholki.push_back({ ind, sekwencja[i].second, podciag });
        }
    }

    vector<vector<int>> graf(wierzcholki.size());
    unordered_map<string, vector<int>> grupy;

    for (int i = 0; i < (int)wierzcholki.size(); ++i) {
        grupy[wierzcholki[i].podciag].push_back(i);
    }

    for (auto const& [podciag, indeksy] : grupy) {
        for (size_t i = 0; i < indeksy.size(); ++i) {
            for (size_t j = i + 1; j < indeksy.size(); ++j) {
                int u = indeksy[i];
                int v = indeksy[j];
                if (wierzcholki[u].nr_sek != wierzcholki[v].nr_sek && 
                    abs(wierzcholki[u].pozycja - wierzcholki[v].pozycja) <= 10 * dl_podciag) {
                    graf[u].push_back(v);
                    graf[v].push_back(u);
                }
            }
        }
    }
    return graf;
}

bool czy_klika(const vector<int>& klika, const vector<vector<int>>& graf, const vector<Wierzcholek>& wierzcholki, int liczba_sekwencji) {
    vector<int> sprawdz(5, 0);
    for (int i = 0; i < klika.size(); i++) {
        sprawdz[wierzcholki[klika[i]].nr_sek] = 1;
    }
    for (int i = 0; i < sprawdz.size(); i++) {
        if (sprawdz[i] == 0) {
            return false;
        }
    }
    for (int i = 0; i < klika.size(); ++i) {
        for (int j = i + 1; j < klika.size(); ++j) {
            if (wierzcholki[i].nr_sek != wierzcholki[j].nr_sek) {
                if (find(graf[klika[i]].begin(), graf[klika[i]].end(), klika[j]) == graf[klika[i]].end()) {
                    return false;
                }
            }
        }
    }
    return true;
}

vector<int> znajdz_klike(vector<vector<int>> graf, const vector<Wierzcholek>& wierzcholki, int liczba_sekwencji) {
    vector<int> kolejnosc;
    for (int i = 0; i < (int)graf.size(); ++i) {
        if (graf[i].size() >= (size_t)liczba_sekwencji - 1) {
            kolejnosc.push_back(i);
        }
    }

    sort(kolejnosc.begin(), kolejnosc.end(), [&](int a, int b) {
        return graf[a].size() > graf[b].size();
    });

    for (int v : kolejnosc) {
        vector<int> kandydaci = graf[v];
        kandydaci.push_back(v);

        if (czy_klika(kandydaci, graf, wierzcholki, liczba_sekwencji)) {
            return kandydaci;
        }
    }
    return {};
}

vector<string> wczytaj_sekwencje(const string& nazwaPliku) {
    ifstream plik(nazwaPliku);
    if (!plik.is_open()) return {};
    string linia;
    vector<string> sekwencje;
    while (getline(plik, linia)) {
        if (linia.empty()) continue; 
        while (!linia.empty() && (linia.back() == '\r' || linia.back() == '\n')) linia.pop_back();
        if (linia[0] != '>') sekwencje.push_back(linia);
    }
    return sekwencje;
}

vector<vector<int>> wczytaj_qual(const string& nazwaPliku) {
    ifstream plik(nazwaPliku);
    if (!plik.is_open()) return {};
    string linia;
    vector<vector<int>> qual;
    while (getline(plik, linia)) {
        if (!linia.empty() && linia[0] == '>') { 
            qual.push_back(vector<int>()); 
        } else if (!qual.empty()) {
            stringstream ss(linia); 
            int liczba;
            while (ss >> liczba) qual.back().push_back(liczba); 
        }
    }
    return qual;
}

void usun_nukleotydy(vector<vector<pair<char, int>>>& sekwencje, vector<vector<pair<int, int>>>& quals, int prog) {
    for (size_t i = 0; i < sekwencje.size(); ++i) {
        auto& seq = sekwencje[i];
        auto& ql = quals[i];
        vector<pair<char, int>> new_seq;
        vector<pair<int, int>> new_ql;
        for (size_t j = 0; j < seq.size(); ++j) {
            if (ql[j].first >= prog) {
                new_seq.push_back(seq[j]);
                new_ql.push_back(ql[j]);
            }
        }
        seq = new_seq;
        quals[i] = new_ql;
    }
}

vector<vector<pair<char, int>>> pozycjonowanie_fasta(const vector<string>& sekwencje_wczytane) {
    vector<vector<pair<char, int>>> wynik; 

    for (const auto& sekwencja : sekwencje_wczytane) {
        vector<pair<char, int>> pozycje; 
        for (int i = 0; i < sekwencja.size(); ++i) { 
            pozycje.emplace_back(sekwencja[i], i+1);   
        }
        wynik.push_back(pozycje);  
    }

    return wynik; 
}

vector<vector<pair<int, int>>> pozycjonowanie_qual(const vector<vector<int>>& qual) {
    vector<vector<pair<int, int>>> wynik; 

    for (const auto& qual_sekwencja : qual) { 
        vector<pair<int, int>> pozycje; 
        for (int i = 0; i < qual_sekwencja.size(); ++i) {
            pozycje.emplace_back(qual_sekwencja[i], i+1); 
        }
        wynik.push_back(pozycje);  
    }

    return wynik; 
}

int main() {
    int prog, oligo;
    cout << "Wprowadz prog jakosci: "; cin >> prog;
    cout << "Wprowadz oczekiwaną długość motywu: "; cin >> oligo;
    string plik_fasta = "sekwencje.fasta";  
    string plik_qual = "sekwencje.qual";    

    vector<string> sekwencje_wczytane = wczytaj_sekwencje(plik_fasta);
    vector<vector<int>> qual_wczytane = wczytaj_qual(plik_qual); 

    if (sekwencje_wczytane.empty() || qual_wczytane.empty()) {
        cout << "Blad: Nie udalo sie wczytac plikow!" << endl;
        return 1;
    }

    if (sekwencje_wczytane.size() != qual_wczytane.size()) {
        cout << "Blad: Liczba sekwencji w pliku FASTA nie odpowiada liczbie sekwencji w pliku QUAL!" << endl;
        return 1;
    }   

    vector<vector<pair<char, int>>> sekwencje = pozycjonowanie_fasta(sekwencje_wczytane);
    vector<vector<pair<int, int>>> qual = pozycjonowanie_qual(qual_wczytane);

    usun_nukleotydy(sekwencje, qual, prog);

    cout << "Wprowadzone sekwencje po usunieciu nukleotydow o jakosci ponizej progu: " << endl;
    for (int i = 0; i < sekwencje.size(); i++) {
        cout << i + 1 << ".";
        for (int j = 0; j < sekwencje[i].size(); j++) {
            cout << sekwencje[i][j].first;
        }
        cout << endl;
        cout << sekwencje[i].size() << endl;
        for (int j = 0; j < qual[i].size(); j++) {
            cout << qual[i][j].first << " ";
        }
        cout << endl;
        cout << qual[i].size() << endl;
    }

    int liczba_sek = sekwencje.size();
    vector<Wierzcholek> wierzcholki;
    vector<vector<int>> graf = tworzenie_grafu(sekwencje, oligo, wierzcholki);
    vector<int> klika = znajdz_klike(graf, wierzcholki, liczba_sek);

    if (klika.empty()){
        cout << "Nie znaleziono motywu we wszystkich sekwencjach" << endl;
    }
    else {
        cout << "Znaleziono motyw: " << endl;
        for (int v : klika) {
            cout << wierzcholki[v].nr_sek + 1 << " " << wierzcholki[v].pozycja << " " << wierzcholki[v].podciag << endl;
        }
        cout << endl;
    }
}