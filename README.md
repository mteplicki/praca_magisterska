# Biblioteka procedur szeregowania zadań w warunkach niepewności budżetowej `MasterThesis`

## Opis

Biblioteka zawiera implementacje algorytmów rozwiązujących problem szeregowania zadań w warunkach niepewności budżetowej (ang. *scheduling under budgeted uncertainty*). W bibliotece zaimplementowane są następujące algorytmy:
- algorytmy dla problemu szeregowania zadań na pojedynczej maszynie przy niepewnych czasach
  - przy kryterium minimalizacji łącznego opóżnienia (ang. *tardiness*) w klasycznym podejściu i przy użyciu reguł dominacji,
  - przy kryterium minimalizacji liczby opóźnionych zadań (ang. *tardy jobs*),
- algorytm dla problemu szeregowania zadań na równoległych niezależnych maszynach przy kryterium minimalizacji momentu zakończenia (ang. *makespan*),
- algorytm dla problemu szeregowania zadań dla permutacyjnego odpornego problemu przepływowego na dwóch maszynach (ang. *two-machine robust permutational flow shop*).

Część z nich to implementacja i usprawnienie algorytmów literaturowych, a jeden z nich (szeregowania zadań na równoległych niezależnych maszynach) jest zaprojektowany przez autora pracy. Wszystkie te algorytmy są oparte na frameworku *column and constraint generation*.


## Instalacja

Biblioteka wymaga zainstalowania języka Julia w wersji 1.11.5 lub nowszej.

Istnieją dwie opcje instalacji biblioteki `MasterThesis`. Pierwsza z nich to instalacja bezpośrednio z repozytorium GitHub:
```julia
import Pkg
Pkg.add("https://github.com/mteplicki/praca_magisterska")
using MasterThesis
```

Druga opcja to zbudowanie biblioteki z plików źródłowych. W tym celu należy pobrać repozytorium i przejść do głównego katalogu projektu. Następnie należy uruchomić REPL Julii i wykonać następujące polecenia:
```bash
julia --project
```
W REPLu należy zaimportować pakiet `Pkg` i wykonać następujące polecenia:
```julia
import Pkg
Pkg.instantiate()
using MasterThesis
```

## Przykładowe użycie

```julia
using MasterThesis
using CPLEX
instance = MultiUnrelatedInstance(50,4,10.0)
model = MultiUnrelatedMakespan(CPLEX.Optimizer, instance)
stats = column_generation(model)
println(stats)
```
