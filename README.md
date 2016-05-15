# feedback-vertex-set

Kompiliert mit: g++ -Og -Wall -std=c++11 -I "/home/fabian/Documents/Progamming/boost_1_60_0/" fvs_solver.cpp main.cpp -o fvs
Eine Warnung.

# An alle:

Es wäre schön wenn ab jetzt nur noch changes gepusht würden die verwendbar und getestet sind. Damit können dann auch alle anderen gleichzeitig oder nach dem push noch weiterarbeiten (da das Projekt noch compiled).

Grundsätzlich wäre es außerdem schön, wenn alle die boost::adjacency_list verstehen würden. Momentan ist in compute_fvs noch ein riesiger Bug (den ich gerade nicht fixen konnte, weil das Repo haywire gegangen ist), der mit der Datenstruktur zu tun hat. Vor allem remove_vertex sollte nicht oder nur sehr vorsichtig verwendet werden, da es alle (!) Indizes und Referenzen auf Knoten vernichtet, die momentan irgendwo außerhalb der Graphenstruktur rumschwirren. Bitte schaut euch an, was für Funktionen ihr verwendet und testet den Code (sry, das ich das nicht auf von Anfang an gemacht habe).

