Für Windows:
--------------------------------------------
2.        Quellcode
--------------------------------------------
-In main.cpp und THDQuantities.cpp sind alle Methoden enthalten
-"THDQuantities.h" ist die Header Datei von "THDQuantities.cpp" . Dort stehen einige
Beschreibungen zu den Methoden 
-"plotting.py" enthält den Python-Skript zur Generierung unserer Plots

------------------------------------------------
1. Starten einer Simulation ohne Kompilierung
------------------------------------------------
-Simulationen können sequenziell über run_simulations.bat (Batch file) gestartet werden

-die Batch files führen "Density1DRods.exe" aus
-Es sind jeweils Beispiele eingetragen, die Parameter haben im batch-file die folgende Reihenfolge:
<.exe> <rod_length> <Lattice_size> <max_manual_steps>

-die Ergebnisse sind dann in "results" mit dem Namen "L<rod_length>_rho.txt" und "L<rod_length>_gamma.txt" zu finden
-es werden 2 verschiedene Files generiert. In einem befinden sich die Dichteprofile im anderen die damit berechnetet Werte

-----------------------------------------------------------
3. Kompilieren und starten des Code über Visual Studio
-----------------------------------------------------------
-es wird Visual Studio (nicht Visual Studio Code) benötigt
-Es muss die Datei "Density1DRods_VS.sln" (Projektmappe) in Visual Studio geöffnet werden welches die Datei beim starten automatisch kompiliert
-Beim starten in können die Parameter beim rechtsklicken von Density1DRods_VS im Projektmappen-Explorer unter 
Eigenschaften -> Konfigurationseigenschaften -> Debugging -> Befehlsargumente angegeben werden:    
 <rod_length> <Lattice_size> <max_manual_steps>