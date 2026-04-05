# ============================================================
#  gnuplot.gp  –  Animation 3D GIF du transport
# ============================================================

# --- Parametres ---
if (!exists("pas")) pas = 4
if (!exists("fps")) fps = 12

# --- Terminal GIF anime ---
set terminal gif animate delay int(100/fps) optimize size 900,650
set output "Gnuplot/transport.gif"

# --- Style ---
set palette defined (0 '#000080', 0.4 '#0080ff', 0.7 '#00ffcc', 0.85 '#ffff00', 1 '#ff2000')
set pm3d
set hidden3d

set xlabel "x"
set ylabel "y"
set zlabel "u"
set title "Transport 2D par differences finies"

set xrange [-1:1]
set yrange [-1:1]
set zrange [-0.1:1.2]
set cbrange [0:1]

set ticslevel 0
set view 35, 310

# --- Recuperer automatiquement le nombre de fichiers ---
# Compte combien de fichiers sol.*.dat existent
N = int(system("ls Sortie/sol.*.dat 2>/dev/null | wc -l"))

print sprintf("Nombre de fichiers detectes : %d", N)

# --- Animation ---
do for [k=1:N:pas] {
    fichier = sprintf("Sortie/sol.%04d.dat", k)

    set label 1 sprintf("Iteration %04d", k) at screen 0.12,0.92

    splot fichier using 1:2:3 with pm3d notitle
}

set output
print "GIF sauvegarde dans Gnuplot/transport.gif"