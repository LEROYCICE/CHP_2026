# ============================================================
#  Makefile  –  Projet CHP : Advection 2D par differences finies
# ============================================================

CXX      = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++11 -I Include

SRC_DIR  = src
INC_DIR  = Include

SRCS = $(SRC_DIR)/main.cpp      \
       $(SRC_DIR)/init.cpp      \
       $(SRC_DIR)/fonctions.cpp \
       $(SRC_DIR)/Prod.cpp      \
       $(SRC_DIR)/BICgSTAB.cpp

OBJS = $(SRCS:.cpp=.o)
TARGET = exec

# ---- Cibles principales ----

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lm
	@echo "Compilation reussie : ./$(TARGET)"

# Compilation de chaque .cpp en .o
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp $(wildcard $(INC_DIR)/*.h)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---- Gnuplot ----
gp:
	gnuplot gnuplot.gp

# ---- Execution ----

run: $(TARGET)
	@mkdir -p Sortie Gnuplots
	./$(TARGET) para.dat

# ---- Nettoyage ----

clean:
	rm -f $(OBJS) $(TARGET)

distclean: clean
	rm -f Sortie/*.dat Gnuplots/*.png

# ---- Aide ----

help:
	@echo "Cibles disponibles :"
	@echo "  make          -> compile le projet"
	@echo "  make run      -> compile et execute avec para.dat"
	@echo "  make gp       -> lance gnuplot"
	@echo "  make clean    -> supprime les objets et l'executable"
	@echo "  make distclean -> supprime egalement les fichiers de sortie"

.PHONY: all run clean distclean help gp