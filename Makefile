# Compiladores e opções
CXX      = g++
CC       = gcc
CXXFLAGS = -O2 -march=native -std=c++17 -Wall
CFLAGS   = -O2 -Wall
LDLIBS   = -llapacke -lopenblas

BUILD  = build
COMMON = $(BUILD)/mmio.o $(BUILD)/matrix.o

# Alvos que são comandos, e não nomes de arquivos
.PHONY: all check parte1 testar parte3 spy spy_rcm clean

all: $(BUILD)/check_env $(BUILD)/parte1 $(BUILD)/testar_metodos $(BUILD)/parte3

# Nas regras abaixo, cada variável automática significa:
#   $@  o arquivo que está sendo gerado (à esquerda do ":")
#   $<  a primeira dependência
#   $^  todas as dependências
# O que vem depois de "|" só precisa existir: a pasta build/ não força recompilação.

$(BUILD)/mmio.o: src/mmio.c src/mmio.h | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/matrix.o: src/matrix.cpp src/matrix.hpp src/mmio.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -Isrc -c $< -o $@

$(BUILD)/check_env: src/check_env.cpp $(BUILD)/mmio.o | $(BUILD)
	$(CXX) $(CXXFLAGS) -Isrc $^ -o $@ $(LDLIBS)

$(BUILD)/parte1: src/parte1.cpp src/matrix.hpp $(COMMON) | $(BUILD)
	$(CXX) $(CXXFLAGS) -Isrc src/parte1.cpp $(COMMON) -o $@ $(LDLIBS)

$(BUILD)/parte2.o: src/parte2.cpp src/parte2.hpp src/matrix.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -Isrc -c $< -o $@

$(BUILD)/testar_metodos: src/testar_metodos.cpp src/parte2.hpp src/matrix.hpp $(BUILD)/parte2.o $(COMMON) | $(BUILD)
	$(CXX) $(CXXFLAGS) -Isrc src/testar_metodos.cpp $(BUILD)/parte2.o $(COMMON) -o $@ $(LDLIBS)

$(BUILD)/rcm.o: src/rcm.cpp src/rcm.hpp src/matrix.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -Isrc -c $< -o $@

$(BUILD)/parte3: src/parte3.cpp src/rcm.hpp src/parte2.hpp src/matrix.hpp $(BUILD)/rcm.o $(BUILD)/parte2.o $(COMMON) | $(BUILD)
	$(CXX) $(CXXFLAGS) -Isrc src/parte3.cpp $(BUILD)/rcm.o $(BUILD)/parte2.o $(COMMON) -o $@ $(LDLIBS)

$(BUILD):
	mkdir -p $(BUILD)

# Verifica o ambiente (BLAS em uma única thread, como pede o enunciado)
check: $(BUILD)/check_env
	OPENBLAS_NUM_THREADS=1 ./$(BUILD)/check_env

# Parte 1: tabela de caracterização em results/parte1.csv
parte1: $(BUILD)/parte1
	OPENBLAS_NUM_THREADS=1 ./$(BUILD)/parte1

# Parte 2: testa todos os métodos contra a LAPACK em matrices/teste_10x10.mtx
testar: $(BUILD)/testar_metodos
	OPENBLAS_NUM_THREADS=1 ./$(BUILD)/testar_metodos

# Parte 3: bandas antes e depois do RCM e tempos densa x banda. Uso: make parte3 MATRIZ=matrices/bcsstk14.mtx
MATRIZ = matrices/teste_10x10.mtx
parte3: $(BUILD)/parte3
	OPENBLAS_NUM_THREADS=1 ./$(BUILD)/parte3 $(MATRIZ)

# Parte 3: padrões de esparsidade lado a lado em figures/ (rode antes make parte3)
spy_rcm:
	conda run -n tp-alc python scripts/spy_rcm.py $(MATRIZ)

# Parte 1: spy plots em figures/
spy:
	conda run -n tp-alc python scripts/spy.py

clean:
	rm -rf $(BUILD)
