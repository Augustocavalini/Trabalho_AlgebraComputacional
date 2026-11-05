# Compiladores e opções
CXX      = g++
CC       = gcc
CXXFLAGS = -O2 -march=native -std=c++17 -Wall
CFLAGS   = -O2 -Wall
LDLIBS   = -llapacke -lopenblas

BUILD  = build
COMMON = $(BUILD)/mmio.o $(BUILD)/matrix.o

# Alvos que são comandos, e não nomes de arquivos
.PHONY: all check parte1 spy clean

all: $(BUILD)/check_env $(BUILD)/parte1 $(BUILD)/metodos.o

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

$(BUILD)/metodos.o: src/metodos.cpp src/metodos.hpp src/matrix.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -Isrc -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

# Verifica o ambiente (BLAS em uma única thread, como pede o enunciado)
check: $(BUILD)/check_env
	OPENBLAS_NUM_THREADS=1 ./$(BUILD)/check_env

# Parte 1: tabela de caracterização em results/parte1.csv
parte1: $(BUILD)/parte1
	OPENBLAS_NUM_THREADS=1 ./$(BUILD)/parte1

# Parte 1: spy plots em figures/
spy:
	conda run -n tp-alc python scripts/spy.py

clean:
	rm -rf $(BUILD)
