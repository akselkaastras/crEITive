-include ../arch.make

SRC = src
INCLUDE = include
OBJ = obj
BIN = cgos

CXXFLAGS += -I$(INCLUDE)
#LDFLAGS = -lgsl -lgslcblas -llapack -lf77blas

SRCS = $(wildcard $(SRC)/*.cpp)
INCLUDES = $(wildcard $(INCLUDE)/*.h)
OBJS = $(addprefix $(OBJ)/,$(patsubst %.cpp,%.o,$(notdir $(SRCS))))

.PHONY: all clean mrproper

all: $(BIN)

$(OBJS): | $(OBJ)

$(OBJ):
	mkdir $(OBJ)

$(OBJ)/cgos.o: $(filter-out $(INCLUDE)/cgos.h,$(INCLUDES))

$(OBJ)/h_zeta_matrix.o: $(INCLUDE)/h_zeta.h

$(OBJ)/cgos_traces_0.o: $(INCLUDE)/zeta_rhs.h $(INCLUDE)/psi_zeta.h

$(OBJ)/cgos_traces.o: $(INCLUDE)/h_zeta_matrix.h $(INCLUDE)/zeta_rhs.h $(INCLUDE)/psi_zeta.h

$(OBJ)/%.o: $(SRC)/%.cpp $(INCLUDE)/cgos.h
	$(CXX) -o $@ -c $< $(CXXFLAGS)

$(BIN): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

clean:
	@rm -f $(OBJ)/*

mrproper: clean
	@rm -f $(BIN)
