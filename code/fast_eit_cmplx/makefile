-include ../arch.make

SRC = src
INCLUDE = include
OBJ = obj
BIN = fast_eit_cmplx

CXXFLAGS += -I$(INCLUDE) -Xclang -fopenmp
LDFLAGS += -lomp

SRCS = $(wildcard $(SRC)/*.cpp)
INCLUDES = $(wildcard $(INCLUDE)/*.h)
OBJS = $(addprefix $(OBJ)/,$(patsubst %.cpp,%.o,$(notdir $(SRCS))))

.PHONY: all clean mrproper

all: $(BIN)

$(OBJS): | $(OBJ)

$(OBJ):
	mkdir $(OBJ)

$(OBJ)/eit.o: $(filter-out $(INCLUDE)/eit.h,$(INCLUDES))

$(OBJ)/h_zeta_matrix.o: $(INCLUDE)/h_zeta.h

$(OBJ)/qhat_zeta_exp.o: $(INCLUDE)/zeta_rhs.h $(INCLUDE)/exp_minus_i_x_xi.h $(INCLUDE)/scattering_transform.h

$(OBJ)/qhat_zeta_0.o: $(INCLUDE)/zeta_rhs.h $(INCLUDE)/exp_minus_i_x_xi.h $(INCLUDE)/psi_zeta.h $(INCLUDE)/scattering_transform.h

$(OBJ)/qhat_zeta.o: $(INCLUDE)/h_zeta_matrix.h $(INCLUDE)/zeta_rhs.h $(INCLUDE)/exp_minus_i_x_xi.h $(INCLUDE)/psi_zeta.h $(INCLUDE)/scattering_transform.h

$(OBJ)/%.o: $(SRC)/%.cpp $(INCLUDE)/eit.h
	$(CXX) -o $@ -c $< $(CXXFLAGS)

$(BIN): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

clean:
	@rm -f $(OBJ)/*

mrproper: clean
	@rm -f $(BIN)
