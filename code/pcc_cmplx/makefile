-include ../arch.make

SRC = src
INCLUDE = include
OBJ = obj
BIN = pcc_cmplx

CXXFLAGS += -I$(INCLUDE)
#LDFLAGS = -llapack -lf77blas

SRCS =$(wildcard $(SRC)/*.cpp)
INCLUDES =$(wildcard $(INCLUDE)/*.h)
OBJS = $(addprefix $(OBJ)/,$(patsubst %.cpp,%.o,$(notdir $(SRCS))))

.PHONY: all clean mrproper

all: $(BIN)

$(OBJS): | $(OBJ)

$(OBJ):
	mkdir $(OBJ)

$(OBJ)/pcc.o: $(filter-out $(INCLUDE)/pcc.h,$(INCLUDES))

$(OBJ)/conductivity.o: $(INCLUDE)/gauss_legendre.h $(INCLUDE)/phi_values.h $(INCLUDE)/associated_legendre_functions.h $(INCLUDE)/exp_i_m_phi.h

$(OBJ)/kp_matrix.o: $(INCLUDE)/associated_legendre_functions.h $(INCLUDE)/exp_i_m_phi.h

$(OBJ)/kp_kbarp_matrix.o: $(INCLUDE)/kp_matrix.h $(INCLUDE)/kpcross_matrix.h $(INCLUDE)/kbarp_matrix.h

$(OBJ)/rhs_matrix.o: $(INCLUDE)/grad_solid_harmonics.h

$(OBJ)/rhs_block_matrix.o: $(INCLUDE)/rhs_matrix.h

$(OBJ)/%.o: $(SRC)/%.cpp $(INCLUDE)/pcc.h
	$(CXX) -o $@ -c $< $(CXXFLAGS)

$(BIN): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

clean:
	@rm -f $(OBJ)/*

mrproper: clean
	@rm -f $(BIN)
