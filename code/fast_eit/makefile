# -include ../arch.make

# SRC = src
# INCLUDE = include
# OBJ = obj
# BIN = fast_eit

# CXXFLAGS += -I$(INCLUDE) -Xclang -fopenmp
# LDFLAGS += -lomp

# SRCS = $(wildcard $(SRC)/*.cpp)
# INCLUDES = $(wildcard $(INCLUDE)/*.h)
# OBJS = $(addprefix $(OBJ)/,$(patsubst %.cpp,%.o,$(notdir $(SRCS))))

# .PHONY: all clean mrproper

# all: $(BIN)

# $(OBJS): | $(OBJ)

# $(OBJ):
# 	mkdir $(OBJ)

# $(OBJ)/eit.o: $(filter-out $(INCLUDE)/eit.h,$(INCLUDES))

# $(OBJ)/h_zeta_matrix.o: $(INCLUDE)/h_zeta.h

# $(OBJ)/qhat_zeta_exp.o: $(INCLUDE)/zeta_rhs.h $(INCLUDE)/exp_minus_i_x_xi.h $(INCLUDE)/scattering_transform.h

# $(OBJ)/qhat_zeta_0.o: $(INCLUDE)/zeta_rhs.h $(INCLUDE)/exp_minus_i_x_xi.h $(INCLUDE)/psi_zeta.h $(INCLUDE)/scattering_transform.h

# $(OBJ)/qhat_zeta.o: $(INCLUDE)/h_zeta_matrix.h $(INCLUDE)/zeta_rhs.h $(INCLUDE)/exp_minus_i_x_xi.h $(INCLUDE)/psi_zeta.h $(INCLUDE)/scattering_transform.h

# $(OBJ)/%.o: $(SRC)/%.cpp $(INCLUDE)/eit.h
# 	$(CXX) -o $@ -c $< $(CXXFLAGS)

# $(BIN): $(OBJS)
# 	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

# clean:
# 	@rm -f $(OBJ)/*

# mrproper: clean
# 	@rm -f $(BIN)



# This is to avoid users overriding the default target
default: all

-include ../arch.make

# Internal variables used for folder names
SRC = src
INCLUDE = include
OBJ = obj
BIN = fast_eit

# Default flags
CXXFLAGS += -I$(INCLUDE)
# It is bad practice to force certain flags in the makefile
# This should be the users responsibility!
# Note that intel compilers prefer -qopenmp, so you can't even do -fopenmp
CXXFLAGS += -Xclang -fopenmp
# Same here... Let users control this, less error prone.
# Also, I don't think you need this. The -fopenmp should do this automatically.
# I also think clang do this. But perhaps it is really needed.
LIBS += -lomp

SRCS = $(wildcard $(SRC)/*.cpp)
INCLUDES = $(wildcard $(INCLUDE)/*.h)
OBJS = $(addprefix $(OBJ)/,$(patsubst %.cpp,%.o,$(notdir $(SRCS))))

.PHONY: all clean mrproper

all: $(BIN)

$(OBJS): | $(OBJ)

# mkdir -p ensures it won't fail if the obj folder already exists
$(OBJ):
	mkdir -p $(OBJ)

$(OBJ)/eit.o: $(filter-out $(INCLUDE)/eit.h,$(INCLUDES))
$(OBJ)/h_zeta_matrix.o: $(INCLUDE)/h_zeta.h
$(OBJ)/qhat_zeta_exp.o: $(INCLUDE)/zeta_rhs.h $(INCLUDE)/exp_minus_i_x_xi.h $(INCLUDE)/scattering_transform.h
$(OBJ)/qhat_zeta_0.o: $(INCLUDE)/zeta_rhs.h $(INCLUDE)/exp_minus_i_x_xi.h $(INCLUDE)/psi_zeta.h $(INCLUDE)/scattering_transform.h
$(OBJ)/qhat_zeta.o: $(INCLUDE)/h_zeta_matrix.h $(INCLUDE)/zeta_rhs.h $(INCLUDE)/exp_minus_i_x_xi.h $(INCLUDE)/psi_zeta.h $(INCLUDE)/scattering_transform.h

$(OBJ)/%.o: $(SRC)/%.cpp $(INCLUDE)/eit.h
	$(CXX) -o $@ -c $< $(CXXFLAGS)

$(BIN): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS) $(LDFLAGS)

clean:
	@rm -f $(OBJ)/*

mrproper: clean
	@rm -f $(BIN)
