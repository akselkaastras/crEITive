CODE = code
BIN = bin
RADSYM = radsym
FAST_RADSYM = fast_radsym
PCC = pcc
FAST_PCC = fast_pcc
PCC_CMPLX = pcc_cmplx
FAST_PCC_CMPLX = fast_pcc_cmplx
EIT = eit
FAST_EIT = fast_eit
EIT_CMPLX = eit_cmplx
FAST_EIT_CMPLX = fast_eit_cmplx
CGOS = cgos

FORWARD = $(RADSYM) $(FAST_RADSYM) $(PCC) $(FAST_PCC) $(PCC_CMPLX) $(FAST_PCC_CMPLX)
INVERSE = $(EIT) $(FAST_EIT) $(EIT_CMPLX) $(FAST_EIT_CMPLX) $(CGOS)
NON_FAST = $(RADSYM) $(PCC) $(PCC_CMPLX) $(EIT) $(EIT_CMPLX) $(CGOS)
FAST = $(FAST_RADSYM) $(FAST_PCC) $(FAST_PCC_CMPLX) $(FAST_EIT) $(FAST_EIT_CMPLX)
NON_CMPLX = $(RADSYM) $(FAST_RADSYM) $(PCC) $(FAST_PCC) $(EIT) $(FAST_EIT) $(CGOS)
CMPLX = $(PCC_CMPLX) $(FAST_PCC_CMPLX) $(EIT_CMPLX) $(FAST_EIT_CMPLX)
PROJECT = $(FORWARD) $(INVERSE)

INSTALL_FORWARD = $(addprefix install_,$(FORWARD))
INSTALL_INVERSE = $(addprefix install_,$(INVERSE))
INSTALL_NON_FAST = $(addprefix install_,$(NON_FAST))
INSTALL_FAST = $(addprefix install_,$(FAST))
INSTALL_NON_CMPLX = $(addprefix install_,$(NON_CMPLX))
INSTALL_CMPLX = $(addprefix install_,$(CMPLX))
INSTALL_PROJECT = $(INSTALL_FORWARD) $(INSTALL_INVERSE)

DEINSTALL_FORWARD = $(addprefix deinstall_,$(FORWARD))
DEINSTALL_INVERSE = $(addprefix deinstall_,$(INVERSE))
DEINSTALL_NON_FAST = $(addprefix deinstall_,$(NON_FAST))
DEINSTALL_FAST = $(addprefix deinstall_,$(FAST))
DEINSTALL_NON_CMPLX = $(addprefix deinstall_,$(NON_CMPLX))
DEINSTALL_CMPLX = $(addprefix deinstall_,$(CMPLX))
DEINSTALL_PROJECT = $(DEINSTALL_FORWARD) $(DEINSTALL_INVERSE)

CLEAN_FORWARD = $(addprefix clean_,$(FORWARD))
CLEAN_INVERSE = $(addprefix clean_,$(INVERSE))
CLEAN_NON_FAST = $(addprefix clean_,$(NON_FAST))
CLEAN_FAST = $(addprefix clean_,$(FAST))
CLEAN_NON_CMPLX = $(addprefix clean_,$(NON_CMPLX))
CLEAN_CMPLX = $(addprefix clean_,$(CMPLX))
CLEAN_PROJECT = $(CLEAN_FORWARD) $(CLEAN_INVERSE)

MRPROPER_FORWARD = $(addprefix mrproper_,$(FORWARD))
MRPROPER_INVERSE = $(addprefix mrproper_,$(INVERSE))
MRPROPER_NON_FAST = $(addprefix mrproper_,$(NON_FAST))
MRPROPER_FAST = $(addprefix mrproper_,$(FAST))
MRPROPER_NON_CMPLX = $(addprefix mrproper_,$(NON_CMPLX))
MRPROPER_CMPLX = $(addprefix mrproper_,$(CMPLX))
MRPROPER_PROJECT = $(MRPROPER_FORWARD) $(MRPROPER_INVERSE)

.PHONY: all install deinstall clean mrproper $(PROJECT) $(INSTALL_PROJECT) $(DEINSTALL_PROJECT) $(CLEAN_PROJECT) $(MRPROPER_PROJECT)

all: $(PROJECT)

install: $(INSTALL_PROJECT)

deinstall: $(DEINSTALL_PROJECT)

clean: $(CLEAN_PROJECT)

mrproper: $(MRPROPER_PROJECT)

forward: $(FORWARD)

inverse: $(INVERSE)

non_fast: $(NON_FAST)

fast: $(FAST)

non_cmplx: $(NON_CMPLX)

cmplx: $(CMPLX)

install_forward: $(INSTALL_FORWARD)

install_inverse:  $(INSTALL_INVERSE)

install_non_fast: $(INSTALL_NON_FAST)

install_fast: $(INSTALL_FAST)

install_non_cmplx: $(INSTALL_NON_CMPLX)

install_cmplx: $(INSTALL_CMPLX)

deinstall_forward: $(DEINSTALL_FORWARD)

deinstall_inverse:  $(DEINSTALL_INVERSE)

deinstall_non_fast: $(DEINSTALL_NON_FAST)

deinstall_fast: $(DEINSTALL_FAST)

deinstall_non_cmplx: $(DEINSTALL_NON_CMPLX)

deinstall_cmplx: $(DEINSTALL_CMPLX)

clean_forward: $(CLEAN_FORWARD)

clean_inverse: $(CLEAN_INVERSE)

clean_non_fast: $(CLEAN_NON_FAST)

clean_fast: $(CLEAN_FAST)

clean_non_cmplx: $(CLEAN_NON_CMPLX)

clean_cmplx: $(CLEAN_CMPLX)

mrproper_forward: $(MRPROPER_FORWARD)

mrproper_inverse: $(MRPROPER_INVERSE)

mrproper_non_fast: $(MRPROPER_NON_FAST)

mrproper_fast: $(MRPROPER_FAST)

mrproper_non_cmplx: $(MRPROPER_NON_CMPLX)

mrproper_cmplx: $(MRPROPER_CMPLX)

$(PROJECT):
	$(MAKE) -C $(CODE)/$@

$(INSTALL_PROJECT):
	cp -f $(CODE)/$(subst install_,,$@)/$(subst install_,,$@) ./$(subst install_,,$@)
	mv ./$(subst install_,,$@) $(BIN)

$(DEINSTALL_PROJECT):
	rm -f $(BIN)/$(subst deinstall_,,$@)

$(CLEAN_PROJECT):
	$(MAKE) -C $(CODE)/$(subst clean_,,$@) clean

$(MRPROPER_PROJECT):
	$(MAKE) -C $(CODE)/$(subst mrproper_,,$@) mrproper
