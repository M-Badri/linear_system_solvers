## Copyright (c) 2019 alberto Otero de la Roza <aoterodelaroza@gmail.com>
## This file is frere software; distributed under GNU/GPL version 3.

## FC=compiler
## FCSYNTAX=syntax-only compilation flag
## FCMODDIR=flag to set read/write .mod and .smod directory
## FCMODREADDIR=flag to set write-only .mod and .smod directory
## MODDIR=.mod and .smod directory (leave blank for root)

PROJECT:=zoo
#FC:=gfortran
FC:=ifort

ifeq ($(FC),gfortran) 
	FCSYNTAX:=-fsyntax-only
	FCMODDIR:=-J
	FCMODREADDIR:=-I
	FCFLAGS:= -cpp -O2 -fno-range-check
else
	FCSYNTAX:=-syntax-only
	FCMODDIR:=-module
	FCMODREADDIR:=-I
	FCFLAGS:= -fpp -O2
endif

MODDIR:=mod/$(FC)
LIBDIR:=lib/$(FC)


#### user input ends here ####

## some tricks for text manipulation
null:=
space:=$(null) $(null)
$(space):=$(space)
define \n


endef

## no implicit rules
.SUFFIXES: 

## auxiliary programs
AWK:=awk
SED:=sed
RM:=rm -f
MKDIR:=mkdir -p
TEST:=test
PROJECT:=lib$(PROJECT)

## known fortran extensions
FORTEXT:=f F fpp FPP for FOR ftn FTN f90 F90 f95 F95 f03 F03 f08 F08

## locate the source files
SOURCES:=$(shell find . -regextype posix-awk -regex '.*\.($(subst $( ),|,$(FORTEXT)))$$')

## compilation and syntax-compilation commands
COMPILE.f08 = $(FC) $(FCFLAGS) $(TARGET_ARCH) -c
MAKEMOD.f08 = $(FC) $(FCFLAGS) $(TARGET_ARCH) $(FCSYNTAX) -c

## create the mod and smod directory; define slashed version of MODDIR
ifneq ($(MODDIR),)
  $(shell $(TEST) -d $(MODDIR) || $(MKDIR) -p $(MODDIR))
  MODDIRSLSH:=$(MODDIR)/
else
  MODDIRSLSH:=./
endif

## create the lib 
ifneq ($(LIBDIR),)
  $(shell $(TEST) -d $(LIBDIR) || $(MKDIR) -p $(LIBDIR))
endif

## create the temporary mod and smod directory
ifneq ($(FCMODREADDIR),)
  MODDIRTMP:=.tmp$(MODDIR)
  $(shell $(TEST) -d $(MODDIRTMP) || $(MKDIR) -p $(MODDIRTMP))
  MAKEMOD.f08+= $(FCMODDIR) $(MODDIR)
  COMPILE.f08+= $(FCMODREADDIR) $(MODDIR) $(FCMODDIR) $(MODDIRTMP)
else
  MAKEMOD.f08+= $(FCMODDIR) $(MODDIR)
  COMPILE.f08+= $(FCMODDIR) $(MODDIR)
endif

## define the anchors and the objects variables
# $(call source-to-extension,source-file-list,new-extension)
define source-to-extension
  $(strip \
    $(foreach ext,$(FORTEXT),\
      $(subst .$(ext),.$2,$(filter %.$(ext),$1))))
endef
OBJECTS:=$(call source-to-extension,$(SOURCES),o)
ANCHORS:=$(call source-to-extension,$(SOURCES),anc)

## default target, main and clean targets
all: $(PROJECT)

$(PROJECT): $(OBJECTS)
	ar -rv $@ $^
	$(shell $(TEST) -f *.i90 || mv *.i90 $(MODDIR))
	mv $(PROJECT) $(LIBDIR) 
	$(RM) $(OBJECTS) $(ANCHORS)
#	$(FC) -o $@ $+

.PHONY: clean
clean:
	-$(RM) *.mod *.smod $(OBJECTS) $(ANCHORS) main *.i90
	-$(RM) $(MODDIR)/*
#	-$(TEST) -d $(MODDIR) && $(RM) -r $(MODDIR)
#	-$(TEST) -d $(MODDIRTMP) && $(RM) -r $(MODDIRTMP)

## syntax-only compilation rule: all anchor files depend on their source
# $(call modsource-pattern-rule,extension)
define modsource-pattern-rule
%.anc: %.$1
	$$(MAKEMOD.f08) $$<
	@touch $$@
endef
$(foreach ext,$(FORTEXT),$(eval $(call modsource-pattern-rule,$(ext))))

## compilation rule: objects depend on their anchor file
%.o: %.anc
	$(COMPILE.f08) $(OUTPUT_OPTION) $(wildcard $(addprefix $*.,$(FORTEXT)))
ifdef MODDIRTMP
	-@$(RM) $(MODDIRTMP)/*.mod $(MODDIRTMP)/*.smod
endif
	@touch $@

## automatically generate the dependency rules
$(eval $(subst $( ),$(\n),$(shell $(AWK) --traditional -f makedepf08.awk $(SOURCES) | sort | uniq | $(SED) -e 's!^.mod/!$(MODDIRSLSH)!' -e 's!:.mod/!:$(MODDIRSLSH)!')))
