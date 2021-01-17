# Makefile

# host toolchain
HCC	= gcc
HCFLAGS	= -O2 -std=c99 -I$(EPIPHANY_HOME)/tools/host/include -Wall

# Use the following flags for the old SDK ver 5.13.09.10
#HLFLAGS	= -L$(EPIPHANY_HOME)/tools/host/lib -le-hal 

# Use the following flags for the new SDK ver 2015.1
HLFLAGS	= -L$(EPIPHANY_HOME)/tools/host/lib -le-hal -le-loader -lm
ECHO	= /bin/echo -e

# target toolchain
ECC	= e-gcc
EOC	= e-objcopy
ECFLAGS	= -Os -std=c99 -falign-loops=8 -falign-functions=8 -Wall
ELFLAGS	= -T$(EPIPHANY_HOME)/bsps/current/fast.ldf -le-lib -lm
EOFLAGS	= -R .shared_dram

# folders
HSRC	= hsrc
HDEST	= hobj
ESRC	= esrc
EDEST	= eobj
DEST	= bin

# applications (host binary and epiphany SREC files) to build
HAPP	= $(DEST)/ep_main
EAPPS	= $(DEST)/testcode.srec

# object files to build
HOBJS	= $(HDEST)/arm16x16.o
#$(HDEST)/matrices.o 
EOBJS	= $(EDEST)/testcode.o
#$(EDEST)/householder.o $(EDEST)/matrices.o $(EDEST)/vectors.o

# === Rules ===============================================================
.SECONDARY:
.PHONY: help all host target folders run clean
EELFS	= $(EAPPS:$(DEST)%srec=$(EDEST)%elf)

help:
	@$(ECHO)
	@$(ECHO) "Epiphany Makefile - Help"
	@$(ECHO) "  help    show this help"
	@$(ECHO) "  host    build host application      ($(HAPP))"
	@$(ECHO) "  target  build epiphany applications ($(EAPPS))"
	@$(ECHO) "  all     build host and target"
	@$(ECHO) "  run     build all, then run host application"
	@$(ECHO) "  clean   remove applications and intermediate files"
	@$(ECHO)

all: host

host: folders $(HAPP)

target: folders $(EAPPS)

folders: $(HDEST) $(EDEST) $(DEST)

run: host
	@$(ECHO) "\tRUN"
	@sudo LD_LIBRARY_PATH=$(LD_LIBRARY_PATH) \
	      EPIPHANY_HDF=$(EPIPHANY_HDF) \
	      $(HAPP)

clean:
	@$(ECHO) "\tCLEAN"
	@rm -v -f $(HOBJS) $(EOBJS) $(EELFS) $(EAPPS) $(HAPP) $(EAPPS)
	@rmdir -v --ignore-fail-on-non-empty $(HDEST) $(EDEST) $(DEST)

$(HDEST):
	@$(ECHO) "\tMKDIR $(HDEST)"
	@mkdir -p $(HDEST)

$(EDEST):
	@$(ECHO) "\tMKDIR $(EDEST)"
	@mkdir -p $(EDEST)

$(DEST):
	@$(ECHO) "\tMKDIR $(DEST)"
	@mkdir -p $(DEST)

# === Host Toolchain ======================================================
$(HAPP): $(HOBJS)
	@$(ECHO) "\t(HOST)   LINK\t$@"
	@$(HCC) -o $@ $^ $(HLFLAGS)

$(HDEST)/%.o: %.c
	@$(ECHO) "\t(HOST)   CC\t$@"
	@$(HCC) $(HCFLAGS) -c -o $@ $^

# === Target Toolchain ====================================================
$(DEST)/%.srec: $(EDEST)/%.elf
	@$(ECHO) "\t(TARGET) OBJCOPY $@"
	@$(EOC) $(EOFLAGS) --output-target srec --srec-forceS3 $^ $@

$(EDEST)/%.elf: $(EDEST)/%.o
	@$(ECHO) "\t(TARGET) LINK\t$@"
	@$(ECC) -o $@ $^ $(ELFLAGS)

#$(EAPPS): $(EOBJS)
#	@$(ECHO) "\t(TARGET~)   LINK\t$@"
#	@$(ECC) -r $(EDEST)/%.o -o created.o

$(EDEST)/%.o: %.c
	@$(ECHO) "\t(TARGET) CC\t$@"
	@$(ECC) $(ECFLAGS) -c -o $@ $^
# =========================================================================

