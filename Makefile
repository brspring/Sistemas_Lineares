PROG = gauss
OBJS = gauss.o utils.o

# Compilador
CC     = gcc

# Diretório onde a biblioteca LIKWID está instalada
LIKWID_DIR = /home/soft/likwid

# Acrescentar onde apropriado as opções para incluir uso da biblioteca LIKWID
CFLAGS = -O0 -DLIKWID_PERFMON -I$(LIKWID_DIR)/include
LFLAGS = -L$(LIKWID_LIB) -llikwid

# Lista de arquivos para distribuição
DISTFILES = gauss.c utils.c utils.h Makefile
DISTDIR = `basename ${PWD}`

.PHONY: all clean purge dist

%.o: %.c %.h
	$(CC) -c $(CFLAGS) -o $@ $<

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

clean:
	@echo "Limpando sujeira ..."
	@rm -f *~ *.bak $(OBJS)

purge:  clean
	@echo "Limpando tudo ..."
	@rm -f $(PROG) core a.out $(DISTDIR) $(DISTDIR).tar

dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tar) ..."
	@ln -s . $(DISTDIR)
	@tar -cvf $(DISTDIR).tar $(addprefix ./$(DISTDIR)/, $(DISTFILES))
	@rm -f $(DISTDIR)
