# PROGRAMA
PROG = perfSL
OBJS = main.o mod1.o mod2.o mod3.o utils.o

# Compilador
CC     = gcc

# Diretório onde a biblioteca LIKWID está instalada
LIKWID_DIR = /home/soft/likwid

# Acrescentar onde apropriado as opções para incluir uso da biblioteca LIKWID
CFLAGS = -O0 -DLIKWID_PERFMON -I$(LIKWID_DIR)/include
LFLAGS = -L$(LIKWID_DIR)/lib -llikwid -lm

# Lista de arquivos para distribuição
DISTFILES = *.c *.h LEIAME* Makefile
DISTDIR = `basename ${PWD}`

.PHONY: all clean purge dist

%.o: %.c %.h
	$(CC) -c $(CFLAGS) -o $@ $<

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

clean:
	@echo "Limpando sujeira ..."
	@rm -f *~ *.bak

purge:  clean
	@echo "Limpando tudo ..."
	@rm -f $(PROG) $(OBJS) core a.out $(DISTDIR) $(DISTDIR).tar

dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tar) ..."
	@ln -s . $(DISTDIR)
	@tar -cvf $(DISTDIR).tar $(addprefix ./$(DISTDIR)/, $(DISTFILES))
	@rm -f $(DISTDIR)
