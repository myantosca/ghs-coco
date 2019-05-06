LATEX=lualatex
BIBER=biber
PA2_REPORT_NAME=cosc6326-pa2-michael-yantosca
FINAL_REPORT_NAME=cosc6326-final-michael-yantosca
INCS=-I/usr/include/openmpi
LIBS=-lmpi
MPICPP=mpic++
CFLAGS=-g -std=c++11
BUILDDIR=./build
SRCDIR=./src

all: exe doc

exe: $(BUILDDIR) $(BUILDDIR)/txt2mpig $(BUILDDIR)/genmpig $(BUILDDIR)/bfs-coco $(BUILDDIR)/ghs-coco

.PHONY: clean clean-exe

$(BUILDDIR):
	@mkdir -p $(BUILDDIR)

doc: $(PA2_REPORT_NAME).pdf $(FINAL_REPORT_NAME).pdf

$(BUILDDIR)/%: $(BUILDDIR)/%.o
	@$(MPICPP) $(CFLAGS) $^ -o $@ $(LIBS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@$(MPICPP) $(CFLAGS) $(INCS) -c $(SRCDIR)/$*.cpp -o $@

%.pdf: %.tex
	@$(LATEX) -shell-escape $*.tex
	@$(BIBER) $*
	@$(LATEX) -shell-escape $*.tex

superclean: clean superclean-doc-$(FINAL_REPORT_NAME) superclean-doc-$(PA2_REPORT_NAME)

clean: clean-doc-$(FINAL_REPORT_NAME) clean-doc-$(PA2_REPORT_NAME) clean-exe
	@find -name '*~' | xargs rm -f

clean-exe:
	@rm -rf $(BUILDDIR)

superclean-doc-%:
	@rm -f $*.pdf

clean-doc-%:
	@rm -f $*.aux
	@rm -f $*.bbl
	@rm -f $*.bcf
	@rm -f $*.log
	@rm -f $*.run.xml
	@rm -f $*.dvi
	@rm -f $*.blg
	@rm -f $*.auxlock
	@rm -f $*.pyg
	@rm -f $*-figure*
	@rm -f $*.toc
	@rm -f $*.out
	@rm -f $*.snm
	@rm -f $*.nav
