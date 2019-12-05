SUBMISSION_DIR= projet_bioinfo_auer_guisnet

all:

export:
	mkdir -p /tmp/$(SUBMISSION_DIR)
	cp README /tmp/$(SUBMISSION_DIR)
	cp AUTHORS /tmp/$(SUBMISSION_DIR)
	cp Makefile /tmp/$(SUBMISSION_DIR)
	cp -r src/ /tmp/$(SUBMISSION_DIR)
	tar cvjf $(SUBMISSION_DIR).tar.bz2 /tmp/$(SUBMISSION_DIR)
	rm -rf /tmp/$(SUBMISSION_DIR)

.PHONY: all export
