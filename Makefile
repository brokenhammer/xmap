# Compile the fortran modules used by python scripts

SUBDIRS = fortran

.PHONY: subdirs clean $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	$(MAKE) -C $(SUBDIRS) $@
