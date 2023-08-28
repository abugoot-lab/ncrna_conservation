.PHONY := build install

PKG_NAME := nocornacod
VERSION := $(shell awk '/^Version/ {print $$2}' $(PKG_NAME)/DESCRIPTION)

build: $(PKG_NAME)
	R -e $(subst PKG_NAME,$(PKG_NAME),'roxygen2::roxygenise("PKG_NAME")')
	R CMD build $<
	R CMD check $<_$(VERSION).tar.gz

install: $(PKG_NAME)_$(VERSION).tar.gz
	R CMD INSTALL $<
