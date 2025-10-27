R_DEFAULT := Rscript
DATASET ?= $(shell $(R_DEFAULT) -e 'cat(yaml::read_yaml("config.yml")$runtime$dataset)' 2>/dev/null)

.PHONY: dbscan h2o compare uncertainty all

dbscan:
	$(R_DEFAULT) scripts/02_dbscan_thin_degrees.R $(DATASET)

h2o:
	$(R_DEFAULT) scripts/03_h2o_train.R $(DATASET)

compare:
	$(R_DEFAULT) scripts/05_h20_vs_ssdm_results.R

uncertainty:
	$(R_DEFAULT) scripts/05_uncertainity.R

all: dbscan h2o compare uncertainty
